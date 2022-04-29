#include <fstream>
#include <mpi.h>
#include <string.h>

#include "Force.h"
#include "Interpret.h"
#include "LowTasks.h"
#include "Pframeforce.h"
#include "SideEuler.h"
#include "SourceEuler.h"
#include "Theo.h"
#include "boundary_conditions.h"
#include "commbound.h"
#include "config.h"
#include "constants.h"
#include "data.h"
#include "fpe.h"
#include "global.h"
#include "init.h"
#include "logging.h"
#include "options.h"
#include "output.h"
#include "parameters.h"
#include "particles.h"
#include "quantities.h"
#include "selfgravity.h"
#include "split.h"
#include "start_mode.h"
#include "time.h"
#include "units.h"
#include "util.h"
#include "viscosity.h"
#include "handle_signals.h"


int TimeToWrite;
int Restart = 0;
static int StillWriteOneOutput;
extern int Corotating;
extern int SelfGravity, SGZeroMode;

// TODO: NRad darf nicht größer sein als MAX1D
// copy operator for t_polargrid
// Handle SIGTERM
// write all polargrids on error

int main(int argc, char *argv[])
{

	register_signal_handlers();

    t_data data;

    N_hydro_iter = 0;

    resize_radialarrays(MAX1D);

    N_outer_loop = 0;

    int CPU_NameLength;
    char CPU_Name[MPI_MAX_PROCESSOR_NAME + 1];

    // initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &CPU_Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &CPU_Number);
    selfgravity::mpi_init();

    // are we master CPU?
    CPU_Master = (CPU_Rank == 0 ? 1 : 0);

    // print some information about program
    logging::print_master(LOG_INFO "fargo: This file was compiled on %s, %s.\n",
			  __DATE__, __TIME__);
#ifdef GIT_COMMIT
    logging::print_master(LOG_INFO "fargo: Last git commit: %s\n", GIT_COMMIT);
#endif
#ifdef GIT_CHANGED
    logging::print_master(
	LOG_INFO "fargo: Files changed since git commit: %s\n", GIT_CHANGED);
#endif
#ifdef _GNU_SOURCE
    logging::print_master(LOG_INFO
			  "fargo: This version of FARGO used _GNU_SOURCE\n");
#endif
#ifdef NDEBUG
    logging::print_master(
	LOG_INFO
	"fargo: This version of FARGO used NDEBUG. So no assertion checks!\n");
#else
    logging::print_master(
	LOG_INFO
	"fargo: This version of FARGO does assertion checks! Compile with NDEBUG to speed up!\n");
#endif

    // print information on which processor we're running
    MPI_Get_processor_name(CPU_Name, &CPU_NameLength);
    CPU_Name[CPU_NameLength] = '\0';
    logging::print(LOG_INFO "fargo: running on %s\n", CPU_Name);

    // control behavoir for floating point exceptions trapping (default is not
    // to do anything)
    setfpe();

    // handle command line parameters
    options::parse(argc, argv);

    ReadVariables(options::parameter_file, data, argc, argv);

    // check if there is enough free space for all outputs (check before any
    // output are files created)
    output::check_free_space(data);

    parameters::summarize_parameters();
    parameters::write_grid_data_to_file();

    units::print_code_units();
    units::write_code_unit_file();
    constants::print_constants();

    SplitDomain();

    TellEverything();

    if (options::memory_usage) {
	data.print_memory_usage(NRadial, NAzimuthal);
	PersonalExit(0);
    }

    if (options::disable)
	PersonalExit(0);

    data.set_size(GlobalNRadial, NAzimuthal, NRadial, NAzimuthal);

    init_radialarrays();

    // Here planets are initialized feeling star potential
    data.get_planetary_system().read_from_file(PLANETCONFIG);
    data.get_planetary_system().compute_dist_to_primary();
    data.get_planetary_system().init_roche_radii();

    VISCOUS_ACCRETION = false;
    if (parameters::boundary_inner ==
	parameters::boundary_condition_viscous_outflow) {
	VISCOUS_ACCRETION = true;
    }
    for (unsigned int k = 0;
	 k < data.get_planetary_system().get_number_of_planets(); ++k) {
	if (data.get_planetary_system().get_planet(k).get_acc() < 0.0) {
	    VISCOUS_ACCRETION = true;
	}
    }

    logging::print_master(LOG_INFO "planets loaded.\n");

    if (VISCOUS_ACCRETION) {
	logging::print_master(
	    LOG_INFO
	    "VISCOUS_ACCRETION is true, recomputing viscosity before accreting mass.\n");
    }

    if ((data.get_planetary_system().get_number_of_planets() <= 1) &&
	(Corotating == YES)) {
	logging::print_master(
	    LOG_ERROR
	    "Error: Corotating frame is not possible with 0 or 1 planets.\n");
	PersonalExit(1);
    }

    boundary_conditions::init_prescribed_time_variable_boundaries(data);
    init_physics(data);

    // update planet velocity due to disk potential
    if (parameters::disk_feedback) {
	ComputeDiskOnNbodyAccel(data);
	data.get_planetary_system().correct_velocity_for_disk_accel();
    }
    logging::print_master(LOG_INFO "planets initialised.\n");

    if (parameters::integrate_particles) {
	particles::init(data);
    }

    // save starting values (needed for damping)
    copy_polargrid(data[t_data::V_RADIAL0], data[t_data::V_RADIAL]);
    copy_polargrid(data[t_data::V_AZIMUTHAL0], data[t_data::V_AZIMUTHAL]);
    copy_polargrid(data[t_data::DENSITY0], data[t_data::DENSITY]);
    copy_polargrid(data[t_data::ENERGY0], data[t_data::ENERGY]);

    bool dont_do_restart_output_at_start = false;
    if (start_mode::mode == start_mode::mode_restart) {
	dont_do_restart_output_at_start = true;

	start_mode::restart_from = output::get_misc(start_mode::restart_from,
						    start_mode::restart_debug);

	// load grids at t = 0
	logging::print_master(LOG_INFO "Loading polargrinds at t = 0...\n");
	data[t_data::DENSITY].read2D((unsigned int)0, false);
	data[t_data::V_RADIAL].read2D((unsigned int)0, false);
	data[t_data::V_AZIMUTHAL].read2D((unsigned int)0, false);
	if (parameters::Adiabatic) {
	    data[t_data::ENERGY].read2D((unsigned int)0, false);
	}

	// save starting values (needed for damping)
	copy_polargrid(data[t_data::V_RADIAL0], data[t_data::V_RADIAL]);
	copy_polargrid(data[t_data::V_AZIMUTHAL0], data[t_data::V_AZIMUTHAL]);
	copy_polargrid(data[t_data::DENSITY0], data[t_data::DENSITY]);
	copy_polargrid(data[t_data::ENERGY0], data[t_data::ENERGY]);

	// recalculate SigmaMed/EnergyMed
	RefillSigma(&data[t_data::DENSITY]);
	if (parameters::Adiabatic)
	    RefillEnergy(&data[t_data::ENERGY]);

	// load grids at t = restart_from
	logging::print_master(LOG_INFO "Loading polargrinds at t = %u...\n",
			      start_mode::restart_from);
	data[t_data::DENSITY].read2D(start_mode::restart_from,
				     start_mode::restart_debug);
	data[t_data::V_RADIAL].read2D(start_mode::restart_from,
				      start_mode::restart_debug);
	data[t_data::V_AZIMUTHAL].read2D(start_mode::restart_from,
					 start_mode::restart_debug);
	if (parameters::Adiabatic) {
	    data[t_data::ENERGY].read2D(start_mode::restart_from,
					start_mode::restart_debug);

	    if (data[t_data::QPLUS].get_write_2D() &&
		data[t_data::QPLUS].file_exists(start_mode::restart_from,
						start_mode::restart_debug)) {
		data[t_data::QPLUS].read2D(start_mode::restart_from,
					   start_mode::restart_debug);
	    } else {
		logging::print_master(
		    LOG_INFO
		    "Cannot read Qplus, no bitwise identical restarting possible!\n");
		compute_heating_cooling_for_CFL(data);
	    }
	    if (data[t_data::QMINUS].get_write_2D() &&
		data[t_data::QMINUS].file_exists(start_mode::restart_from,
						 start_mode::restart_debug)) {
		data[t_data::QMINUS].read2D(start_mode::restart_from,
					    start_mode::restart_debug);
	    } else {
		logging::print_master(
		    LOG_INFO
		    "Cannot read Qminus, no bitwise identical restarting possible!\n");
		compute_heating_cooling_for_CFL(data);
	    }
	}
	if (parameters::integrate_particles) {
	    if (start_mode::restart_debug) {
		die("Debug restart not implemented for particles yet!\n");
	    }
	    particles::restart(start_mode::restart_from);
	}

	// restart planetary system
	logging::print_master(LOG_INFO "Restarting planetary system...\n");
	data.get_planetary_system().restart(start_mode::restart_from,
					    start_mode::restart_debug);

	MPI_Barrier(MPI_COMM_WORLD);
	logging::print_master(LOG_INFO
			      "Finished restarting planetary system.\n");

	recalculate_derived_disk_quantities(data, true);
    } else {
	// create planet files
	dtemp = 0.0;
	data.get_planetary_system().create_planet_files();

	// create 1D info files
	if (CPU_Master) {
	    output::write_1D_info(data);
	}
	MPI_Barrier(MPI_COMM_WORLD);
    }
    PhysicalTimeInitial = PhysicalTime;

    logging::start_timer();

    CommunicateBoundaries(&data[t_data::DENSITY], &data[t_data::V_RADIAL],
			  &data[t_data::V_AZIMUTHAL], &data[t_data::ENERGY]);

    CommunicateBoundaries(&data[t_data::DENSITY0], &data[t_data::V_RADIAL0],
			  &data[t_data::V_AZIMUTHAL0], &data[t_data::ENERGY0]);

    for (; N_outer_loop <= NTOT; ++N_outer_loop) {
	data.get_planetary_system().compute_dist_to_primary();
	data.get_planetary_system().calculate_orbital_elements();
	ComputeCircumPlanetaryMasses(data);
	// write outputs

	bool force_update_for_output = true;
	N_output = (N_outer_loop / NINTERM); // note: integer division
	bool write_complete_output = NINTERM * N_output == N_outer_loop;
	if (dont_do_restart_output_at_start) {
	    write_complete_output = false;
	}
	/// asure planet torques are computed
	if (!parameters::disk_feedback &&
	    (write_complete_output || parameters::write_at_every_timestep)) {
	    ComputeDiskOnNbodyAccel(data);
	}

	if (write_complete_output) {
	    // Outputs are done here
	    TimeToWrite = YES;
	    force_update_for_output = false;

	    // write polar grids
	    output::write_grids(data, N_output, N_hydro_iter, PhysicalTime,
				false);
	    // write planet data
	    data.get_planetary_system().write_planets(N_output, 0);
	    // write misc stuff (important for resuming)
	    output::write_misc(false);
	    // write time info for coarse output
	    output::write_coarse_time(N_output, N_outer_loop);
	    // write particles
	    if (parameters::integrate_particles)
		particles::write(N_output);
	    if (GotoNextOutput && (!StillWriteOneOutput)) {
		PersonalExit(0);
	    }
	    StillWriteOneOutput--;
	} else {
	    TimeToWrite = NO;
	}

	//(void) InnerOutputCounter;
	// InnerOutputCounter++;
	// if (InnerOutputCounter == 1) {
	if ((write_complete_output || parameters::write_at_every_timestep) &&
	    !(dont_do_restart_output_at_start)) {
	    // InnerOutputCounter = 0;
	    data.get_planetary_system().write_planets(N_output, 1);
	    // WriteBigPlanetSystemFile(sys, TimeStep);
	}

	// write disk quantities like eccentricity, ...
	if ((write_complete_output || parameters::write_at_every_timestep) &&
	    !(dont_do_restart_output_at_start) &&
	    parameters::write_disk_quantities) {
	    output::write_quantities(data, force_update_for_output);
	}

	if (write_complete_output && parameters::write_torques) {
	    output::write_torques(data, N_output, force_update_for_output);
	}
	if (parameters::write_lightcurves &&
	    (parameters::write_at_every_timestep || write_complete_output) &&
	    !(dont_do_restart_output_at_start)) {
	    output::write_lightcurves(data, N_output, force_update_for_output);
	}
	dont_do_restart_output_at_start = false;

	// Exit if last timestep reached and last output is written
	if (N_outer_loop == NTOT) {
	    break;
	}

	// do hydro and nbody
	AlgoGas(data);
	dtemp = 0.0;
    }

    logging::print_runtime_final();

    // free up everything
    config::free_config_list();
    DeallocateBoundaryCommunicationBuffers();
    free(OUTPUTDIR);
    free(PLANETCONFIG);
    free(PRESCRIBED_BOUNDARY_OUTER_FILE);
    delete[] options::parameter_file;
    FreeEuler();

    selfgravity::mpi_finalize();
    FreeSplitDomain();

    MPI_Finalize();

    return 0;
}
