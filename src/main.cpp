#include <fstream>
#include <mpi.h>
#include <string.h>
#include <experimental/filesystem>

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
#include "handle_signals.h"
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
	output::write_output_version();

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

	if (parameters::is_damping_initial) {
		// save starting values (needed for damping)
		copy_polargrid(data[t_data::V_RADIAL0], data[t_data::V_RADIAL]);
		copy_polargrid(data[t_data::V_AZIMUTHAL0], data[t_data::V_AZIMUTHAL]);
		copy_polargrid(data[t_data::SIGMA0], data[t_data::SIGMA]);
		copy_polargrid(data[t_data::ENERGY0], data[t_data::ENERGY]);
	}

    bool dont_do_restart_output_at_start = false;
    if (start_mode::mode == start_mode::mode_restart) {
	dont_do_restart_output_at_start = true;
	
	N_outer_loop = 0;
	start_mode::restart_from = output::load_misc();

	if (parameters::is_damping_initial) {
		// load grids at t = 0
		const std::string snapshot_dir_old = snapshot_dir;
		snapshot_dir = std::string(OUTPUTDIR) + "/snapshots/damping";
		if (!std::experimental::filesystem::exists(snapshot_dir)) {
			logging::print_master(LOG_ERROR "Damping zone activated but no snapshot with damping data found. Make sure to copy the 'damping' snapshot!\n");
			PersonalExit(1);
		}

		logging::print_master(LOG_INFO "Loading polargrinds for damping...\n");
		data[t_data::SIGMA].read2D();
		data[t_data::V_RADIAL].read2D();
		data[t_data::V_AZIMUTHAL].read2D();
		if (parameters::Adiabatic) {
			data[t_data::ENERGY].read2D();
		}
		snapshot_dir = snapshot_dir_old;

		// save starting values (needed for damping)
		copy_polargrid(data[t_data::V_RADIAL0], data[t_data::V_RADIAL]);
		copy_polargrid(data[t_data::V_AZIMUTHAL0], data[t_data::V_AZIMUTHAL]);
		copy_polargrid(data[t_data::SIGMA0], data[t_data::SIGMA]);
		copy_polargrid(data[t_data::ENERGY0], data[t_data::ENERGY]);
	}

	// recalculate SigmaMed/EnergyMed
	RefillSigma(&data[t_data::SIGMA]);
	if (parameters::Adiabatic)
	    RefillEnergy(&data[t_data::ENERGY]);

	// load grids at t = restart_from
	logging::print_master(LOG_INFO "Loading polargrinds at t = %u...\n",
			      start_mode::restart_from);
	data[t_data::SIGMA].read2D();
	data[t_data::V_RADIAL].read2D();
	data[t_data::V_AZIMUTHAL].read2D();
	if (parameters::Adiabatic) {
	    data[t_data::ENERGY].read2D();

		if (data[t_data::QPLUS].file_exists()) {
		data[t_data::QPLUS].read2D();
	    } else {
		logging::print_master(
		    LOG_INFO
		    "Cannot read Qplus, no bitwise identical restarting possible!\n");
		compute_heating_cooling_for_CFL(data);
	    }
		if (data[t_data::QMINUS].file_exists()) {
		data[t_data::QMINUS].read2D();
	    } else {
		logging::print_master(
		    LOG_INFO
		    "Cannot read Qminus, no bitwise identical restarting possible!\n");
		compute_heating_cooling_for_CFL(data);
	    }
	}
	if (parameters::integrate_particles) {
	    particles::restart();
	}

	// restart planetary system
	logging::print_master(LOG_INFO "Restarting planetary system...\n");
	data.get_planetary_system().restart();

	ComputeCircumPlanetaryMasses(data);

	MPI_Barrier(MPI_COMM_WORLD);
	logging::print_master(LOG_INFO
			      "Finished restarting planetary system.\n");

	recalculate_derived_disk_quantities(data, true);

	if(parameters::variableGamma){

		// For bitwise exact restarting with PVTE
		if (data[t_data::GAMMAEFF].file_exists()) {
		data[t_data::GAMMAEFF].read2D();
		}
		if (data[t_data::MU].file_exists()) {
		data[t_data::MU].read2D();
		}

		if (data[t_data::GAMMA1].file_exists()) {
		data[t_data::GAMMA1].read2D();
		}

		compute_temperature(data, true);
		compute_sound_speed(data, true);
		compute_scale_height(data, true);
		compute_pressure(data, true);
		viscosity::update_viscosity(data);

	}

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

	CommunicateBoundaries(&data[t_data::SIGMA], &data[t_data::V_RADIAL],
			  &data[t_data::V_AZIMUTHAL], &data[t_data::ENERGY]);

	CommunicateBoundaries(&data[t_data::SIGMA0], &data[t_data::V_RADIAL0],
			  &data[t_data::V_AZIMUTHAL0], &data[t_data::ENERGY0]);

    for (; N_outer_loop <= NTOT; ++N_outer_loop) {
	logging::print_master(LOG_INFO "Start of iteration %u of %u\n", N_outer_loop, NTOT);
	// write outputs

	bool force_update_for_output = true;
	N_output = (N_outer_loop / NINTERM); // note: integer division
	bool write_complete_output = (NINTERM * N_output == N_outer_loop);
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
		last_snapshot_dir = snapshot_dir;
		output::write_full_output(data, std::to_string(N_output));
		output::cleanup_autosave();

		if (N_output == 0 && parameters::damping) {
			// Write damping data as a reference.
			const std::string snapshot_dir_old = snapshot_dir;
			output::write_full_output(data, "damping", false);
			snapshot_dir = snapshot_dir_old;
		}

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
		ComputeCircumPlanetaryMasses(data);
	    data.get_planetary_system().write_planets(1);
	    // WriteBigPlanetSystemFile(sys, TimeStep);
	}

	// write disk quantities like eccentricity, ...
	if ((write_complete_output || parameters::write_at_every_timestep) &&
	    !(dont_do_restart_output_at_start) &&
	    parameters::write_disk_quantities) {
	    output::write_quantities(data, force_update_for_output);
	}

	if (write_complete_output && parameters::write_torques) {
		output::write_torques(data, force_update_for_output);
	}
	if (parameters::write_lightcurves &&
	    (parameters::write_at_every_timestep || write_complete_output) &&
	    !(dont_do_restart_output_at_start)) {
	    output::write_lightcurves(data, N_output, force_update_for_output);
	}
	dont_do_restart_output_at_start = false;

	// Exit if last timestep reached and last output is written
	if (N_outer_loop == NTOT) {
		logging::print_master(LOG_INFO "Reached end of simulation at iteration %u of %u\n", N_outer_loop, NTOT);
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
