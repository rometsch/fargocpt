#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string.h>

#include "Force.h"
#include "Interpret.h"
#include "LowTasks.h"
#include "Pframeforce.h"
#include "SideEuler.h"
#include "SourceEuler.h"
#include "Stockholm.h"
#include "Theo.h"
#include "boundary_conditions.h"
#include "commbound.h"
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
#include "config.h"

int TimeToWrite;
int dimfxy = 11, Restart = 0;
static int StillWriteOneOutput;
extern int Corotating;
extern int SelfGravity, SGZeroMode;

unsigned int nTimeStep;

// TODO: NRad darf nicht größer sein als MAX1D
// copy operator for t_polargrid
// Handle SIGTERM
// write all polargrids on error

int main(int argc, char *argv[])
{	
    t_data data;

    N_iter = 0;

    resize_radialarrays(MAX1D);

    unsigned int timeStepStart = 0;

    int CPU_NameLength;
    Force *force;
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
	force = AllocateForce(dimfxy);

	// Here planets are initialized feeling star potential
	data.get_planetary_system().read_from_file(PLANETCONFIG);
	logging::print_master(LOG_INFO "planets loaded.\n");


    if ((data.get_planetary_system().get_number_of_planets() <= 1) &&
	(Corotating == YES)) {
	logging::print_master(
	    LOG_ERROR
	    "Error: Corotating frame is not possible with 0 or 1 planets.\n");
	PersonalExit(1);
    }

    init_physics(data);

	// update planet velocity due to disk potential
	if (parameters::disk_feedback) {
		ComputeDiskOnNbodyAccel(force, data);
	}
	data.get_planetary_system().correct_velocity_for_disk_accel();
	logging::print_master(LOG_INFO "planets initialised.\n");

    if (parameters::integrate_particles) {
	particles::init(data);
    }


    // save starting values (needed for damping)
    copy_polargrid(data[t_data::V_RADIAL0], data[t_data::V_RADIAL]);
    copy_polargrid(data[t_data::V_AZIMUTHAL0], data[t_data::V_AZIMUTHAL]);
    copy_polargrid(data[t_data::DENSITY0], data[t_data::DENSITY]);
    copy_polargrid(data[t_data::ENERGY0], data[t_data::ENERGY]);

    // Initial Density is used to compute the circumplanetary mass with initial
    // density field
    mdcp0 = CircumPlanetaryMass(data);

    if (start_mode::mode == start_mode::mode_restart) {
	// TODO: fix for case that NINTERM changes (probably add small time step
	// to misc.dat)
	timeStepStart = start_mode::restart_from * NINTERM;
	logging::print_master(LOG_INFO "Restarting planetary system...\n");

	logging::print_master(LOG_INFO "Reading misc data...\n");
	PhysicalTime =
	    output::get_misc(start_mode::restart_from, "physical time");
	OmegaFrame = output::get_misc(start_mode::restart_from, "omega frame");
	FrameAngle = output::get_misc(start_mode::restart_from, "frame angle");

	// load grids at t = 0
	logging::print_master(LOG_INFO "Loading polargrinds at t = 0...\n");
	data[t_data::DENSITY].read2D((unsigned int)0);
	data[t_data::V_RADIAL].read2D((unsigned int)0);
	data[t_data::V_AZIMUTHAL].read2D((unsigned int)0);
	if (parameters::Adiabatic)
	    data[t_data::ENERGY].read2D((unsigned int)0);

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
	data[t_data::DENSITY].read2D(start_mode::restart_from);
	data[t_data::V_RADIAL].read2D(start_mode::restart_from);
	data[t_data::V_AZIMUTHAL].read2D(start_mode::restart_from);
	if (parameters::Adiabatic)
	    data[t_data::ENERGY].read2D(start_mode::restart_from);

	if (parameters::integrate_particles)
	    particles::restart(start_mode::restart_from);

	// restart planetary system
	data.get_planetary_system().restart(start_mode::restart_from);


	CommunicateBoundaries(&data[t_data::DENSITY], &data[t_data::V_RADIAL],
			      &data[t_data::V_AZIMUTHAL],
			      &data[t_data::ENERGY]);
	boundary_conditions::apply_boundary_condition(data, 0.0, false);
	recalculate_derived_disk_quantities(data, true);
    } else {
	// create planet files
	data.get_planetary_system().create_planet_files();

	// create 1D info files
	if (CPU_Master) {
	    output::write_1D_info(data);

	    // create mass flow info file
	    if (parameters::write_massflow) {
		output::write_massflow_info(data);
	    }
	}
	MPI_Barrier(MPI_COMM_WORLD);
    }
    PhysicalTimeInitial = PhysicalTime;

    logging::start_timer();

    for (nTimeStep = timeStepStart; nTimeStep <= NTOT; ++nTimeStep) {
	data.get_planetary_system().calculate_orbital_elements();
	// write outputs

	bool force_update_for_output = true;
	TimeStep = (nTimeStep / NINTERM); // note: integer division
	bool write_complete_output = NINTERM * TimeStep == nTimeStep;
	if (write_complete_output) {
	    // Outputs are done here
	    TimeToWrite = YES;
	    force_update_for_output = false;

	    // write polar grids
	    output::write_grids(data, TimeStep, N_iter, PhysicalTime);
	    // write planet data
	    data.get_planetary_system().write_planets(TimeStep, false);
	    // write misc stuff (important for resuming)
	    output::write_misc(TimeStep);
	    // write time info for coarse output
	    output::write_coarse_time(TimeStep, nTimeStep);
	    // write particles
	    if (parameters::integrate_particles)
		particles::write(TimeStep);
	    if ((OnlyInit) || ((GotoNextOutput) && (!StillWriteOneOutput))) {
		PersonalExit(0);
	    }
	    StillWriteOneOutput--;
	} else {
	    TimeToWrite = NO;
	}

	//(void) InnerOutputCounter;
	// InnerOutputCounter++;
	// if (InnerOutputCounter == 1) {
	if ((write_complete_output || parameters::write_at_every_timestep)) {
	    // InnerOutputCounter = 0;
	    data.get_planetary_system().write_planets(TimeStep, true);
	    // WriteBigPlanetSystemFile(sys, TimeStep);
	    UpdateLog(data, force, TimeStep, PhysicalTime);
	    if (Stockholm)
		UpdateLogStockholm(data, PhysicalTime);
	}

	// write disk quantities like eccentricity, ...
	if (write_complete_output || parameters::write_at_every_timestep) {
	    output::write_quantities(data, TimeStep, nTimeStep,
				     force_update_for_output);
	}
	if (write_complete_output || parameters::write_torques) {
	    output::write_torques(data, TimeStep, force_update_for_output);
	}
	if (parameters::write_lightcurves &&
	    (parameters::write_at_every_timestep || write_complete_output)) {
	    output::write_lightcurves(data, TimeStep, force_update_for_output);
	}
	if (parameters::write_massflow && nTimeStep != timeStepStart) {
	    output::write_massflow(data, TimeStep);
	}

	// Exit if last timestep reached and last output is written
	if (nTimeStep == NTOT) {
	    break;
	}

	// do hydro and nbody
	AlgoGas(nTimeStep, force, data);
    }

    logging::print_runtime_final();

	// free up everything
	config::free_config_list();
	DeallocateBoundaryCommunicationBuffers();
	free(OUTPUTDIR);
	free(PLANETCONFIG);
	delete[] options::parameter_file;
    FreeEuler();
    FreeForce(force);

    selfgravity::mpi_finalize();
    FreeSplitDomain();

    MPI_Finalize();

    return 0;
}
