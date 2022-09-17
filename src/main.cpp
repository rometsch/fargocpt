#ifdef _OPENMP
#include <omp.h>
#else
#warning "Your compiler does not support OpenMP, at least with the flags you're using."
#endif

#include <experimental/filesystem>
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
#include "handle_signals.h"
#include "init.h"
#include "logging.h"
#include "options.h"
#include "output.h"
#include "parameters.h"
#include "quantities.h"
#include "selfgravity.h"
#include "split.h"
#include "start_mode.h"
#include "time.h"
#include "units.h"
#include "util.h"
#include "viscosity/viscosity.h"
#include "particles/particles.h"
#include "random/random.h"
#include "simulation.h"
#include "circumplanetary_mass.h"

// copy paste from https://github.com/jeffhammond/HPCInfo/blob/master/mpi/with-threads/mpi-openmp.c
#define MPI_THREAD_STRING(level)  \
		( level==MPI_THREAD_SERIALIZED ? "THREAD_SERIALIZED" : \
				( level==MPI_THREAD_MULTIPLE ? "THREAD_MULTIPLE" : \
						( level==MPI_THREAD_FUNNELED ? "THREAD_FUNNELED" : \
								( level==MPI_THREAD_SINGLE ? "THREAD_SINGLE" : "THIS_IS_IMPOSSIBLE" ) ) ) )

// TODO: NRad darf nicht größer sein als MAX1D
// copy operator for t_polargrid
// write all polargrids on error

static void init_parallel(int argc, char *argv[]) {
	int CPU_NameLength;
    char CPU_Name[MPI_MAX_PROCESSOR_NAME + 1];

	int provided;
	int requested = MPI_THREAD_FUNNELED;
	MPI_Init_thread(&argc, &argv, requested, &provided);
	if(provided < requested) {
		die("MPI_Init_thread provided %s when %s was requested.  Exiting. \n",
					   MPI_THREAD_STRING(provided), MPI_THREAD_STRING(requested));
	}

    // initialize MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &CPU_Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &CPU_Number);
    selfgravity::mpi_init();

#ifdef _OPENMP
	#pragma omp parallel
	{
		int omp_id  = omp_get_thread_num();
		int omp_num = omp_get_num_threads();
		printf("MPI rank # %2d OpenMP thread # %2d of %2d \n", CPU_Rank, omp_id, omp_num);
		fflush(stdout);
	}
#else
	printf("MPI rank # %2d \n", CPU_Rank);
	fflush(stdout);
#endif

    // are we master CPU?
    CPU_Master = (CPU_Rank == 0 ? 1 : 0);

	// print information on which processor we're running
    MPI_Get_processor_name(CPU_Name, &CPU_NameLength);
    CPU_Name[CPU_NameLength] = '\0';
    logging::print(LOG_INFO "fargo: running on %s\n", CPU_Name);
}

static void print_buildtimeinfo() {
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

}

static void finalize_parallel() {
	DeallocateBoundaryCommunicationBuffers();
    selfgravity::mpi_finalize();
    FreeSplitDomain();
	MPI_Finalize();
}

static void finalize() {
	FreeEuler();
	finalize_parallel();
}



static void init_damping(t_data &data) {
    if (parameters::is_damping_initial) {
	// save starting values (needed for damping)
		copy_polargrid(data[t_data::V_RADIAL0], data[t_data::V_RADIAL]);
		copy_polargrid(data[t_data::V_AZIMUTHAL0], data[t_data::V_AZIMUTHAL]);
		copy_polargrid(data[t_data::SIGMA0], data[t_data::SIGMA]);
		copy_polargrid(data[t_data::ENERGY0], data[t_data::ENERGY]);
    }
}

static void restart_load(t_data &data) {

	sim::N_monitor = 0;
	start_mode::restart_from = output::load_misc();

	if (parameters::is_damping_initial) {
	    // load grids at t = 0
	    const std::string snapshot_dir_old = output::snapshot_dir;
	    output::snapshot_dir = output::outdir + "snapshots/damping";
	    if (!std::experimental::filesystem::exists(output::snapshot_dir)) {
		logging::print_master(
		    LOG_ERROR
		    "Damping zone activated but no snapshot with damping data found. Make sure to copy the 'damping' snapshot!\n");
		PersonalExit(1);
	    }

	    logging::print_master(LOG_INFO
				  "Loading polargrinds for damping...\n");
	    data[t_data::SIGMA].read2D();
	    data[t_data::V_RADIAL].read2D();
	    data[t_data::V_AZIMUTHAL].read2D();
	    if (parameters::Adiabatic) {
		data[t_data::ENERGY].read2D();
	    }
	    output::snapshot_dir = snapshot_dir_old;

	    // save starting values (needed for damping)
	    copy_polargrid(data[t_data::V_RADIAL0], data[t_data::V_RADIAL]);
	    copy_polargrid(data[t_data::V_AZIMUTHAL0],
			   data[t_data::V_AZIMUTHAL]);
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
		compute_heating_cooling_for_CFL(data, sim::PhysicalTime);
	    }
	    if (data[t_data::QMINUS].file_exists()) {
		data[t_data::QMINUS].read2D();
	    } else {
		logging::print_master(
		    LOG_INFO
		    "Cannot read Qminus, no bitwise identical restarting possible!\n");
		compute_heating_cooling_for_CFL(data, sim::PhysicalTime);
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

	recalculate_derived_disk_quantities(data, sim::PhysicalTime);

	if (parameters::variableGamma) {

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

		compute_temperature(data);
		compute_sound_speed(data, sim::PhysicalTime);
		compute_scale_height(data, sim::PhysicalTime);
		compute_pressure(data);
	    viscosity::update_viscosity(data);
	}

}

int main(int argc, char *argv[])
{

    register_signal_handlers();

    t_data data;

	// TODO: discuss why this needs to be done explicitly
    resize_radialarrays(MAX1D);

	init_parallel(argc, argv);
    
	print_buildtimeinfo();

    // control behavior for floating point exceptions trapping (default is not
    // to do anything)
    // setfpe();

    // handle command line parameters
    options::parse(argc, argv);

    ReadVariables(options::parameter_file, data, argc, argv);

    // check if there is enough free space for all outputs (check before any
    // output are files created)
    output::check_free_space(data);

	if (options::memory_usage) {
		data.print_memory_usage(NRadial, NAzimuthal);
		PersonalExit(0);
    }

    parameters::write_grid_data_to_file();

    units::print_code_units();
    units::write_code_unit_file();
    constants::print_constants();
    output::write_output_version();

    TellEverything();

    SplitDomain();

    if (options::disable)
	PersonalExit(0);

    data.set_size(GlobalNRadial, NAzimuthal, NRadial, NAzimuthal);

	fargo_random::init();

    init_radialarrays();

    // Here planets are initialized feeling star potential
    data.get_planetary_system().init_system(options::parameter_file);
	quantities::state_disk_ecc_peri_calculation_center(data);

    parameters::summarize_parameters();

    boundary_conditions::init_prescribed_time_variable_boundaries(data);
    init_physics(data);
	sim::CalculateTimeStep(data);

    // update planet velocity due to disk potential
    if (parameters::disk_feedback) {
	ComputeDiskOnNbodyAccel(data);
	data.get_planetary_system().correct_velocity_for_disk_accel();
    }

    if (parameters::integrate_particles) {
	particles::init(data);
    }

	init_damping(data);

    if (start_mode::mode == start_mode::mode_restart) {
		restart_load(data);
    } else {
		// create planet files
		data.get_planetary_system().create_planet_files();

		// create 1D info files
		output::write_1D_info(data);
		
		MPI_Barrier(MPI_COMM_WORLD);
    }

    sim::PhysicalTimeInitial = sim::PhysicalTime;

    logging::start_timer();

    CommunicateBoundaries(&data[t_data::SIGMA], &data[t_data::V_RADIAL],
			  &data[t_data::V_AZIMUTHAL], &data[t_data::ENERGY]);

    CommunicateBoundaries(&data[t_data::SIGMA0], &data[t_data::V_RADIAL0],
			  &data[t_data::V_AZIMUTHAL0], &data[t_data::ENERGY0]);


	if (start_mode::mode != start_mode::mode_restart) {
		sim::handle_outputs(data);
	}

	sim::run(data);	
    

	finalize();

    return 0;
}
