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
#include "viscosity.h"
#include "particles/particles.h"
#include "random/random.h"
#include "simulation.h"
#include "circumplanetary_mass.h"

// TODO: NRad darf nicht größer sein als MAX1D
// copy operator for t_polargrid
// write all polargrids on error

int main(int argc, char *argv[])
{

    register_signal_handlers();

    t_data data;

    sim::N_hydro_iter = 0;
    sim::N_monitor = 0;

    resize_radialarrays(MAX1D);


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

    // control behavior for floating point exceptions trapping (default is not
    // to do anything)
    // setfpe();

    // handle command line parameters
    options::parse(argc, argv);

    ReadVariables(options::parameter_file, data, argc, argv);

    // check if there is enough free space for all outputs (check before any
    // output are files created)
    output::check_free_space(data);

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

	fargo_random::init();

    init_radialarrays();

    // Here planets are initialized feeling star potential
    data.get_planetary_system().init_system(options::parameter_file);

    // data.get_planetary_system().correct_planet_accretion();

	// TODO: check with Lucas why needed?
	if(parameters::heating_star_enabled){
		data[t_data::ASPECTRATIO].set_do_before_write(nullptr);
	}

    parameters::VISCOUS_ACCRETION = false;
    if (parameters::boundary_inner ==
	parameters::boundary_condition_viscous_outflow) {
	parameters::VISCOUS_ACCRETION = true;
    }
    for (unsigned int k = 0;
	 k < data.get_planetary_system().get_number_of_planets(); ++k) {
	if (data.get_planetary_system().get_planet(k).get_acc() < 0.0) {
	    parameters::VISCOUS_ACCRETION = true;
	}
    }

    logging::print_master(LOG_INFO "planets loaded.\n");

    if (parameters::VISCOUS_ACCRETION) {
	logging::print_master(
	    LOG_INFO
	    "VISCOUS_ACCRETION is true, recomputing viscosity before accreting mass.\n");
    }

    if ((data.get_planetary_system().get_number_of_planets() <= 1) &&
	(parameters::corotating)) {
	logging::print_master(
	    LOG_ERROR
	    "Error: Corotating frame is not possible with 0 or 1 planets.\n");
	PersonalExit(1);
    }

    parameters::summarize_parameters();
	if (config::cfg.get_flag("WriteDefaultValues", "no")) {
		config::cfg.write_default(output::outdir + "default_config.yml");
	}


    boundary_conditions::init_prescribed_time_variable_boundaries(data);
    init_physics(data);
	sim::init(data);

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

    if (start_mode::mode == start_mode::mode_restart) {

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

	    compute_temperature(data, true);
	    compute_sound_speed(data, true);
	    compute_scale_height(data, true);
	    compute_pressure(data, true);
	    viscosity::update_viscosity(data);
	}

    } else {
	// create planet files
	sim::dtemp = 0.0;
	data.get_planetary_system().create_planet_files();

	// create 1D info files
	if (CPU_Master) {
	    output::write_1D_info(data);
	}
	MPI_Barrier(MPI_COMM_WORLD);
    }
    sim::PhysicalTimeInitial = sim::PhysicalTime;

    logging::start_timer();

    CommunicateBoundaries(&data[t_data::SIGMA], &data[t_data::V_RADIAL],
			  &data[t_data::V_AZIMUTHAL], &data[t_data::ENERGY]);

    CommunicateBoundaries(&data[t_data::SIGMA0], &data[t_data::V_RADIAL0],
			  &data[t_data::V_AZIMUTHAL0], &data[t_data::ENERGY0]);


	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	sim::run(data);	
    

    logging::print_runtime_final();

    // free up everything
    DeallocateBoundaryCommunicationBuffers();
    FreeEuler();

    selfgravity::mpi_finalize();
    FreeSplitDomain();

    MPI_Finalize();

    return 0;
}
