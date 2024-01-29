#include "parallel.h"

#include <fstream>
#include <string.h>
#include <filesystem>

#include "Interpret.h"
#include "LowTasks.h"
#include "Pframeforce.h"
#include "SourceEuler.h"
#include "Theo.h"
#include "boundary_conditions/boundary_conditions.h"
#include "commbound.h"
#include "config.h"
#include "constants.h"
#include "data.h"
#include "global.h"
#include "handle_signals.h"
#include "init.h"
#include "logging.h"
#include "options.h"
#include "output.h"
#include "parameters.h"
#include "split.h"
#include "start_mode.h"
#include "units.h"
#include "particles/particles.h"
#include "random/random.h"
#include "simulation.h"
#include "buildtime_info.h"
#include "restart.h"
#include "fld.h"



static void finalize() {
	FreeEuler();
	finalize_parallel();
    boundary_conditions::cleanup_custom();

    if (CPU_Master && options::pidfile != "" && std::filesystem::exists(options::pidfile)) {
        std::filesystem::remove(options::pidfile);
    }
    logging::finalize();
}


int main(int argc, char *argv[])
{
    // handle command line parameters
    options::parse(argc, argv);

    register_signal_handlers();

    t_data data;
    
	init_parallel(argc, argv);

	// TODO: discuss why this needs to be done explicitly
    resize_radialarrays(MAX1D);

    
	print_buildtimeinfo();

    // control behavior for floating point exceptions trapping (default is not
    // to do anything)
    // setfpe();


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
    units::write_code_units_file();
    constants::print_constants();
    constants::write_code_constants_file();
    output::write_output_version();

    TellEverything();

    SplitDomain();

    if (options::disable)
	PersonalExit(0);

    data.set_size(GlobalNRadial, NAzimuthal, NRadial, NAzimuthal);

	fargo_random::init();
    if (fld::radiative_diffusion_enabled) {
		fld::init(data.get_n_radial(), data.get_n_azimuthal());
	}

    init_radialarrays();

    // Here planets are initialized feeling star potential
    data.get_planetary_system().init_system();
	data.get_massflow_tracker().init(data.get_planetary_system());
    init_binary_quadropole_moment(data.get_planetary_system());

    parameters::summarize_parameters();
    if (CPU_Master) {
        config::cfg.exit_on_unknown_key();
    }

    boundary_conditions::init(data);
    init_physics(data);
	sim::CalculateTimeStep(data);

	if(parameters::planet_orbit_disk_test){
		data.get_planetary_system().get_planet(0).set_mass(1.0e-66);
	} else {
    // update planet velocity due to disk potential
    if (parameters::disk_feedback) {
	ComputeDiskOnNbodyAccel(data);
	data.get_planetary_system().correct_velocity_for_disk_accel();
    }
	}

    if (parameters::integrate_particles) {
	particles::init(data);
    }

    if (start_mode::mode == start_mode::mode_restart) {
		restart_load(data);
    } else {
		// create 1D info files
		output::write_1D_info(data);
        output::write_2D_info(data);
		
		MPI_Barrier(MPI_COMM_WORLD);
    }

    sim::timeInitial = sim::time;

    logging::start_timer();

    CommunicateBoundariesAll(data);
    CommunicateBoundariesAllInitial(data);

	if (start_mode::mode != start_mode::mode_restart) {
		sim::handle_outputs(data);
	}

	sim::run(data);	
    

	finalize();

    return 0;
}
