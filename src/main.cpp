#include "parallel.h"

#include <fstream>
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
#include "buildtime_info.h"
#include "restart.h"
#include "fld.h"



static void finalize() {
	FreeEuler();
	finalize_parallel();
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
    if (parameters::radiative_diffusion_enabled) {
		fld::init(data.get_n_radial(), data.get_n_azimuthal());
	}


    init_radialarrays();

    // Here planets are initialized feeling star potential
    data.get_planetary_system().init_system(options::parameter_file);
	quantities::state_disk_ecc_peri_calculation_center(data);
	data.get_massflow_tracker().init(data.get_planetary_system());
    init_binary_quadropole_moment(data.get_planetary_system());


    parameters::summarize_parameters();

    boundary_conditions::init_prescribed_time_variable_boundaries(data);
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

	boundary_conditions::copy_initial_values(data);

    if (start_mode::mode == start_mode::mode_restart) {
		restart_load(data);
    } else {
		// create 1D info files
		output::write_1D_info(data);
		
		MPI_Barrier(MPI_COMM_WORLD);
    }

    sim::PhysicalTimeInitial = sim::PhysicalTime;

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
