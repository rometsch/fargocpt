#include "restart.h"
#include "parallel.h"
#include "viscosity/viscosity.h"
#include "output.h"
#include "logging.h"
#include <string>
#include "simulation.h"
#include "parameters.h"
#include "start_mode.h"
#include "boundary_conditions.h"
#include "Theo.h"
#include <experimental/filesystem>
#include "SourceEuler.h"
#include "particles/particles.h"
#include "circumplanetary_mass.h"


void restart_load(t_data &data) {

	sim::N_monitor = 0;
	start_mode::restart_from = output::load_misc();

	if (boundary_conditions::initial_values_needed()) {
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
		boundary_conditions::copy_initial_values(data);
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

	if(parameters::massoverflow){
	data.get_massflow_tracker().read_from_file();
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