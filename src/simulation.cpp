#include "simulation.h"
#include "logging.h"
#include "parameters.h"
#include "output.h"
#include "circumplanetary_mass.h"
#include "Pframeforce.h"
#include "start_mode.h"
#include "SourceEuler.h"


static void write_snapshot(t_data &data) {
	// Outputs are done here
	output::last_snapshot_dir = output::snapshot_dir;
	output::write_full_output(data, std::to_string(N_output));
	output::cleanup_autosave();

	if (N_output == 0 && parameters::damping) {
	// Write damping data as a reference.
	const std::string snapshot_dir_old = output::snapshot_dir;
	output::write_full_output(data, "damping", false);
	output::snapshot_dir = snapshot_dir_old;
	}
}

static void handle_outputs(t_data &data) {
	bool need_update_for_output = true;
	N_output = (N_outer_loop / parameters::NINTERM); // note: integer division
	bool write_complete_output = (parameters::NINTERM * N_output == N_outer_loop);

	/// asure planet torques are computed
	if (!parameters::disk_feedback &&
	    (write_complete_output || parameters::write_at_every_timestep)) {
	    ComputeDiskOnNbodyAccel(data);
	}

	if (write_complete_output) {
	    need_update_for_output = false;
		write_snapshot(data);
	}

	if (write_complete_output || parameters::write_at_every_timestep) {
	    ComputeCircumPlanetaryMasses(data);
	    data.get_planetary_system().write_planets(1);
	}

	// write disk quantities like eccentricity, ...
	if ((write_complete_output || parameters::write_at_every_timestep) &&
	    parameters::write_disk_quantities) {
	    output::write_quantities(data, need_update_for_output);
	}

	if (write_complete_output && parameters::write_torques) {
	    output::write_torques(data, need_update_for_output);
	}
	if (parameters::write_lightcurves &&
	    (parameters::write_at_every_timestep || write_complete_output)) {
	    output::write_lightcurves(data, N_output, need_update_for_output);
	}

}

void run_simulation(t_data &data) {

	if (start_mode::mode != start_mode::mode_restart) {
		handle_outputs(data);
	}

    for (; N_outer_loop <= parameters::NTOT; N_outer_loop++) {
		
		logging::print_master(LOG_INFO "Start of iteration %u of %u\n",
		N_outer_loop, parameters::NTOT);

		// do hydro and nbody
		AlgoGas(data);
		dtemp = 0.0;

		handle_outputs(data);
    }
	logging::print_master(
			LOG_INFO "Reached end of simulation at iteration %u of %u\n",
			N_outer_loop, parameters::NTOT);
}