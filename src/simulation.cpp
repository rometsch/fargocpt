#include "simulation.h"
#include "logging.h"
#include "parameters.h"
#include "output.h"
#include "circumplanetary_mass.h"
#include "Pframeforce.h"
#include "start_mode.h"
#include "SourceEuler.h"
#include "boundary_conditions.h"
#include "commbound.h"
#include "frame_of_reference.h"
#include "particles/particles.h"
#include "viscosity.h"
#include "TransportEuler.h"
#include "sts.h"
#include "accretion.h"
#include "cfl.h"

namespace sim {

double hydro_dt;
double last_dt;
double dtemp;

double PhysicalTime, PhysicalTimeInitial;
unsigned int N_output;
unsigned int N_outer_loop;
unsigned long int N_hydro_iter;


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

static double CalculateTimeStep(t_data &data, double dt)
{
	last_dt = dt;
	const double cfl_dt = cfl::condition_cfl(data);
	return std::min(parameters::CFL_max_var * last_dt, cfl_dt);;
}


/*
Do one step of the integration.

Returns the time covered.
*/
static void step(t_data &data, const double dt) {

	if (parameters::disk_feedback) {
	    ComputeDiskOnNbodyAccel(data);
	}

	refframe::ComputeIndirectTermDisk(data);
	if (parameters::disk_feedback) {
	    UpdatePlanetVelocitiesWithDiskForce(data, dt);
	}

	/// IndirectTerm is fully completed here (Disk + Nbody)
	refframe::ComputeIndirectTermNbody(data, dt);
	data.get_planetary_system().apply_indirect_term_on_Nbody(
		refframe::IndirectTerm, dt);


	if (parameters::integrate_planets) {
	    data.get_planetary_system().integrate(PhysicalTime, dt);
	    /// Nbody positions and velocities are not updated yet!
	}

	if (parameters::calculate_disk) {
	    /** Gravitational potential from star and planet(s) is computed and
	     * stored here*/
	    if (parameters::body_force_from_potential) {
		CalculateNbodyPotential(data);
	    } else {
		CalculateAccelOnGas(data);
	    }
	}

	if (parameters::integrate_particles) {
	    particles::integrate(data, dt);
	}

	/** Planets' positions and velocities are updated from gravitational
	 * interaction with star and other planets */
	if (parameters::integrate_planets) {
		data.get_planetary_system().copy_data_from_rebound();
		data.get_planetary_system().move_to_hydro_frame_center();

	    /// Needed for Aspectratio mode = 1
	    /// and to correctly compute circumplanetary disk mass
	    data.get_planetary_system().compute_dist_to_primary();

	    /// Needed if they can change and massoverflow or planet accretion
	    /// is on
	    data.get_planetary_system().calculate_orbital_elements();
	}

	/* Below we correct v_azimuthal, planet's position and velocities if we
	 * work in a frame non-centered on the star. Same for dust particles. */
	refframe::handle_corotation(data, dt);

	/* Now we update gas */
	if (parameters::calculate_disk) {
		//HandleCrash(data);

	    update_with_sourceterms(data, dt);

	    if (parameters::EXPLICIT_VISCOSITY) {
		// compute and add acceleartions due to disc viscosity as a
		// source term
		update_with_artificial_viscosity(data, dt);
		if (parameters::Adiabatic) {
		    SetTemperatureFloorCeilValues(data, __FILE__, __LINE__);
		}
		recalculate_derived_disk_quantities(data, true);
		ComputeViscousStressTensor(data);
		viscosity::update_velocities_with_viscosity(data, dt);
	    }

	    if (!parameters::EXPLICIT_VISCOSITY) {
		Sts(data, dt);
	    }

	    if (parameters::Adiabatic) {
		SubStep3(data, dt);
		if (parameters::radiative_diffusion_enabled) {
		    radiative_diffusion(data, dt);
		}
	    }

	    /// TODO moved apply boundaries here
		boundary_conditions::apply_boundary_condition(data, 0.0, false);

		Transport(data, &data[t_data::SIGMA], &data[t_data::V_RADIAL],
			  &data[t_data::V_AZIMUTHAL], &data[t_data::ENERGY],
			  dt);
	}

	PhysicalTime += dt;
	N_hydro_iter = N_hydro_iter + 1;
	logging::print_runtime_info(N_outer_loop / parameters::NINTERM, N_outer_loop,
				    dt);

	if (parameters::calculate_disk) {
	    CommunicateBoundaries(&data[t_data::SIGMA], &data[t_data::V_RADIAL],
				  &data[t_data::V_AZIMUTHAL],
				  &data[t_data::ENERGY]);

	    // We only recompute once, assuming that cells hit by planet
	    // accretion are not also hit by viscous accretion at inner
	    // boundary.
	    if (parameters::VISCOUS_ACCRETION) {
		compute_sound_speed(data, true);
		compute_scale_height(data, true);
		viscosity::update_viscosity(data);
	    }

		// minimum density is assured inside AccreteOntoPlanets
	    accretion::AccreteOntoPlanets(data, dt);

	    boundary_conditions::apply_boundary_condition(data, dt, true);

	    // const double total_disk_mass_new =
	    //  quantities::gas_total_mass(data, 2.0*RMAX);

	    // data[t_data::DENSITY] *=
	    //(total_disk_mass_old / total_disk_mass_new);

	    CalculateMonitorQuantitiesAfterHydroStep(data, N_outer_loop,
						     dt);

	    if (parameters::variableGamma &&
		!parameters::VISCOUS_ACCRETION) { // If VISCOUS_ACCRETION is active,
				      // scale_height is already updated
		// Recompute scale height after Transport to update the 3D
		// density
		compute_sound_speed(data, true);
		compute_scale_height(data, true);
	    }
	    // this must be done after CommunicateBoundaries
	    recalculate_derived_disk_quantities(data, true);
	}
}

static void print_progress() {
	logging::print_master(
	LOG_VERBOSE
	"AlgoGas: Total: %*i/%i (%5.2f %%) - Timestep: %#7f/%#7f (%5.2f %%)\n",
	(int)ceil(log10(parameters::NTOT)), N_outer_loop, parameters::NTOT,
	(double)N_outer_loop / (double)parameters::NTOT * 100.0, dtemp, parameters::DT,
	dtemp / parameters::DT * 100.0);
}

void init(t_data &data) {
	boundary_conditions::apply_boundary_condition(data, 0.0, false); // TODO: move to init
	refframe::init_corotation(data); // TODO: move to init

	if (start_mode::mode != start_mode::mode_restart) {
		handle_outputs(data);
	}
}

void run(t_data &data) {

	init(data);

	// Jumpstart dt
	hydro_dt = last_dt;

	const double t_final = parameters::NTOT * parameters::DT;

	bool towrite = false;

    for (; PhysicalTime < t_final; N_hydro_iter++) {

		// TODO: check whether this is still necessary
		// do hydro and nbody
		// recalculate timestep, even for no_disk = true, so that particle drag has
    	// reasonable timestep size
		if (!parameters::calculate_disk) {
			hydro_dt = CalculateTimeStep(data, last_dt);
		}
					
		
		if (SIGTERM_RECEIVED) {
			output::write_full_output(data, "autosave");
			PersonalExit(0);
		}

		const double cfl_dt = CalculateTimeStep(data, hydro_dt);

		const double time_left_till_write = (N_outer_loop+1)*parameters::DT - PhysicalTime;
		if ( cfl_dt > time_left_till_write) {
			hydro_dt = time_left_till_write;
			N_outer_loop++;
			towrite = true;
		} else {
			hydro_dt = cfl_dt;
			towrite = false;
		}

		step(data, hydro_dt);

		if (towrite) {
			handle_outputs(data);
			print_progress();
		}
    }
	logging::print_master(
			LOG_INFO "Reached end of simulation at iteration %u of %u\n",
			N_outer_loop - 1, parameters::NTOT);
}


} // close namespace sim