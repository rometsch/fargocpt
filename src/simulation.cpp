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
#include "viscosity/viscosity.h"
#include "viscosity/artificial_viscosity.h"
#include "TransportEuler.h"
#include "sts.h"
#include "accretion.h"
#include "cfl.h"
#include "quantities.h"

namespace sim {

double last_dt;

double PhysicalTime, PhysicalTimeInitial;
unsigned int N_snapshot = 0;
unsigned int N_monitor = 0;
unsigned long int N_hydro_iter = 0;

double total_disk_mass_old;


static void write_snapshot(t_data &data) {
	// Outputs are done here
	output::last_snapshot_dir = output::snapshot_dir;
	output::write_full_output(data, std::to_string(N_snapshot));
	output::cleanup_autosave();

	if (N_snapshot == 0 && parameters::damping) {
	// Write damping data as a reference.
	const std::string snapshot_dir_old = output::snapshot_dir;
	output::write_full_output(data, "damping", false);
	output::snapshot_dir = snapshot_dir_old;
	}
}

void handle_outputs(t_data &data) {
	bool need_update_for_output = true;
	N_snapshot = (N_monitor / parameters::NINTERM); // note: integer division
	
	const bool to_write_snapshot = (parameters::NINTERM * N_snapshot == N_monitor);
	const bool to_write_monitor = to_write_snapshot || parameters::write_at_every_timestep;


	/// asure planet torques are computed
	if (!parameters::disk_feedback && to_write_monitor) {
	    ComputeDiskOnNbodyAccel(data);
	}

	if (to_write_snapshot) {
	    need_update_for_output = false;
		write_snapshot(data);
	}

	if (to_write_snapshot && parameters::write_torques) {
	    output::write_torques(data, need_update_for_output);
	}

	if (to_write_monitor) {
		dt_logger.write(N_snapshot, N_monitor);
		if(ECC_GROWTH_MONITOR){
			output::write_ecc_peri_changes(sim::N_snapshot, sim::N_monitor);
		}
		output::write_monitor_time();
	    ComputeCircumPlanetaryMasses(data);
	    data.get_planetary_system().write_planets(1);

		if (parameters::write_lightcurves) {
			output::write_lightcurves(data, N_snapshot, need_update_for_output);
		}
	}

	// write disk quantities like eccentricity, ...
	if ((to_write_monitor) && parameters::write_disk_quantities) {
	    output::write_quantities(data, need_update_for_output);
	}

}

double CalculateTimeStep(t_data &data)
{
	double rv = last_dt;

	if (parameters::calculate_disk) {
		const double cfl_dt = cfl::condition_cfl(data);
		rv = std::min(parameters::CFL_max_var * last_dt, cfl_dt);
		last_dt = cfl_dt;

		if(PRINT_SIG_INFO){
			cfl::condition_cfl(data, cfl_dt);
		}

	}
	dt_logger.update(rv);

	return rv;
}


// static double PlanNextTimestepSize(t_data &data, double dt, double force_calc)
// {

//     if (!SloppyCFL || force_calc) {
// 	sim::last_dt = dt;
// 	const double dt_cfl = cfl::condition_cfl(data, 0.0);

// 	if(PRINT_SIG_INFO){
// 		cfl::condition_cfl(data, dt_cfl);
// 	}

// 	// don't let dt grow too fast
// 	const double limited_dt = std::min(parameters::CFL_max_var * last_dt, dt_cfl);

// 	// Limit dt un such a way, that we precisely end up on DT
// 	const double deltaT = DT - dtemp; // time till full DT
// 	const double inverse_dt_limited = std::max(deltaT / limited_dt, 1.0);
// 	dt = deltaT / inverse_dt_limited;
//     }
//     return dt;
// }

/*
Do one step of the integration.

Returns the time covered.
*/
static void step_Euler(t_data &data, const double dt) {

	const double time = PhysicalTime;

	if (parameters::disk_feedback) {
	    ComputeDiskOnNbodyAccel(data);
	}

	if (parameters::disk_feedback) {
	    UpdatePlanetVelocitiesWithDiskForce(data, dt);
	}

	refframe::ComputeIndirectTermDisk(data);
	refframe::ComputeIndirectTermNbody(data, time, dt);
	refframe::ComputeIndirectTermFully();

	data.get_planetary_system().apply_indirect_term_on_Nbody(
		refframe::IndirectTerm, dt);


	if (parameters::integrate_planets) {
		data.get_planetary_system().integrate(time, dt);
		/// Nbody positions and velocities are not updated yet!
	}

	if (parameters::calculate_disk) {
		/** Gravitational potential from star and planet(s) is computed and
		 * stored here*/
		if (parameters::body_force_from_potential) {
		CalculateNbodyPotential(data, time);
		} else {
		CalculateAccelOnGas(data, time);
		}
	}

	if (parameters::integrate_particles) {
		particles::integrate(data, time, dt);
	}

	/** Planets' positions and velocities are updated from gravitational
	 * interaction with star and other planets */
	if (parameters::integrate_planets) {
		data.get_planetary_system().copy_data_from_rebound();
		data.get_planetary_system().move_to_hydro_center_and_update_orbital_parameters();
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
		art_visc::update_with_artificial_viscosity(data, sim::PhysicalTime, dt);
		if (parameters::Adiabatic) {
		    SetTemperatureFloorCeilValues(data, __FILE__, __LINE__);
		}
		recalculate_viscosity(data, sim::PhysicalTime);
		viscosity::compute_viscous_stress_tensor(data);
		viscosity::update_velocities_with_viscosity(data, dt);
	    }

	    if (!parameters::EXPLICIT_VISCOSITY) {
		Sts(data, time, dt);
	    }

	    if (parameters::Adiabatic) {
		SubStep3(data, time, dt);
		if (parameters::radiative_diffusion_enabled) {
		    radiative_diffusion(data, time, dt);
		}
	    }

	    /// TODO moved apply boundaries here
		boundary_conditions::apply_boundary_condition(data, time, 0.0, false);

		Transport(data, &data[t_data::SIGMA], &data[t_data::V_RADIAL],
			  &data[t_data::V_AZIMUTHAL], &data[t_data::ENERGY],
			  dt);
	}

	// TODO: move outside step
	PhysicalTime += dt;
	N_hydro_iter = N_hydro_iter + 1;
	logging::print_runtime_info();

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

	    boundary_conditions::apply_boundary_condition(data, time, dt, true);

	    // const double total_disk_mass_new =
	    //  quantities::gas_total_mass(data, 2.0*RMAX);

	    // data[t_data::DENSITY] *=
	    //(total_disk_mass_old / total_disk_mass_new);

	    quantities::CalculateMonitorQuantitiesAfterHydroStep(data, N_monitor,
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


/**
	\param data
	\param sys
	compare Bodenheimer et al. Numerical Methods in Astrophyiscs
	Chapter 6.3.4 Operator Splitting
	Gas:	kick drift kick
	Nbody:	kick drift kick
*/
[[maybe_unused]] static void step_LeapFrog(t_data &data, const double step_dt)
{
	const double frog_dt = step_dt/2;
	const double start_time = PhysicalTime;
	const double end_time = PhysicalTime + step_dt;

	//////////////// Sourceterms first half /////////////////////
	if (parameters::disk_feedback) {
	    ComputeDiskOnNbodyAccel(data);
	}
	refframe::ComputeIndirectTermDisk(data);
	/// Indirect term will not be updated for the second leapfrog step
	/// so compute it for the full timestep
	/// It should be recomputed when using euler though
	refframe::ComputeIndirectTermNbody(data, start_time, step_dt);
	refframe::ComputeIndirectTermFully();

	/// Sourceterms Nbody
	if (parameters::integrate_planets) {
		if (parameters::disk_feedback) {
			// This version of disk feedback is exactly the following leapfrog method:
			// v_1/2 = v_0   + dt/2 a(r_0)
			// r_1 =   r_0   + dt v_1/2
			// v_1 =   v_1/2 + dt/2 a(r_1)
			UpdatePlanetVelocitiesWithDiskForce(data, frog_dt);
		}
		data.get_planetary_system().apply_indirect_term_on_Nbody(refframe::IndirectTerm, frog_dt);
	}

	/// Sourceterms dust
	if (parameters::integrate_particles) {
		particles::update_velocities_from_indirect_term(frog_dt);
	}

	/// Sourceterms Gas
	if (parameters::calculate_disk) {
		/// Use Nbody at x_i+1/2 for gas interaction
		if (parameters::body_force_from_potential) {
		CalculateNbodyPotential(data, start_time);
		} else {
		CalculateAccelOnGas(data, start_time);
		}

		update_with_sourceterms(data, frog_dt);

	    if (parameters::EXPLICIT_VISCOSITY) {

		art_visc::update_with_artificial_viscosity(data, start_time, frog_dt);
		if (parameters::Adiabatic) {
			SetTemperatureFloorCeilValues(data, __FILE__, __LINE__);
		}

		recalculate_viscosity(data, start_time);
		viscosity::compute_viscous_stress_tensor(data);
		viscosity::update_velocities_with_viscosity(data, frog_dt);
	    }
	    if (!parameters::EXPLICIT_VISCOSITY) {
		Sts(data, start_time, frog_dt);
	    }

	    if (parameters::Adiabatic) {
		SubStep3(data, start_time, frog_dt);
		if (parameters::radiative_diffusion_enabled) {
			radiative_diffusion(data, start_time, frog_dt);
		}
	    }
	}
		//////////////// END Sourceterms first half /////////////////////


		//////////////// Transport full step  /////////////////////

		if (parameters::integrate_planets) {
			refframe::init_corotation(data);
			data.get_planetary_system().integrate(start_time, step_dt);
			data.get_planetary_system().copy_data_from_rebound();

			data.get_planetary_system().compute_dist_to_primary();
			data.get_planetary_system().calculate_orbital_elements();
		}

		if (parameters::integrate_particles) {
			particles::integrate(data, start_time, frog_dt);
		}

		if (parameters::calculate_disk) {
		boundary_conditions::apply_boundary_condition(data, start_time, 0.0, false);

		Transport(data, &data[t_data::SIGMA], &data[t_data::V_RADIAL],
			  &data[t_data::V_AZIMUTHAL], &data[t_data::ENERGY],
			  step_dt);
		}

		// Below we correct v_azimuthal, planet's position and velocities if we
		// work in a frame non-centered on the star. Same for dust particles.
		refframe::handle_corotation(data, frog_dt);
		//////////////// Transport full step   /////////////////////

	//////////////// Sourceterms second half /////////////////////


	if (parameters::disk_feedback) {
		ComputeDiskOnNbodyAccel(data, true);
	}
	refframe::ComputeIndirectTermDisk(data);

	if(parameters::indirect_term_mode == INDIRECT_TERM_EULER){
	refframe::ComputeIndirectTermNbody(data, end_time, step_dt);
	}
	refframe::ComputeIndirectTermFully();

	/// Sourceterms Gas
	if (parameters::calculate_disk) {

		if (parameters::body_force_from_potential) {
		CalculateNbodyPotential(data, end_time);
		} else {
		CalculateAccelOnGas(data, end_time);
		}

		compute_pressure(data);
		if(parameters::self_gravity){
			compute_scale_height(data, end_time);
		}
		update_with_sourceterms(data, frog_dt);

		if (parameters::EXPLICIT_VISCOSITY) {
		art_visc::update_with_artificial_viscosity(data, end_time, frog_dt);

		recalculate_viscosity(data, end_time);
		viscosity::compute_viscous_stress_tensor(data);
		viscosity::update_velocities_with_viscosity(data, frog_dt);
		}
		if (!parameters::EXPLICIT_VISCOSITY) {
		Sts(data, end_time, frog_dt);
		}

		if (parameters::Adiabatic) {
		SubStep3(data, end_time, frog_dt);
		if (parameters::radiative_diffusion_enabled) {
			radiative_diffusion(data, end_time, frog_dt);
		}
		}
	}

	/// Sourceterms dust
	if (parameters::integrate_particles) {
	particles::update_velocities_from_indirect_term(frog_dt);
	}

	/// Sourceterms Nbody
	if (parameters::integrate_planets) {
		if (parameters::disk_feedback) {
			UpdatePlanetVelocitiesWithDiskForce(data, frog_dt);
		}
		data.get_planetary_system().apply_indirect_term_on_Nbody(refframe::IndirectTerm, frog_dt);
		data.get_planetary_system().move_to_hydro_center_and_update_orbital_parameters();
	}

	///////////// END Nbody update  ///////////////////

	//////////////// END Sourceterms second half   /////////////////////

	PhysicalTime = end_time;
	N_hydro_iter += 1;
	logging::print_runtime_info();

	if (parameters::calculate_disk) {
	    CommunicateBoundaries(
		&data[t_data::SIGMA], &data[t_data::V_RADIAL],
		&data[t_data::V_AZIMUTHAL], &data[t_data::ENERGY]);

	    // We only recompute once, assuming that cells hit by planet
	    // accretion are not also hit by viscous accretion at inner
	    // boundary.
	    if (parameters::VISCOUS_ACCRETION) {
		recalculate_viscosity(data, end_time);
	    }

		// minimum density is assured inside AccreteOntoPlanets
	    accretion::AccreteOntoPlanets(data, step_dt);

		boundary_conditions::apply_boundary_condition(data, end_time, step_dt, true);

		if(parameters::keep_mass_constant){
			const double total_disk_mass_new =
			quantities::gas_total_mass(data, RMAX);
			 data[t_data::SIGMA] *=
			(total_disk_mass_old / total_disk_mass_new);
		}

		quantities::CalculateMonitorQuantitiesAfterHydroStep(data, N_monitor,
							 step_dt);

		if (parameters::variableGamma &&
		!parameters::VISCOUS_ACCRETION) { // If VISCOUS_ACCRETION is active,
					  // scale_height is already updated
		// Recompute scale height after Transport to update the 3D
		// density
		compute_sound_speed(data, end_time);
		compute_scale_height(data, end_time);
		}
		// this must be done after CommunicateBoundaries
		recalculate_derived_disk_quantities(data, end_time);

	}
}

///
/// \brief step_LeapFrog_old
/// Gas:	kick drift kick
/// Nbody:	drift kick drift
/// \param data
/// \param step_dt
///
[[maybe_unused]] static void step_LeapFrog_old(t_data &data, const double step_dt)
{
	const double frog_dt = step_dt/2;
	const double start_time = PhysicalTime;
	const double midstep_time = PhysicalTime + frog_dt;
	const double end_time = PhysicalTime + step_dt;

	//////////////// Leapfrog compute v_i+1/2 /////////////////////

	refframe::ComputeIndirectTermNbody(data, start_time, step_dt);
	if (parameters::integrate_planets) { //// Nbody drift / 2
	refframe::init_corotation(data);
	data.get_planetary_system().integrate(start_time, frog_dt);
	data.get_planetary_system().copy_data_from_rebound();
	data.get_planetary_system().compute_dist_to_primary();
	data.get_planetary_system().calculate_orbital_elements();
	}

	if (parameters::disk_feedback) {
		ComputeDiskOnNbodyAccel(data);
	}
	refframe::ComputeIndirectTermDisk(data);

	/// Indirect term will not be updated for the second leapfrog step
	/// so compute it for the full timestep
	/// It should be recomputed when using euler though
	refframe::ComputeIndirectTermFully();

	if (parameters::integrate_planets) { /// Nbody Kick 1 / 2
		if (parameters::disk_feedback) {
			UpdatePlanetVelocitiesWithDiskForce(data, frog_dt);
		}
		data.get_planetary_system().apply_indirect_term_on_Nbody(refframe::IndirectTerm, frog_dt);
	}

	if (parameters::integrate_particles) {
		particles::integrate(data, start_time, frog_dt);
		particles::update_velocities_from_indirect_term(frog_dt);
	}

	refframe::handle_corotation(data, frog_dt);

	if (parameters::calculate_disk) {
		/// Gas Kick 1/2
		if (parameters::body_force_from_potential) {
		CalculateNbodyPotential(data, start_time);
		} else {
		CalculateAccelOnGas(data, start_time);
		}

		update_with_sourceterms(data, frog_dt);

		if (parameters::EXPLICIT_VISCOSITY) {

		art_visc::update_with_artificial_viscosity(data, start_time, frog_dt);
		if (parameters::Adiabatic) {
			SetTemperatureFloorCeilValues(data, __FILE__, __LINE__);
		}
		//recalculate_viscosity(data, start_time);
		viscosity::compute_viscous_stress_tensor(data);
		viscosity::update_velocities_with_viscosity(data, frog_dt);
		}
		if (!parameters::EXPLICIT_VISCOSITY) {
		Sts(data, start_time, frog_dt);
		}

		if (parameters::Adiabatic) {
		SubStep3(data, start_time, frog_dt);
		if (parameters::radiative_diffusion_enabled) {
			radiative_diffusion(data, start_time, frog_dt);
		}
		}
		//////////////// END /// Gas Kick 1/2 /////////////////////


		//////////////// Gas drift 1/1 /////////////////////
		boundary_conditions::apply_boundary_condition(data, start_time, 0.0, false);

		Transport(data, &data[t_data::SIGMA], &data[t_data::V_RADIAL],
			  &data[t_data::V_AZIMUTHAL], &data[t_data::ENERGY],
			  step_dt);
		//////////////// END Gas drift 1/1   /////////////////////

	}

	//////////////// Gas kick 2/2   /////////////////////
	/// planets positions still at x_i+1/2 for gas interaction
	if (parameters::disk_feedback) {
		ComputeDiskOnNbodyAccel(data, true);
	}
	refframe::ComputeIndirectTermDisk(data);

	if(parameters::indirect_term_mode == INDIRECT_TERM_EULER){
	refframe::ComputeIndirectTermNbody(data, midstep_time, step_dt);
	}
	refframe::ComputeIndirectTermFully();

	/// update gas while Nbody positions are still at x_i+1/2
	if (parameters::calculate_disk) {

		if (parameters::body_force_from_potential) {
		CalculateNbodyPotential(data, midstep_time);
		} else {
		CalculateAccelOnGas(data, midstep_time);
		}

		compute_pressure(data);
		if(parameters::self_gravity){
			compute_scale_height(data, midstep_time);
		}
		update_with_sourceterms(data, frog_dt);

		if (parameters::EXPLICIT_VISCOSITY) {
		art_visc::update_with_artificial_viscosity(data, midstep_time, frog_dt);
		if (parameters::Adiabatic) {
			SetTemperatureFloorCeilValues(data, __FILE__, __LINE__);
		}

		//recalculate_viscosity(data, midstep_time);
		viscosity::compute_viscous_stress_tensor(data);
		viscosity::update_velocities_with_viscosity(data, frog_dt);
		}
		if (!parameters::EXPLICIT_VISCOSITY) {
		Sts(data, midstep_time, frog_dt);
		}

		if (parameters::Adiabatic) {
		SubStep3(data, midstep_time, frog_dt);
		if (parameters::radiative_diffusion_enabled) {
			radiative_diffusion(data, midstep_time, frog_dt);
		}
		}
	}

	/// We update particles with Nbody at x_i+1/2
	/// and gas at x_i/v_i, so we use gas at x_i+1/v_i+1 to finish the update step
	if (parameters::integrate_particles) {
	particles::update_velocities_from_indirect_term(frog_dt);
	particles::integrate(data, midstep_time, frog_dt);
	}

	// Finish timestep of the planets but do not update Nbody system yet //
	if (parameters::integrate_planets) {
		/// Nbody kick 2/2
		if (parameters::disk_feedback) {
			UpdatePlanetVelocitiesWithDiskForce(data, frog_dt);
		}
		data.get_planetary_system().apply_indirect_term_on_Nbody(refframe::IndirectTerm, frog_dt);

		/// Nbody drift 2/2
		refframe::init_corotation(data);
		data.get_planetary_system().integrate(midstep_time, frog_dt);
		data.get_planetary_system().copy_data_from_rebound();
		data.get_planetary_system().move_to_hydro_center_and_update_orbital_parameters();
	}

	/* Below we correct v_azimuthal, planet's position and velocities if we
	 * work in a frame non-centered on the star. Same for dust particles. */
	refframe::handle_corotation(data, frog_dt);
	///////////// END Nbody update  ///////////////////

	//////////////// END Leapfrog compute v_i+1   /////////////////////

	PhysicalTime = end_time;
	N_hydro_iter += 1;
	logging::print_runtime_info();

	if (parameters::calculate_disk) {
		CommunicateBoundaries(
		&data[t_data::SIGMA], &data[t_data::V_RADIAL],
		&data[t_data::V_AZIMUTHAL], &data[t_data::ENERGY]);

		// We only recompute once, assuming that cells hit by planet
		// accretion are not also hit by viscous accretion at inner
		// boundary.
		if (parameters::VISCOUS_ACCRETION) {
		compute_sound_speed(data, end_time);
		compute_scale_height(data, end_time);
		viscosity::update_viscosity(data);
		}

		// minimum density is assured inside AccreteOntoPlanets
		accretion::AccreteOntoPlanets(data, step_dt);

		boundary_conditions::apply_boundary_condition(data, end_time, step_dt, true);

		if(parameters::keep_mass_constant){
			const double total_disk_mass_new =
			quantities::gas_total_mass(data, RMAX);
			 data[t_data::SIGMA] *=
			(total_disk_mass_old / total_disk_mass_new);
		}

		quantities::CalculateMonitorQuantitiesAfterHydroStep(data, N_monitor,
							 step_dt);

		if (parameters::variableGamma &&
		!parameters::VISCOUS_ACCRETION) { // If VISCOUS_ACCRETION is active,
					  // scale_height is already updated
		// Recompute scale height after Transport to update the 3D
		// density
		compute_sound_speed(data, end_time);
		compute_scale_height(data, end_time);
		}
		// this must be done after CommunicateBoundaries
		recalculate_derived_disk_quantities(data, end_time);

	}
}


void init(t_data &data) {
	boundary_conditions::apply_boundary_condition(data, PhysicalTime, 0.0, false); // TODO: move to init
	refframe::init_corotation(data); // TODO: move to init

	if (start_mode::mode != start_mode::mode_restart) {
		CalculateTimeStep(data);
	}

	if (parameters::calculate_disk) {
	CommunicateBoundaries(&data[t_data::SIGMA], &data[t_data::V_RADIAL],
			      &data[t_data::V_AZIMUTHAL],
			      &data[t_data::ENERGY]);
    }

	total_disk_mass_old = 1.0;
	if(parameters::keep_mass_constant){
	total_disk_mass_old =
	quantities::gas_total_mass(data, RMAX);
	}


}

static void step(t_data &data, const double step_dt) {
	if (parameters::leap_frog) {
		step_LeapFrog(data, step_dt);
	} else {
		step_Euler(data, step_dt);
	}

}

static void handle_signals(t_data &data) {
	if (SIGTERM_RECEIVED) {
		output::write_full_output(data, "autosave");
		PersonalExit(0);
	}
}

void run(t_data &data) {

	double step_dt = last_dt;
	double cfl_dt = last_dt;

	init(data);

	const double t_final = parameters::NTOT * parameters::NINTERM * parameters::DT;

	bool towrite = false;

    for (; PhysicalTime < t_final; N_hydro_iter++) {
		
		handle_signals(data);

		cfl_dt = CalculateTimeStep(data);

		const double time_next_monitor = (N_monitor+1)*parameters::DT;
		const double time_left_till_write = time_next_monitor - PhysicalTime;

		if (cfl_dt > time_left_till_write) {
			step_dt = time_left_till_write;
			N_monitor++;
			towrite = true;
		} else {
			step_dt = cfl_dt;
			towrite = false;
		}

		step(data, step_dt);

		if (towrite) {
			handle_outputs(data);
			logging::print_runtime_info();
		}
    }


	logging::print_runtime_final();

}





} // close namespace sim
