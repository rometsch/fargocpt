#include "frame_of_reference.h"
#include "parameters.h"
#include "SideEuler.h"
#include "particles/particles.h"
#include "types.h"
#include "Pframeforce.h"

namespace refframe {

double planet_corot_ref_old_x;
double planet_corot_ref_old_y;
double OmegaFrame;
double FrameAngle;

Pair IndirectTerm;
Pair IndirectTermDisk;
Pair IndirectTermPlanets;

void init_corotation(t_data &data)
{
    if (parameters::corotating == YES) {
	// save old planet positions
	const unsigned int n = parameters::corotation_reference_body;
	const auto &planet = data.get_planetary_system().get_planet(n);
	planet_corot_ref_old_x = planet.get_x();
	planet_corot_ref_old_y = planet.get_y();
    }
}

void handle_corotation(t_data &data, const double dt)
{
    if (parameters::corotating == YES) {
	unsigned int n = parameters::corotation_reference_body;
	auto &planet = data.get_planetary_system().get_planet(n);
	const double x = planet.get_x();
	const double y = planet.get_y();
	const double distance_new = std::sqrt(std::pow(x, 2) + std::pow(y, 2));
	const double distance_old =
	    std::sqrt(std::pow(planet_corot_ref_old_x, 2) + std::pow(planet_corot_ref_old_y, 2));
	const double cross = planet_corot_ref_old_x * y - x * planet_corot_ref_old_y;

	// new = r_new x r_old = distance_new * distance_old * sin(alpha*dt)
	const double OmegaNew =
		std::asin(cross / (distance_new * distance_old)) / dt;

	const double domega = (OmegaNew - OmegaFrame);
	if (parameters::calculate_disk) {
	    correct_v_azimuthal(data[t_data::V_AZIMUTHAL], domega);
	}
	OmegaFrame = OmegaNew;
    }

    data.get_planetary_system().rotate(OmegaFrame * dt);

    if (parameters::integrate_particles) {
	particles::rotate(OmegaFrame * dt);
    }

    FrameAngle += OmegaFrame * dt;
}


/**
 * @brief ComputeIndirectTerm: IndirectTerm is the correction therm that needs
 * to be added to the accelerations.
 * @param force
 * @param data
 */
void ComputeIndirectTermDisk(t_data &data)
{
    IndirectTermDisk.x = 0.0;
    IndirectTermDisk.y = 0.0;

    // compute disk indirect term
    if (parameters::disk_feedback) {
	// add up contributions from disk on all bodies used to calculate the
	// center
	double mass_center = 0.0;
	for (unsigned int n = 0; n < parameters::n_bodies_for_hydroframe_center;
	     n++) {
	    t_planet &planet = data.get_planetary_system().get_planet(n);
	    const double mass = planet.get_mass();
	    const Pair accel = planet.get_disk_on_planet_acceleration();
	    IndirectTermDisk.x -= mass * accel.x;
	    IndirectTermDisk.y -= mass * accel.y;
	    mass_center += mass;
	}
	IndirectTermDisk.x /= mass_center;
	IndirectTermDisk.y /= mass_center;
    }
}

void ComputeIndirectTermNbodyEuler(t_data &data)
{
    IndirectTermPlanets.x = 0.0;
    IndirectTermPlanets.y = 0.0;

    // compute nbody indirect term
    // add up contributions from mutual interactions from all bodies used to
    // calculate the center
    double mass_center = 0.0;
    for (unsigned int n = 0; n < parameters::n_bodies_for_hydroframe_center;
	 n++) {
	t_planet &planet = data.get_planetary_system().get_planet(n);
	const double mass = planet.get_mass();
	const Pair accel = planet.get_nbody_on_planet_acceleration();
	IndirectTermPlanets.x -= mass * accel.x;
	IndirectTermPlanets.y -= mass * accel.y;
	mass_center += mass;
    }
    IndirectTermPlanets.x /= mass_center;
	IndirectTermPlanets.y /= mass_center;
}

void ComputeIndirectTermNbody(t_data &data, const double current_time, const double dt)
{

	if(parameters::n_bodies_for_hydroframe_center == data.get_planetary_system().get_number_of_planets()){
	IndirectTermPlanets.x = 0.0;
	IndirectTermPlanets.y = 0.0;
        return;
	}

	if(parameters::indirect_term_mode == INDIRECT_TERM_EULER){
		ComputeNbodyOnNbodyAccel(data.get_planetary_system());
		ComputeIndirectTermNbodyEuler(data);
	} else {

	if(dt != 0.0){ // Indirect term from Rebound
	/// compute the Indirect term as the effective acceleration from a high order nbody integrator.
	/// this typically leads to vel_center ~ 0/0 but pos_center != 0/0, but shifting the center to 0.0 causes an error
	/// because the gas does not feel the kick
	data.get_planetary_system().copy_data_to_rebound();
	data.get_planetary_system().m_rebound->t = current_time;
	const pair delta_vel = data.get_planetary_system().get_hydro_frame_center_delta_vel_rebound_predictor(dt);
	pair accel{delta_vel.x/dt, delta_vel.y/dt};

	IndirectTermPlanets.x = -accel.x;
	IndirectTermPlanets.y = -accel.y;
	} else {
	IndirectTermPlanets.x = 0.0;
	IndirectTermPlanets.y = 0.0;
	}
	}
}

void ComputeIndirectTermFully(){
	IndirectTerm.x = IndirectTermDisk.x + IndirectTermPlanets.x;
	IndirectTerm.y = IndirectTermDisk.y + IndirectTermPlanets.y;
}



}