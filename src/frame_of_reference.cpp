#include "frame_of_reference.h"
#include "parameters.h"
#include "SideEuler.h"
#include "global.h"
#include "particles/particles.h"
#include "types.h"

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

    if (parameters::integrate_planets) {
	data.get_planetary_system().rotate(OmegaFrame * dt);
    }
    if (parameters::integrate_particles) {
	particles::rotate(OmegaFrame * dt);
    }

    FrameAngle += OmegaFrame * dt;
}

}