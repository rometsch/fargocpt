#include "particle.h"
#include "global.h"
#include <math.h>
#include "parameters.h"

double t_particle::get_squared_distance_to_star()
{
    if (parameters::CartesianParticles) {
	const double &x = r;
	const double &y = phi;
	return std::pow(x, 2) + std::pow(y, 2);
    } else {
	return r * r;
    }
}

double t_particle::get_distance_to_star()
{
    if (parameters::CartesianParticles) {
	return std::sqrt(get_squared_distance_to_star());
    } else {
	return r;
    }
}

double t_particle::get_angle() const
{

    if (parameters::CartesianParticles) {
	const double &x = r;
	const double &y = phi;

	double phi_ = std::atan2(y, x);

	if (phi_ < 0) {
	    phi_ += 2 * M_PI;
	}

	return phi_;
    } else {
	return phi;
    }
}

double t_particle::get_r_dot() const
{

    if (parameters::CartesianParticles) {
	const double &x = r;
	const double &y = phi;
	const double &vx = r_dot;
	const double &vy = phi_dot;

	const double local_r = std::sqrt(x * x + y * y);

	const double local_r_dot = (x * vx + y * vy) / local_r;
	return local_r_dot;
    } else {
	return r_dot;
    }
}

double t_particle::get_phi_dot() const
{

    if (parameters::CartesianParticles) {
	const double &x = r;
	const double &y = phi;
	const double &vx = r_dot;
	const double &vy = phi_dot;

	const double local_r = std::sqrt(x * x + y * y);

	const double local_phi_dot = (x * vy - vx * y) / (local_r * local_r);
	return local_phi_dot;
    } else {
	return phi_dot;
    }
}
