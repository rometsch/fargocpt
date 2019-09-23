#include "particle.h"
#include <math.h>
#include "global.h"

double t_particle::get_squared_distance_to_star()
{
	if(CartesianParticles)
	{
		const double &x = r;
		const double &y = phi;
		return pow2(x) + pow2(y);
	}
	else
	{
		return r*r;
	}
}

double t_particle::get_distance_to_star()
{
	if(CartesianParticles)
	{
		return sqrt(get_squared_distance_to_star());
	}
	else
	{
		return r;
	}
}

double t_particle::get_angle() const
{

	if(CartesianParticles)
	{
		const double &x = r;
		const double &y = phi;

		double phi_ = atan2(y,x);

		if(phi_<0)
		{
			phi_ += 2*PI;
		}

		return phi_;
	}
	else {
		return phi;
	}
}



double t_particle::get_r_dot() const
{

	if(CartesianParticles)
	{
		const double &x  = r;
		const double &y  = phi;
		const double &vx = r_dot;
		const double &vy = phi_dot;

		const double local_r = sqrt(x*x+y*y);

		const double local_r_dot = (x*vx + y*vy) / local_r;
		return local_r_dot;
	}
	else {
		return r_dot;
	}
}

double t_particle::get_phi_dot() const
{

	if(CartesianParticles)
	{
		const double &x  = r;
		const double &y  = phi;
		const double &vx = r_dot;
		const double &vy = phi_dot;

		const double local_r = sqrt(x*x+y*y);

		const double local_phi_dot = (x*vy - vx*y) / (local_r*local_r);
		return local_phi_dot;
	}
	else {
		return phi_dot;
	}
}
