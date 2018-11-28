#ifndef PARTICLE_H
#define PARTICLE_H

#include <math.h>
#include "constants.h"
#include "util.h"

class t_particle
{
	public:
		/// id (e.g. for tracking)
		unsigned int id;
		/// x position
		double x;
		/// y position
		double y;
		/// x velocity
		double vx;
		/// y velocity
		double vy;
		/// x acceleration
		double ax;
		/// y acceleration
		double ay;
		/// mass
		double mass;
		/// radius
		double radius;

	  double get_squared_distance_to_star();
		double get_distance_to_star();
		double get_angle();
};

#endif // PARTICLE_H
