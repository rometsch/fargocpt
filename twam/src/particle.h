#ifndef _PARTICLE_H_
#define _PARTICLE_H_

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

		inline double get_squared_distance_to_star(void) { return pow2(x)+pow2(y); };
		inline double get_distance_to_star(void) { return sqrt(pow2(x)+pow2(y)); };
		inline double get_angle(void) { double phi = atan2(y, x); if (phi < 0) { phi += 2*PI; } return phi; };
};

#endif