#include "particle.h"
#include <math.h>

double t_particle::get_squared_distance_to_star()
{
    return pow2(x) + pow2(y);
}

double t_particle::get_distance_to_star()
{
    return sqrt( get_squared_distance_to_star() );
}

double t_particle::get_angle()
{
    double phi = atan2(y,x);

    if(phi<0)
    {
        phi += 2*PI;
    }

    return phi;
}
