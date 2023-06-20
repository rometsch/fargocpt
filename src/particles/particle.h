#pragma once

#include "../constants.h"
#include "../util.h"
#include <cmath>

class t_particle
{
  public:
    /// id (e.g. for tracking)
    unsigned int id;
    /// r position
    double r;
    /// phi position
    double phi;
    /// r derivative
    double r_dot;
    /// phi derivative
    double phi_dot;
    /// r dot derivative
    double r_ddot;
    /// phi dot derivative
    double phi_ddot;
    /// mass
    double mass;
    /// radius
    double radius;
    /// adaptive timestep
    double timestep;
    /// last error for timestep estimation
    double facold;
    /// stokes number
    double stokes;


    double get_squared_distance_to_star();
    double get_distance_to_star();
    double get_angle(void) const;
    double get_r_dot(void) const;
    double get_phi_dot(void) const;
};
