#pragma once

#include "constants.h"
#include "global.h"
#include "types.h"
#include "nbody/planetary_system.h"

double initial_energy(const double R, const double M);
double initial_viscous_radial_speed(const double R, const double M);
double initial_locally_isothermal_smoothed_v_az(const double R, const double M);
double initial_locally_isothermal_smoothed_v_az_with_quadropole_moment(const double R, const double M);
void init_binary_quadropole_moment(const t_planetary_system &psys);
double initial_locally_isothermal_v_az(const double R, const double M);
double compute_v_kepler(const double R, const double M);
double calculate_omega_kepler(double r);
void RefillSigma(t_polargrid *Density);
void RefillEnergy(t_polargrid *Energy);
double eggleton_1983(const double q, const double r);
double init_l1(const double central_star_mass, const double other_star_mass);
void update_l1(const double central_star_mass, const double m2, double &l1);
