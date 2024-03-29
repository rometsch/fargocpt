#pragma once

#include "nbody/planetary_system.h"
#include "polargrid.h"

double initial_energy(const double R, const double M);
double initial_viscous_radial_speed(const double R, const double M);
double support_azi_quadrupole(const double R);
double support_azi_smoothing_derivative(const double R);
double support_azi_pressure(const double R);
double initial_locally_isothermal_smoothed_v_az(const double R, const double M);
double initial_locally_isothermal_smoothed_v_az_with_quadropole_moment(const double R, const double M);
void init_binary_quadropole_moment(const t_planetary_system &psys);
double initial_locally_isothermal_v_az(const double R, const double M);
double compute_v_kepler(const double R, const double M);
double calculate_omega_kepler(double r);
void compute_azi_avg_Sigma(t_polargrid &Density);
void compute_azi_avg_Energy(t_polargrid &Energy);
double eggleton_1983(const double q, const double r);
double init_l1(const double central_star_mass, const double other_star_mass);
double update_l1(const double central_star_mass, const double other_star_mass, double l1);
