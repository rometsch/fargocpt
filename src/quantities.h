#ifndef QUANTITIES_H
#define QUANTITIES_H

#include "data.h"

namespace quantities
{

double gas_total_mass(t_data &data, const double quantitiy_radius);
double gas_aspect_ratio(t_data &data, const double quantitiy_radius);
double gas_disk_radius(t_data &data, const double total_mass);
double gas_angular_momentum(t_data &data, const double quantitiy_radius);
double gas_internal_energy(t_data &data, const double quantitiy_radius);
double gas_viscous_dissipation(t_data &data, const double quantitiy_radius);
double gas_luminosity(t_data &data, const double quantitiy_radius);
double gas_kinematic_energy(t_data &data, const double quantitiy_radius);
double gas_radial_kinematic_energy(t_data &data, const double quantitiy_radius);
double gas_azimuthal_kinematic_energy(t_data &data,
				      const double quantitiy_radius);
double gas_gravitational_energy(t_data &data, const double quantitiy_radius);

void calculate_disk_quantities(t_data &data);
void calculate_alpha_grav(t_data &data);
void calculate_alpha_grav_mean_sumup(t_data &data,
					 double dt);
void calculate_alpha_reynolds(t_data &data);
void calculate_alpha_reynolds_mean_sumup(t_data &data,
					 double dt);
void calculate_toomre(t_data &data);
void calculate_radial_luminosity(t_data &data);
void calculate_radial_dissipation(t_data &data);
void calculate_massflow(t_data &data);
void compute_aspectratio(t_data &data);
void calculate_advection_torque(t_data &data);
void calculate_gravitational_torque(t_data &data);
void calculate_viscous_torque(t_data &data);
} // namespace quantities

#endif // QUANTITIES_H
