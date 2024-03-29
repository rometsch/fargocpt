#pragma once

#include "data.h"

namespace quantities
{

double gas_total_mass(t_data &data, const double quantitiy_radius);
double gas_quantity_reduce(const t_polargrid& arr, const double quantitiy_radius);
double gas_reduce_mass_average(t_data &data, const t_polargrid& arr, const double quantitiy_radius);
double gas_allreduce_mass_average(t_data &data, const t_polargrid& arr, const double quantitiy_radius);
double gas_disk_radius(t_data &data, const double total_mass);
double gas_angular_momentum(t_data &data, const double quantitiy_radius);
double gas_internal_energy(t_data &data, const double quantitiy_radius);
double gas_viscous_dissipation(t_data &data, const double quantitiy_radius);
double gas_luminosity(t_data &data, const double quantitiy_radius);
double gas_kinematic_energy(t_data &data, const double quantitiy_radius);
double gas_radial_kinematic_energy(t_data &data, const double quantitiy_radius);
double gas_azimuthal_kinematic_energy(t_data &data,
				      const double quantitiy_radius);

void fill_alpha_array(t_data &data, unsigned int timestep,
				   bool force_update);

void calculate_disk_ecc_peri(t_data &data, double &Ecc, double &Per, const bool force_update);
void calculate_disk_ecc_vector(t_data &data, unsigned int timestep, bool force_update);
void calculate_disk_delta_ecc_peri(t_data &data, double& dEcc, double& dPer);

void calculate_alpha_grav(t_data &data, unsigned int timestep,
			  bool force_update);
void calculate_alpha_grav_mean_sumup(t_data &data, unsigned int timestep,
				     double dt);
void calculate_alpha_reynolds(t_data &data, unsigned int timestep,
			      bool force_update);
void calculate_alpha_reynolds_mean_sumup(t_data &data, unsigned int timestep,
					 double dt);
void calculate_toomre(t_data &data, unsigned int timestep, bool force_update);
void calculate_radial_luminosity(t_data &data, unsigned int timestep,
				 bool force_update);
void calculate_radial_dissipation(t_data &data, unsigned int timestep,
				  bool force_update);
void calculate_massflow(t_data &data, unsigned int timestep, bool force_update);
void compute_aspectratio(t_data &data, unsigned int timestep,
			 bool force_update);
void calculate_advection_torque(t_data &data, unsigned int timestep,
				bool force_update);
void calculate_gravitational_torque(t_data &data, unsigned int timestep,
				    bool force_update);
void calculate_viscous_torque(t_data &data, unsigned int timestep,
			      bool force_update);
void CalculateMonitorQuantitiesAfterHydroStep(t_data &data,
				int nTimeStep, double dt);
void CalculateMonitorQuantitiesForOutput(t_data &data, double &tadv, double &tvisc, double &tgrav, const double quantities_limit_radius);
} // namespace quantities
