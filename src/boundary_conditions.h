#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include "data.h"
#include "types.h"

namespace boundary_conditions
{

void init_prescribed_time_variable_boundaries(t_data &data);
void boundary_condition_precribed_time_variable_outer(t_data &data,
							  t_polargrid *densitystar, const double current_time);
void apply_boundary_condition_temperature(t_data &data);
void apply_boundary_condition(t_data &data, const double current_time, const double dt, const bool final);
void open_boundary_inner(t_data &data);
void open_boundary_outer(t_data &data);
void zero_gradient_boundary_inner(t_data &data);
void zero_gradient_boundary_outer(t_data &data);
void reflecting_boundary_inner(t_data &data);
void reflecting_boundary_outer(t_data &data);
void viscous_outflow_boundary_inner(t_data &data);
void damping(t_data &data, double dt);
void mass_overflow(t_data &data, const double current_time);
void mass_overflow_willy(t_data &data, t_polargrid *densitystar,
			 bool transport);
void boundary_layer_inner_boundary(t_data &data);
void boundary_layer_outer_boundary(t_data &data);
void keplerian2d_boundary_inner(t_data &data);
void keplerian2d_boundary_outer(t_data &data);

void initial_center_of_mass_boundary_inner(t_data &data);
void damping_initial_center_of_mass_inner(t_data &data, double dt);
void initial_center_of_mass_boundary_outer(t_data &data);
void damping_initial_center_of_mass_outer(t_data &data, double dt);

void jibin_boundary_inner(t_data &data);
void jibin_boundary_outer(t_data &data);

void damping_single_inner(t_polargrid &quantity, t_polargrid &quantity0,
			  double dt);
void damping_single_outer(t_polargrid &quantity, t_polargrid &quantity0,
			  double dt);
void damping_single_inner_zero(t_polargrid &quantity, t_polargrid &quantity0,
			       double dt);
void damping_single_outer_zero(t_polargrid &quantity, t_polargrid &quantity0,
			       double dt);
void damping_single_inner_mean(t_polargrid &quantity, t_polargrid &quantity0,
			       double dt);
void damping_single_outer_mean(t_polargrid &quantity, t_polargrid &quantity0,
			       double dt);
void damping_vradial_inner_visc(t_polargrid &vrad, t_polargrid &viscosity,
				double dt);
} // namespace boundary_conditions

#endif // BOUNDARY_CONDITIONS_H
