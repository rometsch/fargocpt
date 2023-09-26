#pragma once

#include "../data.h"
#include "../types.h"

namespace boundary_conditions
{

/****************************************
/// Types of boundary conditions.
*****************************************/


extern void (*sigma_inner_func)(t_polargrid &, t_polargrid &);
extern void (*sigma_outer_func)(t_polargrid &, t_polargrid &);
extern void (*energy_inner_func)(t_polargrid &, t_polargrid &);
extern void (*energy_outer_func)(t_polargrid &, t_polargrid &);
extern void (*vrad_inner_func)(t_polargrid &, t_polargrid &);
extern void (*vrad_outer_func)(t_polargrid &, t_polargrid &);
extern void (*vaz_inner_func)(t_polargrid &, t_polargrid &);
extern void (*vaz_outer_func)(t_polargrid &, t_polargrid &);

extern std::string sigma_inner_name;
extern std::string sigma_outer_name;
extern std::string energy_inner_name;
extern std::string energy_outer_name;
extern std::string vrad_inner_name;
extern std::string vrad_outer_name;
extern std::string vaz_inner_name;
extern std::string vaz_outer_name;

void parse_config();

void old_apply_boundary_condition(t_data &data, const double current_time, const double dt, const bool final);
void apply_boundary_condition(t_data &data, const double current_time, const double dt, const bool final);

void zero_gradient_boundary_inner_single(t_polargrid &x);


/****************************************
/// Individual variable boundaries
****************************************/

void zero_gradient_inner(t_polargrid &x, t_polargrid &dummy);
void zero_gradient_outer(t_polargrid &x, t_polargrid &dummy);

void outflow_inner(t_polargrid &x, t_polargrid &dummy);
void outflow_outer(t_polargrid &x, t_polargrid &dummy);

void reference_value_inner(t_polargrid &x, t_polargrid &x0);
void reference_value_outer(t_polargrid &x, t_polargrid &x0);


/****************************************
/// Basic
****************************************/

void reference_value_boundary_inner(t_data &data);
void reference_value_boundary_outer(t_data &data);

void open_boundary_inner(t_data &data);
void open_boundary_outer(t_data &data);

void zero_gradient_boundary_inner(t_data &data);
void zero_gradient_boundary_outer(t_data &data);

void reflecting_boundary_inner(t_data &data);
void reflecting_boundary_outer(t_data &data);

void keplerian2d_boundary_inner(t_data &data);
void keplerian2d_boundary_outer(t_data &data);

void viscous_outflow_boundary_inner(t_data &data);

/****************************************
/// Azimuthal velocity
****************************************/

void ApplyKeplerianBoundaryInner(t_polargrid &v_azimuthal);
void ApplySubKeplerianBoundaryOuter(t_polargrid &v_azimuthal,
				    const bool did_sg);
void zero_shear_boundary(t_polargrid &v_azimuthal);

/****************************************
/// Prescribed value boundaries
****************************************/
void init_prescribed_time_variable_boundaries(t_data &data);
void boundary_condition_precribed_time_variable_outer(t_data &data,
							  t_polargrid *densitystar, const double current_time);

/****************************************
/// Inflow boundaries
****************************************/
void mass_overflow(t_data &data, const double current_time);
void mass_overflow_willy(t_data &data, t_polargrid *densitystar,
			 bool transport);
void boundary_layer_inner_boundary(t_data &data);
void boundary_layer_outer_boundary(t_data &data);


/****************************************
/// BC centered on center of mass
****************************************/
void initial_center_of_mass_boundary_inner(t_data &data);
void damping_initial_center_of_mass_inner(t_data &data, double dt);
void initial_center_of_mass_boundary_outer(t_data &data);
void damping_initial_center_of_mass_outer(t_data &data, double dt);


/****************************************
/// Spreading ring test
****************************************/
void spreading_ring_inner(t_data &data);
void spreading_ring_outer(t_data &data);


/****************************************
/// Damping
****************************************/
void damping(t_data &data, double dt);
bool initial_values_needed();
void copy_initial_values(t_data &data);
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


/****************************************
/// Legacy
****************************************/
void NonReflectingBoundary_inner(t_data &data, t_polargrid *VRadial,
				 t_polargrid *Density, t_polargrid *Energy);
void NonReflectingBoundary_outer(t_data &data, t_polargrid *VRadial,
				 t_polargrid *Density, t_polargrid *Energy);

void EvanescentBoundary(t_data &data, double step);

void ApplyOuterSourceMass(t_polargrid *Density, PolarGrid *VRadial);


} // namespace boundary_conditions