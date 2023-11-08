#pragma once

#include "../data.h"

namespace boundary_conditions
{

/****************************************
/// Types of boundary conditions.
*****************************************/


extern void (*sigma_inner_func)(t_polargrid &, t_polargrid &, t_data &d);
extern void (*sigma_outer_func)(t_polargrid &, t_polargrid &, t_data &d);
extern void (*energy_inner_func)(t_polargrid &, t_polargrid &, t_data &d);
extern void (*energy_outer_func)(t_polargrid &, t_polargrid &, t_data &d);
extern void (*vrad_inner_func)(t_polargrid &, t_polargrid &, t_data &d);
extern void (*vrad_outer_func)(t_polargrid &, t_polargrid &, t_data &d);
extern void (*vaz_inner_func)(t_polargrid &, t_polargrid &, t_data &d);
extern void (*vaz_outer_func)(t_polargrid &, t_polargrid &, t_data &d);
extern void (*special_inner_func)(t_polargrid &, t_polargrid &, t_data &d);
extern void (*special_outer_func)(t_polargrid &, t_polargrid &, t_data &d);


extern std::string sigma_inner_name;
extern std::string sigma_outer_name;
extern std::string energy_inner_name;
extern std::string energy_outer_name;
extern std::string vrad_inner_name;
extern std::string vrad_outer_name;
extern std::string vaz_inner_name;
extern std::string vaz_outer_name;
extern std::string composite_inner_name;
extern std::string composite_outer_name;

extern double keplerian_azimuthal_outer_factor;
extern double keplerian_azimuthal_inner_factor;
extern double keplerian_radial_outer_factor;
extern double keplerian_radial_inner_factor;
	

/****************************************
/// Parameters
****************************************/

// speed of the viscous boundary inflow
extern double viscous_outflow_speed;

extern bool rochlobe_overflow;
// Rochlobe overflow: number of Nbody object from which the stream originates
extern unsigned int rof_planet;
// Rochlobe overflow, sets initial temperature and width of stream
extern double rof_temperature;
// Rochlobe overflow accretion rate through stream
extern double rof_mdot;
// Rochlobe overflow ramping time
extern double rof_rampingtime;
// Rochlobe overflow variable transfer
extern bool rof_variableTransfer;
// Rochlobe overflow variable transfer averaging time
extern double rof_averaging_time;
// Rochlobe overflow variable transfer gamma
extern double rof_gamma;

/// enable different damping types
enum t_damping_type {
    damping_none,
    damping_reference,
    damping_mean,
    damping_zero,
    damping_visc
};

/// Struct for handling damping at boundaries
struct t_DampingType {
    void (*inner_damping_function)(t_polargrid &, t_polargrid &, double);
    void (*outer_damping_function)(t_polargrid &, t_polargrid &, double);
    t_data::t_polargrid_type array_to_damp;
    t_data::t_polargrid_type array_with_damping_values;
    t_damping_type type_inner;
    t_damping_type type_outer;
};
extern int damping_energy_id;

extern bool damping_enabled;
/// is at least one variable damped to initial values
extern bool is_damping_reference;
/// inner damping limit
extern double damping_inner_limit;
/// outer damping limit
extern double damping_outer_limit;
/// damping time factor
extern double damping_time_factor;
extern double damping_time_radius_outer;
/// vector to handle damping structs
extern std::vector<t_DampingType> damping_vector;


/****************************************
/// Functions
****************************************/

void parse_config();

void apply_boundary_condition(t_data &data, const double current_time, const double dt, const bool final);

bool reference_values_needed();
void init(t_data &data);

/****************************************
/// Individual variable boundaries
****************************************/

void zero_gradient_inner(t_polargrid &x, t_polargrid &dummy, t_data &ddummy);
void zero_gradient_outer(t_polargrid &x, t_polargrid &dummy, t_data &ddummy);

void outflow_inner(t_polargrid &x, t_polargrid &dummy, t_data &ddummy);
void outflow_outer(t_polargrid &x, t_polargrid &dummy, t_data &ddummy);

void reference_inner(t_polargrid &x, t_polargrid &x0, t_data &ddummy);
void reference_outer(t_polargrid &x, t_polargrid &x0, t_data &ddummy);

void viscous_outflow_inner(t_polargrid &vrad, t_polargrid &dummy, t_data &data);
void viscous_inflow_outer(t_polargrid &vrad, t_polargrid &dummy, t_data &data);

void zero_shear_inner(t_polargrid &vaz, t_polargrid &dummy, t_data &ddummy);
void zero_shear_outer(t_polargrid &vaz, t_polargrid &dummy, t_data &ddummy);

void reflecting_inner(t_polargrid &vrad, t_polargrid &dummy, t_data &ddummy);
void reflecting_outer(t_polargrid &vrad, t_polargrid &dummy, t_data &ddummy);

void keplerian_azimuthal_inner(t_polargrid &vaz, t_polargrid &dummy, t_data &ddummy);
void keplerian_azimuthal_outer(t_polargrid &vaz, t_polargrid &dummy, t_data &ddummy);

void keplerian_radial_inner(t_polargrid &vrad, t_polargrid &dummy, t_data &ddummy);
void keplerian_radial_outer(t_polargrid &vrad, t_polargrid &dummy, t_data &ddummy);

void balanced_inner(t_polargrid &vaz, t_polargrid &dummy, t_data &ddummy);
void balanced_outer(t_polargrid &vaz, t_polargrid &dummy, t_data &ddummy);

void diskmodel_inner_sigma(t_polargrid &x, t_polargrid &dummy, t_data &ddummy);
void diskmodel_outer_sigma(t_polargrid &x, t_polargrid &dummy, t_data &ddummy);
void diskmodel_inner_energy(t_polargrid &x, t_polargrid &dummy, t_data &ddummy);
void diskmodel_outer_energy(t_polargrid &x, t_polargrid &dummy, t_data &ddummy);


/****************************************
/// No operation for using custom 
/// or special boundaries
****************************************/
void no_operation(t_polargrid &x, t_polargrid &dummy, t_data &ddummy);

/****************************************
/// Custom boundary conditions
****************************************/
// Modify these to implement you own custom boundary conditions.

void custom_inner(t_data &data, const double t, const double dt);
void custom_outer(t_data &data, const double t, const double dt);
void init_custom(t_data& data);
void cleanup_custom();

/****************************************
/// Inflow boundaries
****************************************/
void rochelobe_overflow_boundary(t_data &data, t_polargrid *densitystar,
			 bool transport);

/****************************************
/// BC centered on center of mass
****************************************/

void init_center_of_mass(t_data &data);
// void diskmodel_center_of_mass_inner(t_data &data, const double dt, const bool final);
// void diskmodel_center_of_mass_outer(t_data &data, const double dt, const bool final);

void diskmodel_center_of_mass_boundary_inner(t_data &data);
void damping_diskmodel_center_of_mass_inner(t_data &data, double dt);
void diskmodel_center_of_mass_boundary_outer(t_data &data);
void damping_diskmodel_center_of_mass_outer(t_data &data, double dt);

/****************************************
/// Damping
****************************************/
void damping(t_data &data, const double dt);
void damping_config();
void describe_damping();
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



} // namespace boundary_conditions