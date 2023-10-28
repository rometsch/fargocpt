/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"

#include "../global.h"
#include "../logging.h"
#include "../parameters.h"
#include "../quantities.h"

extern boolean OuterSourceMass;

namespace boundary_conditions
{

void init(t_data &data) {
	init_custom(data);

	if (composite_inner_name == "centerofmass" || composite_outer_name == "centerofmass") {
		init_center_of_mass(data);
	}
}

bool reference_values_needed() {
		bool needed = false;
		needed = needed || sigma_inner_name == "reference";
		needed = needed || sigma_outer_name == "reference";
		needed = needed || energy_inner_name == "reference";
		needed = needed || energy_outer_name == "reference";
		needed = needed || vrad_inner_name == "reference";
		needed = needed || vrad_outer_name == "reference";
		needed = needed || vaz_inner_name == "reference";
		needed = needed || vaz_outer_name == "reference";
		return needed;
}


static void handle_damping(t_data &data, const double dt, const bool final) {
	// if this is the final call of the boundary condition function during one timestep and damping is enable, do it
    if (final && damping_enabled) {
		damping(data, dt);

		if(ECC_GROWTH_MONITOR){
			quantities::calculate_disk_delta_ecc_peri(data, delta_ecc_damp, delta_peri_damp);
		}
    }
}

void apply_boundary_condition(t_data &data, const double current_time, const double dt, const bool final) {
	
	handle_damping(data, dt, final);

	{
		t_polargrid &x = data[t_data::SIGMA];
		t_polargrid &x0 = data[t_data::SIGMA0];
		sigma_inner_func(x, x0, data);
		sigma_outer_func(x, x0, data);
	}

	{
		t_polargrid &x = data[t_data::ENERGY];
		t_polargrid &x0 = data[t_data::ENERGY0];
		energy_inner_func(x, x0, data);
		energy_outer_func(x, x0, data);
	}

	{
		t_polargrid &x = data[t_data::V_RADIAL];
		t_polargrid &x0 = data[t_data::V_RADIAL0];
		vrad_inner_func(x, x0, data);
		vrad_outer_func(x, x0, data);
	}

	{
		t_polargrid &x = data[t_data::V_AZIMUTHAL];
		t_polargrid &x0 = data[t_data::V_AZIMUTHAL0];
		vaz_inner_func(x, x0, data);
		vaz_outer_func(x, x0, data);
	}

	if (special_name == "custom") {
		custom(data, current_time, dt);
	}

}


} // namespace boundary_conditions
