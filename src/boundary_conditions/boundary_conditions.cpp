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




void apply_boundary_condition(t_data &data, const double current_time, const double dt, const bool final) {
	// if this is the final call of the boundary condition function during one timestep and damping is enable, do it
    if (final && parameters::damping) {
		damping(data, dt);

		if(ECC_GROWTH_MONITOR){
			quantities::calculate_disk_delta_ecc_peri(data, delta_ecc_damp, delta_peri_damp);
		}
    }

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
		custom(data);
	}

}

void old_apply_boundary_condition(t_data &data, const double current_time, const double dt, const bool final)
{
    // if this is the final call of the boundary condition function during one timestep and damping is enable, do it
    if (final && parameters::damping) {
		damping(data, dt);

		if(ECC_GROWTH_MONITOR){
			quantities::calculate_disk_delta_ecc_peri(data, delta_ecc_damp, delta_peri_damp);
		}
    }

    // inner boundary
    switch (parameters::boundary_inner) {
		case parameters::boundary_condition_reference:
			reference_value_boundary_inner(data);
		break;
		case parameters::boundary_condition_open:
			open_boundary_inner(data);
		break;
		case parameters::boundary_condition_reflecting:
			reflecting_boundary_inner(data);
		break;
		case parameters::boundary_condition_zero_gradient:
			zero_gradient_boundary_inner(data);
		break;
		case parameters::boundary_condition_boundary_layer:
			boundary_layer_inner_boundary(data);
		break;
		case parameters::boundary_condition_nonreflecting:
			NonReflectingBoundary_inner(data, &data[t_data::V_RADIAL],
						&data[t_data::SIGMA],
						&data[t_data::ENERGY]);
		break;
		case parameters::boundary_condition_viscous_outflow:
			viscous_outflow_boundary_inner(data);
		break;
		case parameters::boundary_condition_center_of_mass_initial:
			if (final) {
				damping_initial_center_of_mass_inner(data, dt);
			}
			initial_center_of_mass_boundary_inner(data);
		break;
		case parameters::boundary_condition_jibin_spreading_ring:
			spreading_ring_inner(data);
		break;
		case parameters::boundary_condition_evanescent: // evanescent works only for
								// inner and outer together
								// until now
			if (parameters::boundary_outer ==
				parameters::boundary_condition_evanescent) {
				EvanescentBoundary(data, dt);
			} else {
				logging::print_master(
				LOG_ERROR
				"Different EvanescentBoundary Parameters. Old parameter file?\n");
				die("inner/outer evanescent boundary not implemented yet");
			}
		break;
		case parameters::boundary_condition_keplerian:
			keplerian2d_boundary_inner(data);
		break;
		case parameters::boundary_condition_precribed_time_variable:
			die("Inner precribed time variable boundary condition is not implemented yet!\n");
		break;
    }

    // outer boundary
    switch (parameters::boundary_outer) {
		case parameters::boundary_condition_reference:
			reference_value_boundary_outer(data);
		break;
		case parameters::boundary_condition_open:
			open_boundary_outer(data);
		break;
		case parameters::boundary_condition_reflecting:
			reflecting_boundary_outer(data);
		break;
		case parameters::boundary_condition_center_of_mass_initial: {
			if (final) {
				damping_initial_center_of_mass_outer(data, dt);
			}
			initial_center_of_mass_boundary_outer(data);
		break;
		}
		case parameters::boundary_condition_zero_gradient:
			zero_gradient_boundary_outer(data);
		break;
		case parameters::boundary_condition_boundary_layer:
			boundary_layer_outer_boundary(data);
		break;
		case parameters::boundary_condition_nonreflecting:
			NonReflectingBoundary_outer(data, &data[t_data::V_RADIAL],
							&data[t_data::SIGMA],
							&data[t_data::ENERGY]);
		break;
		case parameters::boundary_condition_jibin_spreading_ring:
			spreading_ring_outer(data);
		break;
		case parameters::boundary_condition_precribed_time_variable: {
			boundary_condition_precribed_time_variable_outer(data,
									&data[t_data::SIGMA], current_time);
		} break;
		case parameters::boundary_condition_viscous_outflow:
			die("outer viscous outflow boundary not implemented");
		break;
		case parameters::boundary_condition_evanescent:
			// EvanescentBoundary_outer(VRadial, VAzimuthal, Density, Energy, dt,
			// sys);
		break;
		case parameters::boundary_condition_keplerian:
			keplerian2d_boundary_outer(data);
		break;
    }

    /// d Omega / dr = 0 has really bad effects on the massflow test
    /// not recommended
    if (parameters::domegadr_zero) {
		zero_shear_boundary(data[t_data::V_RADIAL]);
    }

    if (OuterSourceMass) {
		ApplyOuterSourceMass(&data[t_data::SIGMA], &data[t_data::V_RADIAL]);
    }

    if (parameters::massoverflow) {
		boundary_conditions::mass_overflow_willy(data, nullptr, false);
    }
}

} // namespace boundary_conditions
