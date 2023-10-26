/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"
#include "../parameters.h"
#include "../find_cell_id.h"
#include "../global.h"
#include "../Theo.h"
#include "../util.h"
#include "../SourceEuler.h"

#include <vector>
#include <cassert>
extern std::vector<parameters::t_DampingType> parameters::damping_vector;

namespace boundary_conditions
{


// Determine whether initial values must be stored.
bool initial_values_needed() {
	const bool damping = parameters::is_damping_reference;
	const bool betacooling = parameters::cooling_beta_reference;
	const bool boundary = boundary_conditions::reference_values_needed();
	return damping || betacooling || boundary;
}


void copy_initial_values(t_data &data) {
    if (initial_values_needed()) {
	// save starting values (needed for damping)
		copy_polargrid(data[t_data::V_RADIAL0], data[t_data::V_RADIAL]);
		copy_polargrid(data[t_data::V_AZIMUTHAL0], data[t_data::V_AZIMUTHAL]);
		copy_polargrid(data[t_data::SIGMA0], data[t_data::SIGMA]);
		copy_polargrid(data[t_data::ENERGY0], data[t_data::ENERGY]);
    }
}

void damping_single_inner(t_polargrid &quantity, t_polargrid &quantity0,
			  double dt)
{
    // use the correct radius array corresponding to quantity
    t_radialarray &radius = quantity.is_scalar() ? Rb : Ra;
    const bool is_density = strcmp(quantity.get_name(), "Sigma") == 0;

    // is this CPU in the inner damping domain?
    if ((parameters::damping_inner_limit > 1.0) &&
	(radius[0] < RMIN * parameters::damping_inner_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMIN * parameters::damping_inner_limit), false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMIN * parameters::damping_inner_limit), true);
	}

    const double tau = parameters::damping_time_factor * 2.0 * M_PI /
		     calculate_omega_kepler(RMIN);

	// Needed for OpenMP to work
	double &InnerWaveDampingMassCreation = MassDelta.InnerWaveDampingMassCreation;
	double &InnerWaveDampingMassRemoval = MassDelta.InnerWaveDampingMassRemoval;

	#pragma omp parallel for reduction (+ : InnerWaveDampingMassCreation, InnerWaveDampingMassRemoval)
	for (unsigned int n_radial = 0; n_radial <= limit; ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMIN * parameters::damping_inner_limit) /
		    (RMIN - RMIN * parameters::damping_inner_limit),
		2);
	    double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
        const double X = quantity(n_radial, n_azimuthal);
        const double X0 = quantity0(n_radial, n_azimuthal);
        const double Xnew = (X - X0) * exp_factor + X0;
		quantity(n_radial, n_azimuthal) = Xnew;
        const double delta = Xnew - X;
		if (is_density) {
		    if (delta > 0) {
			sum_without_ghost_cells(InnerWaveDampingMassCreation,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(InnerWaveDampingMassRemoval,
						-delta * Surf[n_radial],
						n_radial);
		    }
		}
	    }
	}
    }
}

void damping_single_outer(t_polargrid &quantity, t_polargrid &quantity0,
			  double dt)
{
    // use the correct radius array corresponding to quantity
    t_radialarray &radius = quantity.is_scalar() ? Rb : Ra;
    const bool is_density = strcmp(quantity.get_name(), "Sigma") == 0;

    // is this CPU in the outer damping domain?
    if ((parameters::damping_outer_limit < 1.0) &&
	(radius[quantity.get_max_radial()] >
	 RMAX * parameters::damping_outer_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMAX * parameters::damping_outer_limit) + 1, false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMAX * parameters::damping_outer_limit) + 1, true);
	}

	double tau = parameters::damping_time_factor * 2.0 * M_PI /
			 calculate_omega_kepler(parameters::damping_time_radius_outer);

	// Needed for OpenMP to work
	double &OuterWaveDampingMassCreation = MassDelta.OuterWaveDampingMassCreation;
	double &OuterWaveDampingMassRemoval = MassDelta.OuterWaveDampingMassRemoval;

	#pragma omp parallel for reduction (+ : OuterWaveDampingMassCreation, OuterWaveDampingMassRemoval)
	for (unsigned int n_radial = limit;
	     n_radial < quantity.get_size_radial(); ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMAX * parameters::damping_outer_limit) /
		    (RMAX - RMAX * parameters::damping_outer_limit),
		2);
	    double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
        const double X = quantity(n_radial, n_azimuthal);
        const double X0 = quantity0(n_radial, n_azimuthal);
        const double Xnew = (X - X0) * exp_factor + X0;
		quantity(n_radial, n_azimuthal) = Xnew;
        const double delta = Xnew - X;
		if (is_density) {
		    if (delta > 0) {
			sum_without_ghost_cells(OuterWaveDampingMassCreation,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(OuterWaveDampingMassRemoval,
						-delta * Surf[n_radial],
						n_radial);
		    }
		}
	    }
	}
    }
}

void damping_single_inner_zero(t_polargrid &quantity, t_polargrid &quantity0,
			       double dt)
{
    (void)quantity0;
    // use the correct radius array corresponding to quantity
    t_radialarray &radius = quantity.is_scalar() ? Rb : Ra;
    const bool is_density = strcmp(quantity.get_name(), "Sigma") == 0;

    // is this CPU in the inner damping domain?
    if ((parameters::damping_inner_limit > 1.0) &&
	(radius[0] < RMIN * parameters::damping_inner_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMIN * parameters::damping_inner_limit), false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMIN * parameters::damping_inner_limit), true);
	}

    const double tau = parameters::damping_time_factor * 2.0 * M_PI /
		     calculate_omega_kepler(RMIN);

	// Needed for OpenMP to work
	double &InnerWaveDampingMassCreation = MassDelta.InnerWaveDampingMassCreation;
	double &InnerWaveDampingMassRemoval = MassDelta.InnerWaveDampingMassRemoval;

	#pragma omp parallel for reduction (+ : InnerWaveDampingMassCreation, InnerWaveDampingMassRemoval)
	for (unsigned int n_radial = 0; n_radial <= limit; ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMIN * parameters::damping_inner_limit) /
		    (RMIN - RMIN * parameters::damping_inner_limit),
		2);
	    double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
		const double X = quantity(n_radial, n_azimuthal);
		double X0;
		if (is_density) {
		    X0 = parameters::sigma_floor * parameters::sigma0;
		} else {
		    X0 = 0.0;
		}
		const double Xnew = (X - X0) * exp_factor + X0;
		quantity(n_radial, n_azimuthal) = Xnew;
		const double delta = Xnew - X;
		if (is_density) {
		    if (delta > 0) {
			sum_without_ghost_cells(InnerWaveDampingMassCreation,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(InnerWaveDampingMassRemoval,
						-delta * Surf[n_radial],
						n_radial);
		    }
		}
	    }
	}
    }
}

void damping_single_outer_zero(t_polargrid &quantity, t_polargrid &quantity0,
			       double dt)
{
    (void)quantity0;
    // use the correct radius array corresponding to quantity
    t_radialarray &radius = quantity.is_scalar() ? Rb : Ra;

    const bool is_density = strcmp(quantity.get_name(), "Sigma") == 0;

    // is this CPU in the outer damping domain?
    if ((parameters::damping_outer_limit < 1.0) &&
	(radius[quantity.get_max_radial()] >
	 RMAX * parameters::damping_outer_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMAX * parameters::damping_outer_limit) + 1, false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMAX * parameters::damping_outer_limit) + 1, true);
	}

    const double tau = parameters::damping_time_factor * 2.0 * M_PI /
			 calculate_omega_kepler(parameters::damping_time_radius_outer);

	// Needed for OpenMP to work
	double &OuterWaveDampingMassCreation = MassDelta.OuterWaveDampingMassCreation;
	double &OuterWaveDampingMassRemoval = MassDelta.OuterWaveDampingMassRemoval;

	#pragma omp parallel for reduction (+ : OuterWaveDampingMassCreation, OuterWaveDampingMassRemoval)
	for (unsigned int n_radial = limit;
	     n_radial < quantity.get_size_radial(); ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMAX * parameters::damping_outer_limit) /
		    (RMAX - RMAX * parameters::damping_outer_limit),
		2);
	    double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
		const double X = quantity(n_radial, n_azimuthal);
		const double X0 = 0.0;
		const double Xnew = (X - X0) * exp_factor + X0;
		quantity(n_radial, n_azimuthal) = Xnew;
		const double delta = Xnew - X;
		if (is_density) {
		    if (delta > 0) {
			sum_without_ghost_cells(OuterWaveDampingMassCreation,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(OuterWaveDampingMassRemoval,
						-delta * Surf[n_radial],
						n_radial);
		    }
		}
	    }
	}
    }
}

void damping_single_inner_mean(t_polargrid &quantity, t_polargrid &quantity0,
			       double dt)
{
    // use the correct radius array corresponding to quantity
    t_radialarray &radius = quantity.is_scalar() ? Rb : Ra;
    const bool is_density = strcmp(quantity.get_name(), "Sigma") == 0;

    // is this CPU in the inner damping domain?
    if ((parameters::damping_inner_limit > 1.0) &&
	(radius[0] < RMIN * parameters::damping_inner_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMIN * parameters::damping_inner_limit), false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMIN * parameters::damping_inner_limit), true);
	}

    const double tau = parameters::damping_time_factor * 2.0 * M_PI /
		     calculate_omega_kepler(RMIN);

	// get mean quantity
	for (unsigned int n_radial = 0; n_radial <= limit; ++n_radial) {
	    quantity0(n_radial, 0) = 0.0;
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
		quantity0(n_radial, 0) += quantity(n_radial, n_azimuthal);
	    }
	    quantity0(n_radial, 0) /= quantity.get_size_azimuthal();
	}

	// Needed for OpenMP to work, OpenMP want a local variable for reduction
	double &InnerWaveDampingMassCreation = MassDelta.InnerWaveDampingMassCreation;
	double &InnerWaveDampingMassRemoval = MassDelta.InnerWaveDampingMassRemoval;

	#pragma omp parallel for reduction (+ : InnerWaveDampingMassCreation, InnerWaveDampingMassRemoval)
	for (unsigned int n_radial = 0; n_radial <= limit; ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMIN * parameters::damping_inner_limit) /
		    (RMIN - RMIN * parameters::damping_inner_limit),
		2);
        const double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
        const double X = quantity(n_radial, n_azimuthal);
        const double X0 = quantity0(n_radial, 0);
        const double Xnew = (X - X0) * exp_factor + X0;
		quantity(n_radial, n_azimuthal) = Xnew;
        const double delta = Xnew - X;
		if (is_density) {
		    if (delta > 0) {
			sum_without_ghost_cells(InnerWaveDampingMassCreation,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(InnerWaveDampingMassRemoval,
						-delta * Surf[n_radial],
						n_radial);
		    }
		}
	    }
	}
    }
}

void damping_vradial_inner_visc(t_polargrid &vrad, t_polargrid &viscosity,
				double dt)
{
    bool is_vrad = strcmp(vrad.get_name(), "vrad") == 0;
    assert(is_vrad);
    (void)is_vrad; /// removes IDE 'unused' warning

    t_radialarray &rinf = Ra;

    // is this CPU in the inner damping domain?
    if ((parameters::damping_inner_limit > 1.0) &&
	(rinf[0] < RMIN * parameters::damping_inner_limit)) {
	// find range

	const static unsigned int limit =
	    get_rinf_id(RMIN * parameters::damping_inner_limit);

	const double tau = parameters::damping_time_factor * 2.0 * M_PI /
			   calculate_omega_kepler(RMIN);

	const double s = parameters::viscous_outflow_speed;

	#pragma omp parallel for
	for (unsigned int n_radial = Zero_no_ghost; n_radial <= limit;
	     ++n_radial) {
	    double factor = std::pow(
		(rinf[n_radial] - RMIN * parameters::damping_inner_limit) /
		    (RMIN - RMIN * parameters::damping_inner_limit),
		2);
        const double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < vrad.get_size_azimuthal(); ++n_azimuthal) {
		const double viscosity_above = viscosity(n_radial, n_azimuthal);

		double Nu;

		if (n_radial == 0) {
		    Nu = viscosity_above;
		} else {
		    const double viscosity_below =
			viscosity(n_radial - 1, n_azimuthal);
		    Nu = 0.5 * (viscosity_above + viscosity_below);
		}

		// V_rad =  - 1.5 / r * Nu (Kley, Papaloizou and Ogilvie, 2008)
		const double V_rad = -1.5 * s / Rinf[n_radial] * Nu;

		const double X = vrad(n_radial, n_azimuthal);
		const double X0 = V_rad;
		const double Xnew = (X - X0) * exp_factor + X0;
		vrad(n_radial, n_azimuthal) = Xnew;
	    }
	}
    }
}

void damping_single_outer_mean(t_polargrid &quantity, t_polargrid &quantity0,
			       double dt)
{
    // use the correct radius array corresponding to quantity
    t_radialarray &radius = quantity.is_scalar() ? Rb : Ra;
    bool is_density = strcmp(quantity.get_name(), "Sigma") == 0;

    // is this CPU in the outer damping domain?
    if ((parameters::damping_outer_limit < 1.0) &&
	(radius[quantity.get_max_radial()] >
	 RMAX * parameters::damping_outer_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMAX * parameters::damping_outer_limit) + 1, false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMAX * parameters::damping_outer_limit) + 1, true);
	}

    const double tau = parameters::damping_time_factor * 2.0 * M_PI /
			 calculate_omega_kepler(parameters::damping_time_radius_outer);

	for (unsigned int n_radial = limit;
	     n_radial < quantity.get_size_radial(); ++n_radial) {
	    quantity0(n_radial, 0) = 0.0;
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
		quantity0(n_radial, 0) += quantity(n_radial, n_azimuthal);
	    }
	    quantity0(n_radial, 0) /= quantity.get_size_azimuthal();
	}

	// Needed for OpenMP to work
	double &OuterWaveDampingMassCreation = MassDelta.OuterWaveDampingMassCreation;
	double &OuterWaveDampingMassRemoval = MassDelta.OuterWaveDampingMassRemoval;

	#pragma omp parallel for reduction (+ : OuterWaveDampingMassCreation, OuterWaveDampingMassRemoval)
	for (unsigned int n_radial = limit;
	     n_radial < quantity.get_size_radial(); ++n_radial) {
        const double factor = std::pow(
		(radius[n_radial] - RMAX * parameters::damping_outer_limit) /
		    (RMAX - RMAX * parameters::damping_outer_limit),
		2);
        const double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
        const double X = quantity(n_radial, n_azimuthal);
        const double X0 = quantity0(n_radial, 0);
        const double Xnew = (X - X0) * exp_factor + X0;
		quantity(n_radial, n_azimuthal) = Xnew;
        const double delta = Xnew - X;
		if (is_density) {
		    if (delta > 0) {
			sum_without_ghost_cells(OuterWaveDampingMassCreation,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(OuterWaveDampingMassRemoval,
						-delta * Surf[n_radial],
						n_radial);
		    }
		}
	    }
	}
    }
}


/**
	damping of all selected quantities
*/
void damping(t_data &data, const double dt)
{
	// Iterate over all function pointers previously stored
    if (parameters::damping) {
	const unsigned int number_of_quantities_to_damp =
	    parameters::damping_vector.size();

	for (unsigned int i = 0; i < number_of_quantities_to_damp; ++i) {
	    const parameters::t_DampingType *damper =
		&parameters::damping_vector[i];
	    if (damper->inner_damping_function != nullptr)
		(damper->inner_damping_function)(
		    data[damper->array_to_damp],
		    data[damper->array_with_damping_values], dt);
	    if (damper->outer_damping_function != nullptr)
		(damper->outer_damping_function)(
		    data[damper->array_to_damp],
		    data[damper->array_with_damping_values], dt);
	}
    }
}

} // namespace boundary_conditions
