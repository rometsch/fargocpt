/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"

#include "../Theo.h"
#include "../find_cell_id.h"
#include "../global.h"
#include "../logging.h"
#include "../parameters.h"
#include "../util.h"
#include "../frame_of_reference.h"
#include "../simulation.h"
#include "../constants.h"
#include "../quantities.h"
#include "../axilib.h"
#include "../selfgravity.h"

#include <algorithm>
#include <cstring>
#include <cmath>
#include <vector>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include "viscosity/viscous_radial_speed.h"


// temporary
#include "SourceEuler.h"
#include "LowTasks.h"
#include "SideEuler.h"
extern std::vector<parameters::t_DampingType> parameters::damping_vector;

namespace boundary_conditions
{


// Determine whether initial values must be stored.
bool initial_values_needed() {
	const bool damping = parameters::is_damping_initial;
	const bool betacooling = parameters::cooling_beta_initial;
	const bool boundary = parameters::boundary_inner==parameters::boundary_condition_initial || parameters::boundary_outer==parameters::boundary_condition_initial;
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

	// Needed for OpenMP to work
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

void damping_initial_center_of_mass_outer(t_data &data, double dt)
{

	// use the correct radius array corresponding to quantity
	t_polargrid &vrad_arr = data[t_data::V_RADIAL];
    t_polargrid &vphi_arr = data[t_data::V_AZIMUTHAL];

    const unsigned int np = data.get_planetary_system().get_number_of_planets();
    const Pair com_pos = data.get_planetary_system().get_center_of_mass(np);
    const Pair com_vel =
	data.get_planetary_system().get_center_of_mass_velocity(np);
    const double com_mass = data.get_planetary_system().get_mass(np);

    // is this CPU in the outer damping domain?
    if ((parameters::damping_outer_limit < 1.0) &&
	(Rinf[vrad_arr.get_max_radial()] >
	 RMAX * parameters::damping_outer_limit)) {

	const unsigned int clamped_vrad_id = clamp_r_id_to_radii_grid(
		get_rinf_id(RMAX * parameters::damping_outer_limit), vrad_arr.is_vector())+1;

    const double tau = parameters::damping_time_factor * 2.0 * M_PI /
			 calculate_omega_kepler(parameters::damping_time_radius_outer);

	#pragma omp parallel for
	for (unsigned int n_radial = clamped_vrad_id;
		 n_radial < MaxMo_no_ghost_vr; ++n_radial) {
	    double factor = std::pow(
		(Rinf[n_radial] - RMAX * parameters::damping_outer_limit) /
		    (RMAX - RMAX * parameters::damping_outer_limit),
		2);
        const double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < vrad_arr.get_size_azimuthal(); ++n_azimuthal) {

		const double phi = (double)n_azimuthal * dphi;
		const double rinf = Rinf[n_radial];

		const double cell_x = rinf * std::cos(phi);
		const double cell_y = rinf * std::sin(phi);

		// Position in center of mass frame
		const double x_com = cell_x - com_pos.x;
		const double y_com = cell_y - com_pos.y;
		const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

		// pressure support correction
		double vr_init;
		double vphi_init;
		if (parameters::initialize_pure_keplerian) {
			vphi_init = compute_v_kepler(r_com, com_mass);
			if(parameters::initialize_vradial_zero){
				vr_init = 0.0;
			} else {
			vr_init = initial_viscous_radial_speed(r_com, com_mass);
			}
		} else {
            if(parameters::v_azimuthal_with_quadropole_support){
            vphi_init = initial_locally_isothermal_smoothed_v_az_with_quadropole_moment(r_com, com_mass);
            } else { // no quadropole support
            vphi_init = initial_locally_isothermal_smoothed_v_az(r_com, com_mass);
            }
			vr_init = viscous_speed::lookup_initial_vr_outer(r_com);
		}

		// Velocity in center of mass frame
		const double cell_vphi_com = vphi_init;
		const double cell_vr_com = vr_init;

		const double cell_vx_com =
		    (cell_vr_com * x_com - cell_vphi_com * y_com) / r_com;
		const double cell_vy_com =
		    (cell_vr_com * y_com + cell_vphi_com * x_com) / r_com;

		// shift velocity from center of mass frame to primary frame
		const double cell_vx = cell_vx_com + com_vel.x;
		const double cell_vy = cell_vy_com + com_vel.y;

		const double vr0 = (cell_x * cell_vx + cell_y * cell_vy) / rinf;

		const double vr = vrad_arr(n_radial, n_azimuthal);
		const double vr_new = (vr - vr0) * exp_factor + vr0;
		vrad_arr(n_radial, n_azimuthal) = vr_new;
	    }
	}

	const unsigned int clamped_vphi_id = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMAX * parameters::damping_outer_limit),
		vphi_arr.is_vector()) + 1;

	#pragma omp parallel for
	for (unsigned int n_radial = clamped_vphi_id;
		 n_radial < Max_no_ghost; ++n_radial) {
	    double factor = std::pow(
		(Rmed[n_radial] - RMAX * parameters::damping_outer_limit) /
		    (RMAX - RMAX * parameters::damping_outer_limit),
		2);
        const double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < vphi_arr.get_size_azimuthal(); ++n_azimuthal) {

		const double phi = ((double)n_azimuthal - 0.5) * dphi;
		const double rmed = Rmed[n_radial];

		const double cell_x = rmed * std::cos(phi);
		const double cell_y = rmed * std::sin(phi);

		// Position in center of mass frame
		const double x_com = cell_x - com_pos.x;
		const double y_com = cell_y - com_pos.y;
		const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

		// pressure support correction
		double vphi0;
		double vr0;
		if (parameters::initialize_pure_keplerian) {
			vphi0 = compute_v_kepler(r_com, com_mass);
			if(parameters::initialize_vradial_zero){
				vr0 = 0.0;
			} else {
			vr0 = initial_viscous_radial_speed(r_com, com_mass);
			}
		} else {
            if(parameters::v_azimuthal_with_quadropole_support){
            vphi0 = initial_locally_isothermal_smoothed_v_az_with_quadropole_moment(r_com, com_mass);
            } else { // no quadropole support
            vphi0 = initial_locally_isothermal_smoothed_v_az(r_com, com_mass);
            }
			vr0 = viscous_speed::lookup_initial_vr_outer(r_com);
		}

		// Velocity in center of mass frame
		const double cell_vphi_com = vphi0;
		const double cell_vr_com = vr0;

		const double cell_vx_com =
		    (cell_vr_com * x_com - cell_vphi_com * y_com) / r_com;
		const double cell_vy_com =
		    (cell_vr_com * y_com + cell_vphi_com * x_com) / r_com;

		// shift velocity from center of mass frame to primary frame
		const double cell_vx = cell_vx_com + com_vel.x;
		const double cell_vy = cell_vy_com + com_vel.y;

		const double vp0 = (cell_x * cell_vy - cell_vx * cell_y) / rmed -
				refframe::OmegaFrame * rmed;

		const double vp = vphi_arr(n_radial, n_azimuthal);
		const double vp_new = (vp - vp0) * exp_factor + vp0;
		vphi_arr(n_radial, n_azimuthal) = vp_new;
	    }
	}

	if (parameters::Adiabatic){
	t_polargrid &energy = data[t_data::ENERGY];
	//t_polargrid &sigma = data[t_data::DENSITY];
	#pragma omp parallel for
	for (unsigned int nr = clamped_vphi_id;
		 nr < Max_no_ghost; ++nr) {
		double factor = std::pow(
		(Rmed[nr] - RMAX * parameters::damping_outer_limit) /
			(RMAX - RMAX * parameters::damping_outer_limit),
		2);
		double exp_factor = std::exp(-dt * factor / tau);

		for (unsigned int naz = 0;
		 naz < energy.get_size_azimuthal(); ++naz) {
	const double cell_x = (*CellCenterX)(nr, naz);
	const double cell_y = (*CellCenterY)(nr, naz);

	// Position in center of mass frame
	const double x_com = cell_x - com_pos.x;
	const double y_com = cell_y - com_pos.y;
	const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

	/// Initial profile temperature
	const double cell_energy_profile = initial_energy(r_com, com_mass);

	const double cell_energy0 = cell_energy_profile;

	const double cell_energy = energy(nr, naz);
	const double energy_new = (cell_energy - cell_energy0) * exp_factor + cell_energy0;

	energy(nr, naz)  = energy_new;
		}
	}
	}
    }
}

void damping_initial_center_of_mass_inner(t_data &data, double dt)
{

	t_polargrid &vrad_arr = data[t_data::V_RADIAL];
	t_polargrid &vphi_arr = data[t_data::V_AZIMUTHAL];

	const unsigned int np = parameters::n_bodies_for_hydroframe_center;
	const Pair com_pos = data.get_planetary_system().get_center_of_mass(np);
	const Pair com_vel =
	data.get_planetary_system().get_center_of_mass_velocity(np);
	const double com_mass = data.get_planetary_system().get_mass(np);

	// is this CPU in the inner damping domain?
	if ((parameters::damping_inner_limit > 1.0) &&
	(Rmed[0] < RMIN * parameters::damping_inner_limit)) {

	const unsigned int clamped_vrad_id = clamp_r_id_to_radii_grid(
		get_rinf_id(RMIN * parameters::damping_inner_limit),
		vrad_arr.is_vector());

	double tau = parameters::damping_time_factor * 2.0 * M_PI /
			 calculate_omega_kepler(RMIN);

	#pragma omp parallel for
	for (unsigned int n_radial = One_no_ghost_vr;
		 n_radial <= clamped_vrad_id; ++n_radial) {
		double factor = std::pow(
		(Rinf[n_radial] - RMIN * parameters::damping_inner_limit) /
			(RMIN - RMIN * parameters::damping_inner_limit),
		2);
		const double exp_factor = std::exp(-dt * factor / tau);

		for (unsigned int n_azimuthal = 0;
		 n_azimuthal < vrad_arr.get_size_azimuthal(); ++n_azimuthal) {

		const double phi = (double)n_azimuthal * dphi;
		const double rinf = Rinf[n_radial];

		const double cell_x = rinf * std::cos(phi);
		const double cell_y = rinf * std::sin(phi);

		// Position in center of mass frame
		const double x_com = cell_x - com_pos.x;
		const double y_com = cell_y - com_pos.y;
		const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

		// pressure support correction
		double vr_init;
		double vphi_init;
		if (parameters::initialize_pure_keplerian) {
			vphi_init = compute_v_kepler(r_com, com_mass);
			if(parameters::initialize_vradial_zero){
				vr_init = 0.0;
			} else {
			vr_init = initial_viscous_radial_speed(r_com, com_mass);
			}
		} else {
			vphi_init = initial_locally_isothermal_smoothed_v_az(r_com, com_mass);
			vr_init = viscous_speed::lookup_initial_vr_inner(r_com);
		}

		// Velocity in center of mass frame
		const double cell_vphi_com = vphi_init;
		const double cell_vr_com = vr_init;

		const double cell_vx_com =
			(cell_vr_com * x_com - cell_vphi_com * y_com) / r_com;
		const double cell_vy_com =
			(cell_vr_com * y_com + cell_vphi_com * x_com) / r_com;

		// shift velocity from center of mass frame to primary frame
		const double cell_vx = cell_vx_com + com_vel.x;
		const double cell_vy = cell_vy_com + com_vel.y;

		const double vr0 = (cell_x * cell_vx + cell_y * cell_vy) / rinf;

		const double vr = vrad_arr(n_radial, n_azimuthal);
		const double vr_new = (vr - vr0) * exp_factor + vr0;
		vrad_arr(n_radial, n_azimuthal) = vr_new;
		}
	}

	const unsigned int clamped_vphi_id = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMIN * parameters::damping_inner_limit),
		vphi_arr.is_vector());

	#pragma omp parallel for
	for (unsigned int n_radial = Zero_no_ghost;
		 n_radial <= clamped_vphi_id; ++n_radial) {
		double factor = std::pow(
		(Rmed[n_radial] - RMAX * parameters::damping_outer_limit) /
			(RMAX - RMAX * parameters::damping_outer_limit),
		2);
		const double exp_factor = std::exp(-dt * factor / tau);

		for (unsigned int n_azimuthal = 0;
		 n_azimuthal < vphi_arr.get_size_azimuthal(); ++n_azimuthal) {

		const double phi = ((double)n_azimuthal - 0.5) * dphi;
		const double rmed = Rmed[n_radial];

		const double cell_x = rmed * std::cos(phi);
		const double cell_y = rmed * std::sin(phi);

		// Position in center of mass frame
		const double x_com = cell_x - com_pos.x;
		const double y_com = cell_y - com_pos.y;
		const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

		// pressure support correction
		double vphi0;
		double vr0;
		if (parameters::initialize_pure_keplerian) {
			vphi0 = compute_v_kepler(r_com, com_mass);
			if(parameters::initialize_vradial_zero){
				vr0 = 0.0;
			} else {
			vr0 = initial_viscous_radial_speed(r_com, com_mass);
			}
		} else {
			vphi0 = initial_locally_isothermal_smoothed_v_az(r_com, com_mass);
			vr0 = viscous_speed::lookup_initial_vr_inner(r_com);
		}

		// Velocity in center of mass frame
		const double cell_vphi_com = vphi0;
		const double cell_vr_com = vr0;

		const double cell_vx_com =
			(cell_vr_com * x_com - cell_vphi_com * y_com) / r_com;
		const double cell_vy_com =
			(cell_vr_com * y_com + cell_vphi_com * x_com) / r_com;

		// shift velocity from center of mass frame to primary frame
		const double cell_vx = cell_vx_com + com_vel.x;
		const double cell_vy = cell_vy_com + com_vel.y;

		const double vp0 = (cell_x * cell_vy - cell_vx * cell_y) / rmed -
				refframe::OmegaFrame * rmed;

		const double vp = vphi_arr(n_radial, n_azimuthal);
		const double vp_new = (vp - vp0) * exp_factor + vp0;
		vphi_arr(n_radial, n_azimuthal) = vp_new;
		}
	}

	if (parameters::Adiabatic){
	t_polargrid &energy = data[t_data::ENERGY];
	//t_polargrid &sigma = data[t_data::DENSITY];
	#pragma omp parallel for
	for (unsigned int nr = Zero_no_ghost;
		 nr <= clamped_vphi_id; ++nr) {
		double factor = std::pow(
		(Rmed[nr] - RMAX * parameters::damping_outer_limit) /
			(RMAX - RMAX * parameters::damping_outer_limit),
		2);
		const double exp_factor = std::exp(-dt * factor / tau);

		for (unsigned int naz = 0;
		 naz < energy.get_size_azimuthal(); ++naz) {
	const double cell_x = (*CellCenterX)(nr, naz);
	const double cell_y = (*CellCenterY)(nr, naz);

	// Position in center of mass frame
	const double x_com = cell_x - com_pos.x;
	const double y_com = cell_y - com_pos.y;
	const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

	/// Initial profile temperature
	const double cell_energy_profile = initial_energy(r_com, com_mass);

	const double cell_energy0 = cell_energy_profile;

	const double cell_energy = energy(nr, naz);
	const double energy_new = (cell_energy - cell_energy0) * exp_factor + cell_energy0;

	energy(nr, naz)  = energy_new;
		}
	}
	}

	}
}

/**
	damping of all selected quantities
*/
void damping(t_data &data, double dt)
{
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
