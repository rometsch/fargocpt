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
extern boolean OuterSourceMass;

namespace boundary_conditions
{

void apply_boundary_condition(t_data &data, const double current_time, const double dt, const bool final)
{
    // if this is the final boundary condition and damping is enable, do it
    if (final && parameters::damping) {
	damping(data, dt);

	if(ECC_GROWTH_MONITOR){
		quantities::calculate_disk_delta_ecc_peri(data, delta_ecc_damp, delta_peri_damp);
	}
    }

    // inner boundary
    switch (parameters::boundary_inner) {
	case parameters::boundary_condition_initial:
	initial_boundary_inner(data);
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
	case parameters::boundary_condition_initial:
	initial_boundary_outer(data);
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

void mass_overflow(t_data &data, const double current_time)
{
    if (CPU_Rank != CPU_Highest) {
	return;
    }

    static double last_PhysicalTime = 0.0;

    double dt = current_time - last_PhysicalTime;

    if (dt == 0.0) {
	return;
    }

    last_PhysicalTime = current_time;

    // get location of binary star
    if (parameters::mof_planet + 1 >
	data.get_planetary_system().get_number_of_planets()) {
	logging::print_master(
	    LOG_ERROR
	    "Wrong Planet/Star for Mass Overflow specified! Old parameter file?\n");
	die("Wrong Planet/Star for Mass Overflow specified! Old parameter file?");
    }
    const t_planet &planet =
	data.get_planetary_system().get_planet(parameters::mof_planet);
    const double xplanet = planet.get_x();
    const double yplanet = planet.get_y();

    const double omega_planet = planet.get_omega();
    // get grid cell where binary star is nearest
    // atan2(y,x) is from -PI to PI
    double angle = std::atan2(yplanet, xplanet) * 0.5 * M_1_PI;
    if (angle < 0) {
	angle += 1.0;
    }

    const unsigned int Nrad = data[t_data::SIGMA].get_max_radial();
    const unsigned int Nphi = data[t_data::SIGMA].get_size_azimuthal();
    const double r_cell = Rmed[Nrad];

    const double vr_fraction = 0.002;
    const double vr_stream =
	-omega_planet * r_cell *
	vr_fraction; // radial inflow velocity is a fraction of v_Kepler, (e.g.
		     // vr_fraction = 0.002)
    const double vphi_stream =
	(omega_planet - refframe::OmegaFrame) *
	r_cell; // set angular velocity to zero (in co-rotating frame)

	const double mdot_code = parameters::mof_value;

    const double Sigma_stream = mdot_code * dt / Surf[Nrad - 2];

    int nearest_grid_cell =
	(int)((double)Nphi * angle + 0.5); // nearest gridcell to binary star
    nearest_grid_cell = nearest_grid_cell % Nphi;

    // Calculate sigma from temperature according to
    // Meyer & Meyer-Hofmeister (1983a) equation 17
    const double Porb =
	2.0 * M_PI / omega_planet * units::time.get_cgs_factor() / 3600.0; // TODO: check unit!
    // cross section Q
    const double Q = 2.4e13 * parameters::mof_temperature * Porb * Porb;
    // stream radius W
    const double W = std::sqrt(Q / M_PI);
    // circumference circ
    const double circ = 2.0 * M_PI * r_cell * units::length.get_cgs_factor(); // TODO: check unit!

    double sigma = 2 * W / circ;
    int number_of_cells =
	(double)Nphi * 3.0 *
	sigma; // walk through 3 sigma, so we get 99.7% accuracy
    double sigmabar = Nphi * sigma;

    double check = 0.0;
    for (int i = -number_of_cells; i <= number_of_cells; i++) {

	const double t_ramp =
	    parameters::mof_rampingtime * planet.get_orbital_period();

	double ramp_factor;
	if (current_time < t_ramp) {
	    ramp_factor = std::pow(std::sin(current_time * M_PI_2 / t_ramp), 4);
	} else {
	    ramp_factor = 1.0;
	}

	// adapt gauss profile
	double weight_factor;
	int gridcell;
	if (number_of_cells == 0) {
	    gridcell = nearest_grid_cell;
	    weight_factor = 1.0;
	} else {
	    gridcell = (nearest_grid_cell + i + Nphi) % Nphi;
	    weight_factor = 1.0 / (sigmabar * std::sqrt(2.0 * M_PI)) *
			    std::exp(-1.0 / 2.0 * std::pow(i / sigmabar, 2));
	}
	check += weight_factor;

	int gridcell_r = gridcell + 1;
	gridcell_r = gridcell_r % Nphi;

	double dens = ramp_factor * weight_factor * Sigma_stream;
	if (dens < parameters::sigma0 * parameters::sigma_floor) {
	    dens = parameters::sigma0 * parameters::sigma_floor;
	}

	if (parameters::Adiabatic) {
	    const double T_stream =
		parameters::mof_temperature; // Stream Temperature in Kelvin
	    const double e_stream =
		T_stream * units::temperature.get_inverse_cgs_factor() * dens /
		parameters::MU * constants::R /
		(parameters::ADIABATICINDEX - 1.0); // energy density equivalent to T_stream

	    data[t_data::ENERGY](Nrad - 2, gridcell) = e_stream;
	}

	data[t_data::SIGMA](Nrad - 2, gridcell) += dens;

	data[t_data::V_RADIAL](Nrad - 2, gridcell) = vr_stream;
	data[t_data::V_RADIAL](Nrad - 1, gridcell) = vr_stream;
	data[t_data::V_AZIMUTHAL](Nrad - 2, gridcell) = vphi_stream;
	data[t_data::V_AZIMUTHAL](Nrad - 2, gridcell_r) = vphi_stream;

#ifndef NDEBUG
	logging::print(
	    LOG_VERBOSE
	    "dens %lE, WF %lE , angle %lf, nearest_grid_cell %i, mass_stream %lE, gridcell %i, Nphi %i , noc %i , i %i \n",
	    data[t_data::SIGMA](Nrad, gridcell), 1.0, angle, gridcell,
	    Sigma_stream, gridcell, Nrad, Nphi * 2 + 1, 0);
#endif
    }

    // check if mass overflow is constant
    if (!(check > 0.99) || !(check < 1.01)) {
	logging::print_master(LOG_ERROR "weight_factor %lf should be 0.997 \n",
			      check);
    }
}

/**
 *	Mass overflow function copied from RH2D
 *	added gaussian profile to stream
 */
void mass_overflow_willy(t_data &data, t_polargrid *densitystar, bool transport)
{
    if (CPU_Rank != CPU_Highest) {
	return;
    }

    // get location of binary star
    if (parameters::mof_planet + 1 >
	data.get_planetary_system().get_number_of_planets()) {
	logging::print_master(
	    LOG_ERROR
	    "Wrong Planet/Star for Mass Overflow specified! Old parameter file?\n");
	die("Wrong Planet/Star for Mass Overflow specified! Old parameter file?");
    }

	const double mdot_avg = data.get_massflow_tracker().get_mdot();

    const t_planet &planet =
	data.get_planetary_system().get_planet(parameters::mof_planet);
    const double xplanet = planet.get_x();
	const double yplanet = planet.get_y();

    const double omega_planet = planet.get_omega();
    // get grid cell where binary star is nearest
    // atan2(y,x) is from -PI to PI
    double angle = std::atan2(yplanet, xplanet) * 0.5 * M_1_PI;
    if (angle < 0) {
	angle += 1.0;
    }

    const unsigned int Nrad = data[t_data::SIGMA].get_max_radial();
    const unsigned int Nphi = data[t_data::SIGMA].get_size_azimuthal();
    const double r_cell = Rmed[Nrad];

    const double vr_fraction = 0.002;
    const double vr_stream =
	-omega_planet * r_cell *
	vr_fraction; // radial inflow velocity is a fraction of v_Kepler, (e.g.
		     // vr_fraction = 0.002)
    const double vphi_stream =
	(omega_planet - refframe::OmegaFrame) *
	r_cell; // set angular velocity to zero (in co-rotating frame)

	double mdot_transfer;
	if (parameters::variableTransfer){
        // variable Mdot following Hameury, Lasota & Warner 1999  eq. 4
        mdot_transfer = std::max(parameters::mof_value, parameters::mof_gamma * mdot_avg);
	}else{
		mdot_transfer = parameters::mof_value;
	}

	const double stream_mdot_code = mdot_transfer;

    
    const double Sigma_stream =
	fabs(stream_mdot_code / (dphi * Rinf[Nrad] * vr_stream));

    const int nearest_grid_cell = ((int)((double)Nphi * angle + 0.5)) %
				  Nphi; // nearest gridcell to binary star

    // Calculate sigma from temperature according to
    // Meyer & Meyer-Hofmeister (1983a) equation 17
    const double Porb =
	2.0 * M_PI / omega_planet * units::time.get_cgs_factor() / 3600.0; // TODO: check unit!
    // cross section Q
    const double Q = 2.4e13 * parameters::mof_temperature * Porb * Porb;
    // stream radius W
    const double W = std::sqrt(Q * M_1_PI);
    // circumference circ
    const double circ = 2.0 * M_PI * r_cell * units::length.get_cgs_factor(); // TODO: check unit!

    double sigma = 2.0 * W / circ;
    int number_of_cells =
	(double)Nphi * 3.0 *
	sigma; // walk through 3 sigma, so we get 99.7% accuracy
    double sigmabar = Nphi * sigma;

    const double t_ramp =
	parameters::mof_rampingtime * planet.get_orbital_period();

    double ramp_factor;
    if (sim::PhysicalTime < t_ramp) {
	ramp_factor = std::pow(std::sin(sim::PhysicalTime * M_PI_2 / t_ramp), 6);
    } else {
	ramp_factor = 1.0;
    }

    double check = 0.0;
    for (int i = -number_of_cells; i <= number_of_cells; i++) {

	// adapt gauss profile
	double weight_factor;
	int gridcell;
	if (number_of_cells == 0) {
	    gridcell = nearest_grid_cell;
	    weight_factor = 1.0;
	} else {
	    gridcell = (nearest_grid_cell + i + Nphi) % Nphi;
	    weight_factor = 1.0 / (sigmabar * std::sqrt(2.0 * M_PI)) *
			    std::exp(-1.0 / 2.0 * std::pow(i / sigmabar, 2));
	}
	check += weight_factor;

	int gridcell_r = gridcell + 1;
	gridcell_r = gridcell_r % Nphi;

	double dens = ramp_factor * weight_factor * Sigma_stream;
	if (dens < parameters::sigma0 * parameters::sigma_floor) {
	    dens = parameters::sigma0 * parameters::sigma_floor;
	}

	if (parameters::Adiabatic) {
	    const double T_stream =
		parameters::mof_temperature; // Stream Temperature in Kelvin
	    const double e_stream =
		T_stream * units::temperature.get_inverse_cgs_factor() * dens /
		parameters::MU * constants::R /
		(parameters::ADIABATICINDEX - 1.0); // energy density equivalent to T_stream

	    data[t_data::ENERGY](Nrad, gridcell) = e_stream;
	}

	if (transport && (densitystar != nullptr)) {
	    (*densitystar)(Nrad, gridcell) = dens;
	} else {
	    data[t_data::SIGMA](Nrad, gridcell) = dens;
	}

	data[t_data::V_RADIAL](Nrad, gridcell) = vr_stream;
	data[t_data::V_RADIAL](Nrad + 1, gridcell) = vr_stream;
	data[t_data::V_AZIMUTHAL](Nrad, gridcell) = vphi_stream;
	data[t_data::V_AZIMUTHAL](Nrad, gridcell_r) = vphi_stream;

#ifndef NDEBUG
	logging::print(
	    LOG_VERBOSE
	    "dens %lE, WF %lE , angle %lf, nearest_grid_cell %i, mass_stream %lE, gridcell %i, Nphi %i , noc %i , i %i \n",
	    data[t_data::SIGMA](Nrad, gridcell), 1.0, angle, gridcell,
	    Sigma_stream, gridcell, Nrad, Nphi * 2 + 1, 0);
#endif
    }

    // check if mass overflow is constant
    if (!(check > 0.99) || !(check < 1.01)) {
	logging::print_master(LOG_ERROR "weight_factor %lf should be 0.997 \n",
			      check);
    }
}

void apply_boundary_condition_temperature(t_data &data)
{
    if (CPU_Rank == 0) {
	#pragma omp parallel for
	for (unsigned int n_radial = 0; n_radial <= GHOSTCELLS_B; ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal();
		 ++n_azimuthal) {
		data[t_data::ENERGY](n_radial, n_azimuthal) =
		    data[t_data::TEMPERATURE](n_radial, n_azimuthal) *
		    data[t_data::SIGMA](n_radial, n_azimuthal) /
		    (parameters::ADIABATICINDEX - 1.0) / parameters::MU * constants::R;
	    }
	}

	switch (parameters::boundary_inner) {
	case parameters::boundary_condition_open:
	    open_boundary_inner(data);
	    break;
	case parameters::boundary_condition_reflecting:
	    reflecting_boundary_inner(data);
	    break;
	default:
	    logging::print_master(
		LOG_ERROR
		"Boundary condition at inner boundary not supported for temperature!\n");
	    break;
	}

	#pragma omp parallel for
	for (unsigned int n_radial = 0; n_radial <= GHOSTCELLS_B; ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal();
		 ++n_azimuthal) {
		data[t_data::TEMPERATURE](n_radial, n_azimuthal) =
		    data[t_data::ENERGY](n_radial, n_azimuthal) /
		    data[t_data::SIGMA](n_radial, n_azimuthal) *
		    (parameters::ADIABATICINDEX - 1.0) * parameters::MU / constants::R;
	    }
	}
    }

    if (CPU_Rank == CPU_Highest) {
	#pragma omp parallel for
	for (unsigned int n_radial =
		 data[t_data::ENERGY].get_max_radial() - GHOSTCELLS_B;
	     n_radial <= data[t_data::ENERGY].get_max_radial(); ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal();
		 ++n_azimuthal) {
		data[t_data::ENERGY](n_radial, n_azimuthal) =
		    data[t_data::TEMPERATURE](n_radial, n_azimuthal) *
		    data[t_data::SIGMA](n_radial, n_azimuthal) /
		    (parameters::ADIABATICINDEX - 1.0) / parameters::MU * constants::R;
	    }
	}

	// outer boundary
	switch (parameters::boundary_outer) {
	case parameters::boundary_condition_open:
	    open_boundary_outer(data);
	    break;
	case parameters::boundary_condition_reflecting:
	    reflecting_boundary_outer(data);
	    break;
	default:
	    logging::print_master(
		LOG_ERROR
		"Boundary condition at outer boundary not supported for temperature!\n");
	    break;
	}

	#pragma omp parallel for
	for (unsigned int n_radial =
		 data[t_data::ENERGY].get_max_radial() - GHOSTCELLS_B;
	     n_radial <= data[t_data::ENERGY].get_max_radial(); ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal();
		 ++n_azimuthal) {
		data[t_data::TEMPERATURE](n_radial, n_azimuthal) =
		    data[t_data::ENERGY](n_radial, n_azimuthal) /
		    data[t_data::SIGMA](n_radial, n_azimuthal) *
		    (parameters::ADIABATICINDEX - 1.0) * parameters::MU / constants::R;
	    }
	}
    }
}


} // namespace boundary_conditions
