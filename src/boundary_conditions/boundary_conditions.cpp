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
	jibin_boundary_inner(data);
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
	jibin_boundary_outer(data);
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

/**
	Boundary conditions for calculation of the boundary layer (BL) starting
   here:
	TODO: Stellar radiative flux into disk via implicit routine!
*/

/**
	Inner boundary: zero gradient & fixed velocities
*/

void boundary_layer_inner_boundary(t_data &data)
{

    if (CPU_Rank != 0)
	return;

	#pragma omp parallel for
    for (unsigned int n_azimuthal = 0;
	 n_azimuthal <= data[t_data::SIGMA].get_max_azimuthal();
	 ++n_azimuthal) {
	// zero gradient
	data[t_data::SIGMA](0, n_azimuthal) =
	    data[t_data::SIGMA](1, n_azimuthal);
	data[t_data::ENERGY](0, n_azimuthal) =
	    data[t_data::ENERGY](1, n_azimuthal);

	// set vrad to fraction of Keplerian velocity
	data[t_data::V_RADIAL](1, n_azimuthal) =
	    -1. * parameters::vrad_fraction_of_kepler * std::sqrt(1. / Ra[1]);
	data[t_data::V_RADIAL](0, n_azimuthal) =
	    data[t_data::V_RADIAL](1, n_azimuthal);

	// set vphi to stellar rotation rate 
	data[t_data::V_AZIMUTHAL](0, n_azimuthal) =
	    parameters::stellar_rotation_rate * std::sqrt(1. / Rb[0]);
    }
}

/**
	Outer boundary: floating boundary conditions & pressure correction for
   Omega
*/

void boundary_layer_outer_boundary(t_data &data)
{

    if (CPU_Rank != CPU_Highest)
	return;

	#pragma omp parallel for
    for (unsigned int n_azimuthal = 0;
	 n_azimuthal <= data[t_data::SIGMA].get_max_azimuthal();
	 ++n_azimuthal) {
	// floating BCs
	data[t_data::SIGMA](data[t_data::SIGMA].get_max_radial(), n_azimuthal) =
	    data[t_data::SIGMA](data[t_data::SIGMA].get_max_radial() - 1,
				n_azimuthal) *
	    std::sqrt(Ra[data[t_data::SIGMA].get_max_radial() - 1] /
		      Ra[data[t_data::SIGMA].get_max_radial()]);
	data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial(),
			     n_azimuthal) =
	    data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial() - 1,
				 n_azimuthal) *
	    std::pow(Ra[data[t_data::ENERGY].get_max_radial() - 1] /
			 Ra[data[t_data::ENERGY].get_max_radial()],
		     1.25);

	data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial() - 1,
			       n_azimuthal) =
	    -1. *
	    fabs(data[t_data::V_RADIAL](
		data[t_data::V_RADIAL].get_max_radial() - 2, n_azimuthal)) *
	    std::sqrt(Ra[data[t_data::V_RADIAL].get_max_radial() - 2] /
		      Ra[data[t_data::V_RADIAL].get_max_radial() - 1]);

	data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial(),
			       n_azimuthal) =
	    -1. *
	    fabs(data[t_data::V_RADIAL](
		data[t_data::V_RADIAL].get_max_radial() - 2, n_azimuthal)) *
	    std::sqrt(Ra[data[t_data::V_RADIAL].get_max_radial() - 2] /
		      Ra[data[t_data::V_RADIAL].get_max_radial()]);

	// Omega at outer boundary equals calculate_omega_kepler (plus leading
	// order pressure correction)
	data[t_data::V_AZIMUTHAL](data[t_data::V_AZIMUTHAL].get_max_radial(),
				  n_azimuthal) =
	    1. / std::sqrt(Rb[data[t_data::SIGMA].get_max_radial()]);
	// TODO: Include pressure correction, like in uphi[*jN]
	// = 1./sqrt(Rb[*jN]) +
	// 0.5/Sigma[*jN]*sqrt(pow3(Rb[*jN])*pow2(Rb[*jN]))*.5/Rb[*jN]*(P[*jN+1]-P[*jN-1])/DeltaRa[*jN+1];
    }
}

/**
	End BL Boundary Conditions!
*/


/**
 * @brief initial_center_of_mass_boundary: sets the outer boundary
 *  to the initial profile in the center of mass and then shifts it to the
 * primary center. Lucas thinks this is important when simulating a circumbinary
 * disk when the coordination system is centered on the primary.
 * @param data
 */
void initial_center_of_mass_boundary_outer(t_data &data)
{

    if (CPU_Rank != CPU_Highest)
	return;

    const unsigned int np = data.get_planetary_system().get_number_of_planets();
    const Pair com_pos = data.get_planetary_system().get_center_of_mass(np);
    const Pair com_vel =
	data.get_planetary_system().get_center_of_mass_velocity(np);
    const double com_mass = data.get_planetary_system().get_mass(np);

    auto &sigma = data[t_data::SIGMA];
    auto &energy = data[t_data::ENERGY];
    auto &vrad = data[t_data::V_RADIAL];
    auto &vaz = data[t_data::V_AZIMUTHAL];

    const unsigned int nr = data[t_data::SIGMA].get_max_radial();
	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= data[t_data::SIGMA].get_max_azimuthal();
	 ++naz) {

	{ /// V_PHI
	    const double phi = ((double)naz - 0.5) * dphi;
	    const double rmed = Rmed[nr];

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
			vr0 = initial_viscous_radial_speed(r_com, com_mass);
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

	    const double cell_vphi =
		(cell_x * cell_vy - cell_vx * cell_y) / rmed;
		vaz(nr, naz) = cell_vphi - refframe::OmegaFrame * rmed;
	} /// END V_PHI

	{ /// V_R
	    const double phi = (double)naz * dphi;
	    const double rinf = Rinf[nr];

	    const double cell_x = rinf * std::cos(phi);
	    const double cell_y = rinf * std::sin(phi);

	    // Position in center of mass frame
	    const double x_com = cell_x - com_pos.x;
	    const double y_com = cell_y - com_pos.y;
	    const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

	    // pressure support correction
		double vphi0;
	    double vr0;
	    if (parameters::initialize_pure_keplerian) {
			vphi0 = compute_v_kepler(r_com, com_mass);
			vr0 = initial_viscous_radial_speed(r_com, com_mass);
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

	    const double cell_vr = (cell_x * cell_vx + cell_y * cell_vy) / rinf;
	    vrad(nr, naz) = cell_vr;
	} /// END V_R

	{ /// V_R GHOST CELL
	    const double phi = (double)naz * dphi;
	    const double rinf = Rsup[nr];

	    const double cell_x = rinf * std::cos(phi);
	    const double cell_y = rinf * std::sin(phi);

	    // Position in center of mass frame
	    const double x_com = cell_x - com_pos.x;
	    const double y_com = cell_y - com_pos.y;
	    const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

	    // pressure support correction
		double vphi0;
	    double vr0;
	    if (parameters::initialize_pure_keplerian) {
			vphi0 = compute_v_kepler(r_com, com_mass);
			vr0 = initial_viscous_radial_speed(r_com, com_mass);
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

	    const double cell_vr = (cell_x * cell_vx + cell_y * cell_vy) / rinf;
	    vrad(nr + 1, naz) = cell_vr;
	} /// END V_R GHOST CELL

	{ /// DENSITY and ENERGY
	    const double cell_x = (*CellCenterX)(nr, naz);
	    const double cell_y = (*CellCenterY)(nr, naz);

	    // Position in center of mass frame
	    const double x_com = cell_x - com_pos.x;
	    const double y_com = cell_y - com_pos.y;
	    const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

	    const double cell_sigma =
		parameters::sigma0 *
		std::pow(r_com,
			 -parameters::SIGMASLOPE); // we assume the floor is not reached.
	    sigma(nr, naz) = cell_sigma;

	    /// Initial profile temperature
	const double cell_energy = initial_energy(r_com, com_mass);

	const double temperature_floor =
	    parameters::minimum_temperature *
	    units::temperature.get_inverse_cgs_factor();

	const double energy_floor = temperature_floor * cell_sigma /
				    parameters::MU * constants::R /
				    (parameters::ADIABATICINDEX - 1.0);

	energy(nr, naz) = std::max(cell_energy, energy_floor);

	    /// dP / dr = 0
		/// energy(nr, naz) = energy(nr - 1, naz);
	} /// END DENSITY and ENERGY
    }
}


void initial_center_of_mass_boundary_inner(t_data &data)
{

	if (CPU_Rank != 0)
	return;

	const unsigned int np = parameters::n_bodies_for_hydroframe_center;
	const Pair com_pos = data.get_planetary_system().get_center_of_mass(np);
	const Pair com_vel =
	data.get_planetary_system().get_center_of_mass_velocity(np);
	const double com_mass = data.get_planetary_system().get_mass(np);

	auto &sigma = data[t_data::SIGMA];
	auto &energy = data[t_data::ENERGY];
	auto &vrad = data[t_data::V_RADIAL];
	auto &vaz = data[t_data::V_AZIMUTHAL];

	const unsigned int nr = 0;
	#pragma omp parallel for
	for (unsigned int naz = 0; naz <= data[t_data::SIGMA].get_max_azimuthal();
	 ++naz) {

	{ /// V_PHI
		const double phi = ((double)naz - 0.5) * dphi;
		const double rmed = Rmed[nr];

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
			vr0 = initial_viscous_radial_speed(r_com, com_mass);
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

		const double cell_vphi =
		(cell_x * cell_vy - cell_vx * cell_y) / rmed;
		vaz(nr, naz) = cell_vphi - refframe::OmegaFrame * rmed;
	} /// END V_PHI

	{ /// V_R GHOST CELL
		const double phi = (double)naz * dphi;
		const double rinf = Rinf[nr];

		const double cell_x = rinf * std::cos(phi);
		const double cell_y = rinf * std::sin(phi);

		// Position in center of mass frame
		const double x_com = cell_x - com_pos.x;
		const double y_com = cell_y - com_pos.y;
		const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

		// pressure support correction
		double vphi0;
		double vr0;
		if (parameters::initialize_pure_keplerian) {
			vphi0 = compute_v_kepler(r_com, com_mass);
			vr0 = initial_viscous_radial_speed(r_com, com_mass);
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

		const double cell_vr = (cell_x * cell_vx + cell_y * cell_vy) / rinf;
		vrad(nr, naz) = cell_vr;
	} /// END V_R

	{ /// V_R
		const double phi = (double)naz * dphi;
		const double rinf = Rsup[nr];

		const double cell_x = rinf * std::cos(phi);
		const double cell_y = rinf * std::sin(phi);

		// Position in center of mass frame
		const double x_com = cell_x - com_pos.x;
		const double y_com = cell_y - com_pos.y;
		const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

		// pressure support correction
		double vphi0;
		double vr0;
		if (parameters::initialize_pure_keplerian) {
			vphi0 = compute_v_kepler(r_com, com_mass);
			vr0 = initial_viscous_radial_speed(r_com, com_mass);
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

		const double cell_vr = (cell_x * cell_vx + cell_y * cell_vy) / rinf;
		vrad(nr + 1, naz) = cell_vr;
	} /// END V_R GHOST CELL

	{ /// DENSITY and ENERGY
		const double cell_x = (*CellCenterX)(nr, naz);
		const double cell_y = (*CellCenterY)(nr, naz);

		// Position in center of mass frame
		const double x_com = cell_x - com_pos.x;
		const double y_com = cell_y - com_pos.y;
		const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

		const double cell_sigma =
		parameters::sigma0 *
		std::pow(r_com,
			 -parameters::SIGMASLOPE); // we assume the floor is not reached.
		sigma(nr, naz) = cell_sigma;

		/// Initial profile temperature
	const double cell_energy = initial_energy(r_com, com_mass);

	const double temperature_floor =
		parameters::minimum_temperature *
		units::temperature.get_inverse_cgs_factor();

	const double energy_floor = temperature_floor * cell_sigma /
					parameters::MU * constants::R /
					(parameters::ADIABATICINDEX - 1.0);

	energy(nr, naz) = std::max(cell_energy, energy_floor);

		/// dP / dr = 0
		/// energy(nr, naz) = energy(nr - 1, naz);
	} /// END DENSITY and ENERGY
	}
}

/**
	for viscous spreading ring comparison simulations for Jibin
*/
void jibin_boundary_inner(t_data &data)
{
    if (CPU_Rank != 0)
	return;

    const double R = Rmed[0];

	const double vaz = initial_locally_isothermal_smoothed_v_az(R, 1.0) - R * refframe::OmegaFrame;


	#pragma omp parallel for
    for (unsigned int n_azimuthal = 0;
	 n_azimuthal <= data[t_data::SIGMA].get_max_azimuthal();
	 ++n_azimuthal) {

	// copy first ring into ghost ring
	data[t_data::SIGMA](0, n_azimuthal) =
	    data[t_data::SIGMA](1, n_azimuthal);

	data[t_data::V_AZIMUTHAL](0, n_azimuthal) = vaz;

	if (data[t_data::V_RADIAL](2, n_azimuthal) <= 0.0) { // outflow
	    data[t_data::V_RADIAL](1, n_azimuthal) =
		data[t_data::V_RADIAL](2, n_azimuthal);
	    data[t_data::V_RADIAL](0, n_azimuthal) =
		data[t_data::V_RADIAL](2, n_azimuthal);
	} else { // reflective
	    data[t_data::V_RADIAL](1, n_azimuthal) =
		-data[t_data::V_RADIAL](2, n_azimuthal);
	    data[t_data::V_RADIAL](0, n_azimuthal) =
		-data[t_data::V_RADIAL](2, n_azimuthal);
	}
    }
}

/**
	for viscous spreading ring comparison simulations for Jibin
*/
void jibin_boundary_outer(t_data &data)
{
    if (CPU_Rank == CPU_Highest) {

	const double R = Rmed[data[t_data::V_AZIMUTHAL].get_max_radial()];

	const double vaz = initial_locally_isothermal_smoothed_v_az(R, 1.0) - R * refframe::OmegaFrame;


	#pragma omp parallel for
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::SIGMA].get_max_azimuthal();
	     ++n_azimuthal) {
	    data[t_data::V_AZIMUTHAL](
		data[t_data::V_AZIMUTHAL].get_max_radial(), n_azimuthal) = vaz;
	}
    }
}

} // namespace boundary_conditions
