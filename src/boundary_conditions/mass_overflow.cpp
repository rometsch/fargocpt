/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"
#include "../global.h"
#include "../parameters.h"
#include "../logging.h"
#include "../frame_of_reference.h"
#include "../simulation.h"
#include "../constants.h"
#include "../LowTasks.h"


namespace boundary_conditions
{

/**
 *	Mass overflow function copied from RH2D following Kley et al. 2008
 *	added gaussian profile to stream
 */
void rochelobe_overflow_boundary(t_data &data, t_polargrid *densitystar, bool transport)
{
    if (CPU_Rank != CPU_Highest) {
	return;
    }

    // get location of binary star
    if (rof_planet + 1 >
	data.get_planetary_system().get_number_of_planets()) {
	logging::print_master(
	    LOG_ERROR
	    "Wrong Planet/Star for Mass Overflow specified! Old parameter file?\n");
	die("Wrong Planet/Star for Mass Overflow specified! Old parameter file?");
    }

	const double mdot_avg = data.get_massflow_tracker().get_mdot();

    const t_planet &planet =
	data.get_planetary_system().get_planet(rof_planet);
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
	if (rof_variableTransfer){
        // variable Mdot following Hameury, Lasota & Warner 1999  eq. 4
        mdot_transfer = std::max(rof_mdot, rof_gamma * mdot_avg);
	}else{
		mdot_transfer = rof_mdot;
	}

	const double stream_mdot_code = mdot_transfer;

    
    const double Sigma_stream =
	fabs(stream_mdot_code / (dphi * Rinf[Nrad] * vr_stream));

    const int nearest_grid_cell = ((int)((double)Nphi * angle + 0.5)) %
				  Nphi; // nearest gridcell to binary star

    // Calculate sigma from temperature according to
    // Meyer & Meyer-Hofmeister (1983a) equation 17
    const double Porb =
	2.0 * M_PI / omega_planet * units::time.get_code_to_cgs_factor() / 3600.0; // TODO: check unit!
    // cross section Q
    const double Q = 2.4e13 * rof_temperature * Porb * Porb;
    // stream radius W
    const double W = std::sqrt(Q * M_1_PI);
    // circumference circ
    const double circ = 2.0 * M_PI * r_cell * units::length.get_code_to_cgs_factor(); // TODO: check unit!

    double sigma = 2.0 * W / circ;
    int number_of_cells =
	(double)Nphi * 3.0 *
	sigma; // walk through 3 sigma, so we get 99.7% accuracy
    double sigmabar = Nphi * sigma;

    const double t_ramp =
	rof_rampingtime * planet.get_orbital_period();

    double ramp_factor;
    if (sim::time < t_ramp) {
	ramp_factor = std::pow(std::sin(sim::time * M_PI_2 / t_ramp), 6);
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
		rof_temperature; // Stream Temperature in Kelvin
	    const double e_stream =
		T_stream * units::temperature.get_cgs_to_code_factor() * dens /
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

} // namespace boundary_conditions
