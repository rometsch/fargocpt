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

} // namespace boundary_conditions
