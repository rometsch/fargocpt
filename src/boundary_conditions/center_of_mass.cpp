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
#include "../Theo.h"
#include "../viscosity/viscous_radial_speed.h"
#include "../frame_of_reference.h"
#include "../find_cell_id.h"
#include "../constants.h"



namespace boundary_conditions
{

void init_center_of_mass(t_data& data) {
	viscous_speed::init_vr_table_boundary(data);
}


/**
 * @brief initial_center_of_mass_boundary: sets the outer boundary
 *  to the initial profile in the center of mass and then shifts it to the
 * primary center. Lucas thinks this is important when simulating a circumbinary
 * disk when the coordination system is centered on the primary.
 * @param data
 */
void diskmodel_center_of_mass_boundary_outer(t_data &data)
{

    if (CPU_Rank != CPU_Highest) {
		return;
	}

    const unsigned int np = data.get_planetary_system().get_number_of_planets();
    const Pair com_pos = data.get_planetary_system().get_center_of_mass(np);
    const Pair com_vel = data.get_planetary_system().get_center_of_mass_velocity(np);
    const double com_mass = data.get_planetary_system().get_mass(np);

    auto &sigma = data[t_data::SIGMA];
    auto &energy = data[t_data::ENERGY];
    auto &vrad = data[t_data::V_RADIAL];
    auto &vaz = data[t_data::V_AZIMUTHAL];


	const unsigned int Naz = data[t_data::SIGMA].get_max_azimuthal();
    const unsigned int nr = data[t_data::SIGMA].get_max_radial();
	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {

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
			 -parameters::sigma_slope); // we assume the floor is not reached.
	    sigma(nr, naz) = cell_sigma;

	    /// Initial profile temperature
	const double cell_energy = initial_energy(r_com, com_mass);

	const double temperature_floor =
	    parameters::minimum_temperature *
	    units::temperature.get_cgs_to_code_factor();

	const double energy_floor = temperature_floor * cell_sigma /
				    parameters::MU * constants::R /
				    (parameters::ADIABATICINDEX - 1.0);

	energy(nr, naz) = std::max(cell_energy, energy_floor);

	    /// dP / dr = 0
		/// energy(nr, naz) = energy(nr - 1, naz);
	} /// END DENSITY and ENERGY
    }
}


void diskmodel_center_of_mass_boundary_inner(t_data &data)
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
			 -parameters::sigma_slope); // we assume the floor is not reached.
		sigma(nr, naz) = cell_sigma;

		/// Initial profile temperature
	const double cell_energy = initial_energy(r_com, com_mass);

	const double temperature_floor =
		parameters::minimum_temperature *
		units::temperature.get_cgs_to_code_factor();

	const double energy_floor = temperature_floor * cell_sigma /
					parameters::MU * constants::R /
					(parameters::ADIABATICINDEX - 1.0);

	energy(nr, naz) = std::max(cell_energy, energy_floor);

		/// dP / dr = 0
		/// energy(nr, naz) = energy(nr - 1, naz);
	} /// END DENSITY and ENERGY
	}
}


void damping_diskmodel_center_of_mass_outer(t_data &data, double dt)
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
    if ((damping_outer_limit < 1.0) &&
	(Rinf[vrad_arr.get_max_radial()] >
	 RMAX * damping_outer_limit)) {

	const unsigned int clamped_vrad_id = clamp_r_id_to_radii_grid(
		get_rinf_id(RMAX * damping_outer_limit), vrad_arr.is_vector())+1;

    const double tau = damping_time_factor * 2.0 * M_PI /
			 calculate_omega_kepler(damping_time_radius_outer);

	#pragma omp parallel for
	for (unsigned int n_radial = clamped_vrad_id;
		 n_radial < MaxMo_no_ghost_vr; ++n_radial) {
	    double factor = std::pow(
		(Rinf[n_radial] - RMAX * damping_outer_limit) /
		    (RMAX - RMAX * damping_outer_limit),
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
		get_rmed_id(RMAX * damping_outer_limit),
		vphi_arr.is_vector()) + 1;

	#pragma omp parallel for
	for (unsigned int n_radial = clamped_vphi_id;
		 n_radial < Max_no_ghost; ++n_radial) {
	    double factor = std::pow(
		(Rmed[n_radial] - RMAX * damping_outer_limit) /
		    (RMAX - RMAX * damping_outer_limit),
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
		(Rmed[nr] - RMAX * damping_outer_limit) /
			(RMAX - RMAX * damping_outer_limit),
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

void damping_diskmodel_center_of_mass_inner(t_data &data, double dt)
{

	t_polargrid &vrad_arr = data[t_data::V_RADIAL];
	t_polargrid &vphi_arr = data[t_data::V_AZIMUTHAL];

	const unsigned int np = parameters::n_bodies_for_hydroframe_center;
	const Pair com_pos = data.get_planetary_system().get_center_of_mass(np);
	const Pair com_vel =
	data.get_planetary_system().get_center_of_mass_velocity(np);
	const double com_mass = data.get_planetary_system().get_mass(np);

	// is this CPU in the inner damping domain?
	if ((damping_inner_limit > 1.0) &&
	(Rmed[0] < RMIN * damping_inner_limit)) {

	const unsigned int clamped_vrad_id = clamp_r_id_to_radii_grid(
		get_rinf_id(RMIN * damping_inner_limit),
		vrad_arr.is_vector());

	double tau = damping_time_factor * 2.0 * M_PI /
			 calculate_omega_kepler(RMIN);

	#pragma omp parallel for
	for (unsigned int n_radial = One_no_ghost_vr;
		 n_radial <= clamped_vrad_id; ++n_radial) {
		double factor = std::pow(
		(Rinf[n_radial] - RMIN * damping_inner_limit) /
			(RMIN - RMIN * damping_inner_limit),
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
		get_rmed_id(RMIN * damping_inner_limit),
		vphi_arr.is_vector());

	#pragma omp parallel for
	for (unsigned int n_radial = Zero_no_ghost;
		 n_radial <= clamped_vphi_id; ++n_radial) {
		double factor = std::pow(
		(Rmed[n_radial] - RMAX * damping_outer_limit) /
			(RMAX - RMAX * damping_outer_limit),
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
		(Rmed[nr] - RMAX * damping_outer_limit) /
			(RMAX - RMAX * damping_outer_limit),
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


} // namespace boundary_conditions
