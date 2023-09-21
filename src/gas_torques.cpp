#include "gas_torques.h"
#include "constants.h"
#include "data.h"
#include "global.h"
#include "parameters.h"
#include "util.h"
#include <cmath>

namespace gas_torques
{

void calculate_advection_torque(t_data &data, t_polargrid &t_adv, const double dt)
{

    t_polargrid &vphi = data[t_data::V_AZIMUTHAL];
    t_polargrid &vr = data[t_data::V_RADIAL];
    t_polargrid &dens = data[t_data::SIGMA];

    for (unsigned int n_r = 0; n_r < t_adv.get_size_radial(); ++n_r) {

	const double r = Rmed[n_r];
	const double inv_dr = InvDiffRsup[n_r];

	for (unsigned int n_az = 0; n_az < t_adv.get_size_azimuthal(); ++n_az) {
	    const double sigma_cell = dens(n_r, n_az);

	    const double vr_top = vr(n_r + 1, n_az);
	    const double vr_bot = vr(n_r, n_az);

	    // average to cell center
	    double vr_cell = (Rmed[n_r] - Rinf[n_r]) * vr_top +
			     (Rsup[n_r] - Rmed[n_r]) * vr_bot;
	    vr_cell *= inv_dr;

	    const double vphi_cell_left = vphi(n_r, n_az);
	    const double vphi_cell_right =
		vphi(n_r, get_next_azimuthal_id(n_az));
	    const double vphi_cell = 0.5 * (vphi_cell_left + vphi_cell_right);

	    t_adv(n_r, n_az) +=
		-std::pow(r, 2.0) * sigma_cell * vr_cell * vphi_cell * dt;
	}
    }
}

void calculate_viscous_torque(t_data &data, t_polargrid &t_visc, const double dt)
{
    t_polargrid &visc = data[t_data::VISCOSITY];
    t_polargrid &vphi = data[t_data::V_AZIMUTHAL];
    t_polargrid &vr = data[t_data::V_RADIAL];
    t_polargrid &dens = data[t_data::SIGMA];

    for (unsigned int n_r = 1; n_r < t_visc.get_max_radial(); ++n_r) {

	const double r = Rmed[n_r];
	const double inv_dr = InvDiffRsup[n_r];
	const double inv_dr_med_top = InvDiffRmed[n_r + 1];
	const double inv_dr_med_bot = InvDiffRmed[n_r];

	for (unsigned int n_az = 0; n_az < t_visc.get_size_azimuthal();
	     ++n_az) {
	    const double sigma_cell = dens(n_r, n_az);
	    const double viscosity_cell = visc(n_r, n_az);

	    const double vr_cellm1_top =
		vr(n_r + 1, get_prev_azimuthal_id(n_az));
	    const double vr_cellp1_top =
		vr(n_r + 1, get_next_azimuthal_id(n_az));

	    const double dvr_dphi_top =
		(vr_cellp1_top - vr_cellm1_top) * 0.5 * invdphi;

	    const double vr_cellm1 = vr(n_r, get_prev_azimuthal_id(n_az));
	    const double vr_cellp1 = vr(n_r, get_next_azimuthal_id(n_az));

	    const double dvr_dphi_bot = (vr_cellp1 - vr_cellm1) * 0.5 * invdphi;

	    // average to cell center
	    double dvr_dphi = (Rmed[n_r] - Rinf[n_r]) * dvr_dphi_top +
			      (Rsup[n_r] - Rmed[n_r]) * dvr_dphi_bot;
	    dvr_dphi *= inv_dr;

	    const double vphi_cellm1_top = vphi(n_r + 1, n_az);
	    const double vphi_cellp1_top =
		vphi(n_r + 1, get_next_azimuthal_id(n_az));
	    const double phi_dot_top =
		0.5 * (vphi_cellp1_top + vphi_cellm1_top) / Rmed[n_r + 1];

	    const double vphi_cellm1 = vphi(n_r, n_az);
	    const double vphi_cellp1 = vphi(n_r, get_next_azimuthal_id(n_az));
	    const double phi_dot =
		0.5 * (vphi_cellp1 + vphi_cellm1) / Rmed[n_r];

	    const double vphi_cellm1_bot = vphi(n_r - 1, n_az);
	    const double vphi_cellp1_bot =
		vphi(n_r - 1, get_next_azimuthal_id(n_az));
	    const double phi_dot_bot =
		0.5 * (vphi_cellp1_bot + vphi_cellm1_bot) / Rmed[n_r - 1];

	    // derivatives at the center of top and bottom interface
	    const double dphi_dot_dr_top =
		(phi_dot_top - phi_dot) * inv_dr_med_top;
	    const double dphi_dot_dr_bot =
		(phi_dot - phi_dot_bot) * inv_dr_med_bot;

	    // average to cell center
	    double dphi_dot_dr = (Rmed[n_r] - Rinf[n_r]) * dphi_dot_dr_top +
				 (Rsup[n_r] - Rmed[n_r]) * dphi_dot_dr_bot;
	    dphi_dot_dr *= inv_dr;

	    t_visc(n_r, n_az) +=
		-std::pow(r, 3.0) * viscosity_cell * sigma_cell *
		(dphi_dot_dr + 1.0 / (std::pow(r, 2.0)) * dvr_dphi) * dt;
	}
    }
}

/**
 * @brief calculate_gravitational_torque, see Eq. (32) in Miranda 2017.
 * @param data
 * @param dt
 */
void calculate_gravitational_torque(t_data &data,  t_polargrid &t_grav, const double dt)
{

    t_polargrid &pot_grav = data[t_data::POTENTIAL];

    for (unsigned int n_r = 0; n_r < t_grav.get_size_radial(); ++n_r) {

	const double r = Rmed[n_r];

	for (unsigned int n_az = 0; n_az < t_grav.get_size_azimuthal();
	     ++n_az) {
	    const double sigma_cell = data[t_data::SIGMA](n_r, n_az);

	    // dPhi/dphi
	    // const double gradphi_minus = (pot_grav(n_r, n_az) - pot_grav(n_r,
	    // n_az == 0 ? pot_grav.get_max_azimuthal() : n_az - 1))*invdphi;
	    // const double gradphi = (pot_grav(n_r, n_az ==
	    // pot_grav.get_max_azimuthal() ? 0 : n_az+1) - pot_grav(n_r,
	    // n_az))*invdphi;
	    double gradphi;
	    if (parameters::body_force_from_potential) {
		gradphi = (pot_grav(n_r, get_next_azimuthal_id(n_az)) -
			   pot_grav(n_r, get_prev_azimuthal_id(n_az))) *
			  invdphi * 0.5;
	    } else {
		gradphi = -data[t_data::ACCEL_AZIMUTHAL](n_r, n_az) * r;
	    }
	    //t_grav(n_r, n_az) += -r * sigma_cell * gradphi * dr * dphi * dt;
	    // Surf = r * dr * dphi
	    t_grav(n_r, n_az) += -sigma_cell * gradphi * Surf[n_r] * dt;
	}
    }
}

} // namespace gas_torques
