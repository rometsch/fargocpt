#include "artificial_viscosity.h"

#include "../parameters.h"
#include "../SideEuler.h"
#include "../global.h"
#include "../SourceEuler.h"
#include "../quantities.h"
#include <cassert>

namespace art_visc{

void update_with_artificial_viscosity(t_data &data, const double dt){
	if (parameters::artificial_viscosity ==
	 parameters::artificial_viscosity_TW) {
	art_visc::update_with_artificial_viscosity_TW(data, dt);
	} else { // SN or TW (just for Vazimuthal ghost cells)
	art_visc::update_with_artificial_viscosity_SN(data, dt);
	}

	if (parameters::Adiabatic && parameters::artificial_viscosity_dissipation) {
		SetTemperatureFloorCeilValues(data, __FILE__, __LINE__);
	}

	if(ECC_GROWTH_MONITOR){
		quantities::calculate_disk_delta_ecc_peri(data, delta_ecc_art_visc, delta_peri_art_visc);
	}
}


/**
 * @brief update_with_artificial_viscosity_TW: artificial viscosity as described in Tscharnuter & Winkler 1973
 * off diagonal tensor elements are set to zero to not inflict artificial angular momentum transfer (compare D'angelo et al. 2003)
 * @param data
 * @param dt
 */
void update_with_artificial_viscosity_TW(t_data &data, const double dt)
{

	/// Do not Apply sub keplerian boundary for boundary conditions that set
	/// Vphi themselves
	const bool add_kep_inner =
	(parameters::boundary_inner !=
	parameters::boundary_condition_center_of_mass_initial) &&
	(parameters::boundary_inner !=
	 parameters::boundary_condition_evanescent) &&
	(parameters::boundary_inner !=
	 parameters::boundary_condition_boundary_layer) &&
	(parameters::boundary_inner !=
	 parameters::boundary_condition_precribed_time_variable) &&
	(!parameters::domegadr_zero);

	if (add_kep_inner) {
	ApplyKeplerianBoundaryInner(data[t_data::V_AZIMUTHAL]);
	}

	if ((parameters::boundary_outer !=
	 parameters::boundary_condition_center_of_mass_initial) &&
	(parameters::boundary_outer !=
	 parameters::boundary_condition_zero_gradient) &&
	(parameters::boundary_outer !=
	 parameters::boundary_condition_evanescent) &&
	(parameters::boundary_outer !=
	 parameters::boundary_condition_boundary_layer) &&
	(parameters::boundary_outer !=
	 parameters::boundary_condition_precribed_time_variable) &&
	(!parameters::massoverflow) && (!parameters::domegadr_zero)) {
	ApplySubKeplerianBoundaryOuter(data[t_data::V_AZIMUTHAL],
					   add_kep_inner);
	}

	t_polargrid &energy = data[t_data::ENERGY];
	t_polargrid &density = data[t_data::SIGMA];
	t_polargrid &vr = data[t_data::V_RADIAL];
	t_polargrid &vphi = data[t_data::V_AZIMUTHAL];

	t_polargrid &Q_rr = data[t_data::Q_R];
	t_polargrid &Q_pp = data[t_data::Q_PHI];

	const unsigned int Nr = density.get_size_radial();
	const unsigned int Nphi = density.get_size_azimuthal();

	// calculate div(v)
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const double naz_next = naz == vphi.get_max_azimuthal() ? 0 : naz + 1;

		// div(v) = 1/r d(r v_r)/dr + 1/r d(v_phi)/dphi
		//  	 == d(v_r)/dr + 1/r [ d(v_phi)/dphi + v_r]
		const double eps_rr = (vr(nr+1, naz) - vr(nr, naz)) * InvDiffRsup[nr];
		const double eps_pp =  InvRb[nr] * ((vphi(nr, naz_next) - vphi(nr, naz)) * invdphi + 0.5*(vr(nr + 1, naz) + vr(nr, naz)));

		const double div_V =  std::min(eps_rr + eps_pp, 0.0);

		const double Dr = Ra[nr+1] - Ra[nr];
		const double rDphi = Rmed[nr] * dphi;
		double dx_sq;
		if (Nphi <= 16) {
			// pseudo 1D simulation, we don't care about the azimuthal extent then and need this fix
			// TODO maybe: replace the arbitrary hardcoded threshold of 16
			dx_sq = std::pow(std::min(Dr, rDphi), 2);
		} else {
        	dx_sq = std::pow(std::max(Dr, rDphi), 2); // taking max cell size breaks the ShockTube test due to rDphi being huge (Nphi = 2)
		}
		const double l_sq = std::pow(parameters::artificial_viscosity_factor, 2) * dx_sq;

		const double q_rr = l_sq * density(nr, naz) * -div_V * (eps_rr - 1.0/3.0 * div_V);
		const double q_pp = l_sq * density(nr, naz) * -div_V * (eps_pp - 1.0/3.0 * div_V);

		Q_rr(nr, naz) = q_rr;
		Q_pp(nr, naz) = q_pp;

		if (parameters::Adiabatic && parameters::artificial_viscosity_dissipation) {
		if(nr > Zero_no_ghost && nr < Max_no_ghost){

			const double Qplus = - l_sq * div_V * density(nr, naz) * 1.0/3.0 * (std::pow(eps_rr, 2) + std::pow(eps_pp, 2) + std::pow((eps_rr - eps_pp), 2));
			energy(nr, naz) += Qplus * dt;
		}

	}
	}
	}

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nr-1; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {

		const int naz_prev = (naz == 0 ? Nphi-1 : naz - 1);

		const double sigma_phi_avg = 0.5 * (density(nr, naz) + density(nr, naz_prev));

		// See D'Angelo et al. 2002 Nested-grid calculations of disk-planet
		// interaction It is important to use the conservative form here and
		// not the one from Fargo / Baruteau. a_phi = 1/(r*Sigma) ( 1/r
		// d(r^2 * tau_r_phi)/dr + d(tau_phi_phi)/dphi )


		// dVp / dt = 1/rho 1/r w_q
		// w_q = dQ_pp / d phi
		/*double dVp =
		dt / (Rmed[nr] * sigma_phi_avg) *
		(Q_pp(nr, naz) - Q_pp(nr, naz_prev)) * invdphi;*/


		// Conservative volume integral formulation
		double dVp =
		2.0 * dt / ((Rsup[nr] + Rinf[nr]) * sigma_phi_avg) *
		(Q_pp(nr, naz) - Q_pp(nr, naz_prev)) * invdphi;

		vphi(nr, naz) += dVp;
	}}

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = One_no_ghost_vr; nr < MaxMo_no_ghost_vr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		// a_r = 1/(r*Sigma) ( d(r*tau_r_r)/dr + d(tau_r_phi)/dphi -
		// tau_phi_phi )
		const double sigma_r_avg = 0.5 * (density(nr, naz) + density(nr - 1, naz));

		// dVr / dt = 1/rho u_q
		// u_q = dQ_rr / dr + 1/r Q_rr - 1/r Q_pp

		// Conservative volume integral formulation
		// TODO: interpolation
		double dVr =
				parameters::radial_viscosity_factor * dt / sigma_r_avg *
				2.0 / (std::pow(Rmed[nr], 2) - std::pow(Rmed[nr-1], 2)) *
				((Q_rr(nr, naz)*Rmed[nr] - Q_rr(nr - 1, naz)*Rmed[nr-1])
				 - 0.5 * (Q_pp(nr, naz) + Q_pp(nr - 1, naz))*(Rmed[nr] - Rmed[nr-1]));

		vr(nr, naz) += dVr;
	}
	}
}

/**
	In this substep we add the articifial viscous pressure source terms.
	Shocks are spread over CVNR zones: von Neumann-Richtmyer viscosity
   constant; Beware of misprint in Stone and Norman's paper : use C2^2 instead
   of C2
*/
void update_with_artificial_viscosity_SN(t_data &data, const double dt)
{

	/// Do not Apply sub keplerian boundary for boundary conditions that set
	/// Vphi themselves
	const bool add_kep_inner =
	(parameters::boundary_inner !=
	parameters::boundary_condition_center_of_mass_initial) &&
	(parameters::boundary_inner !=
	 parameters::boundary_condition_evanescent) &&
	(parameters::boundary_inner !=
	 parameters::boundary_condition_boundary_layer) &&
	(parameters::boundary_inner !=
	 parameters::boundary_condition_precribed_time_variable) &&
	(!parameters::domegadr_zero);

	if (add_kep_inner) {
	ApplyKeplerianBoundaryInner(data[t_data::V_AZIMUTHAL]);
	}

	if ((parameters::boundary_outer !=
	 parameters::boundary_condition_center_of_mass_initial) &&
	(parameters::boundary_outer !=
	 parameters::boundary_condition_zero_gradient) &&
	(parameters::boundary_outer !=
	 parameters::boundary_condition_evanescent) &&
	(parameters::boundary_outer !=
	 parameters::boundary_condition_boundary_layer) &&
	(parameters::boundary_outer !=
	 parameters::boundary_condition_precribed_time_variable) &&
	(!parameters::massoverflow) && (!parameters::domegadr_zero)) {
	ApplySubKeplerianBoundaryOuter(data[t_data::V_AZIMUTHAL],
					   add_kep_inner);
	}

	if (parameters::artificial_viscosity ==
		parameters::artificial_viscosity_SN) {

	t_polargrid &energy = data[t_data::ENERGY];
	t_polargrid &density = data[t_data::SIGMA];
	t_polargrid &vr = data[t_data::V_RADIAL];
	t_polargrid &vphi = data[t_data::V_AZIMUTHAL];

	t_polargrid &Qr = data[t_data::Q_R];
	t_polargrid &Qphi = data[t_data::Q_PHI];

	// calculate q_r and q_phi
	const unsigned int Nr = Qr.get_size_radial();
	const unsigned int Nphi = Qr.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {
		double dv_r = vr(nr+1, naz) - vr(nr, naz);

		if (dv_r < 0.0) {
			Qr(nr, naz) =
			std::pow(parameters::artificial_viscosity_factor, 2) *
			density(nr, naz) * std::pow(dv_r, 2);
		} else {
			Qr(nr, naz) = 0.0;
		}

		const double naz_next = naz == Nphi-1 ? 0 : naz + 1;
		double dv_phi =	vphi(nr, naz_next) - vphi(nr, naz);

		if (dv_phi < 0.0) {
			Qphi(nr, naz) =
			std::pow(parameters::artificial_viscosity_factor, 2) *
			density(nr, naz) * std::pow(dv_phi, 2);
		} else {
			Qphi(nr, naz) = 0.0;
		}
		}
	}

	// If gas disk is adiabatic, we add artificial viscosity as a source
	// term for advection of thermal energy polargrid
	// perform this update before the velocities are updated
	if (parameters::Adiabatic) {
		if (parameters::artificial_viscosity_dissipation) {
		#pragma omp parallel for
		for (unsigned int nr = Zero_no_ghost; nr < Max_no_ghost; ++nr) {

			const double dxtheta = dphi * Rmed[nr];
			const double invdxtheta = 1.0 / dxtheta;

			for (unsigned int naz = 0; naz < Nphi; ++naz) {

			const double naz_next = naz == Nphi-1 ? 0 : naz + 1;

			const double dv_r =	vr(nr+1, naz) - vr(nr, naz);
			const double dv_phi = vphi(nr,	naz_next) -	vphi(nr, naz);

			const double energy_new =
				energy(nr, naz) -
				dt * Qr(nr, naz) * dv_r * InvDiffRsup[nr] -
				dt * Qphi(nr, naz) * dv_phi * invdxtheta;

			energy(nr, naz) = energy_new;
			}
		}
		}
	}

	// add artificial viscous pressure source term to v_radial
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = One_no_ghost_vr; nr < MaxMo_no_ghost_vr; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {
		// 1/Sigma dq_r/dr : Sigma is calculated as a mean value between
		// the neightbour cells
		vr(nr, naz) = vr(nr, naz) -
			dt * 2.0 / (density(nr, naz) + density(nr-1, naz)) *
			(Qr(nr, naz) - Qr(nr - 1, naz)) * InvDiffRmed[nr];
		}
	}

	// add artificial viscous pressure source term to v_azimuthal
	#pragma omp parallel for
	for (unsigned int nr = Zero_no_ghost; nr < Max_no_ghost; ++nr) {

		const double dxtheta = dphi * Rmed[nr];
		const double invdxtheta = 1.0 / dxtheta;

		for (unsigned int naz = 0; naz < Nphi; ++naz) {
		// 1/Sigma 1/r dq_phi/dphi : Sigma is calculated as a mean value
		// between the neightbour cells
		const unsigned int naz_prev = (naz == 0 ? Nphi-1 : naz - 1);

		vphi(nr, naz) =	vphi(nr, naz) -
			dt * 2.0 / (density(nr, naz) + density(nr, naz_prev)) *
			(Qphi(nr, naz) - Qphi(nr, naz_prev)) * invdxtheta;
		}
	}
	}
}


}
