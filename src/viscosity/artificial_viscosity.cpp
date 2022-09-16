#include "artificial_viscosity.h"

#include "../parameters.h"
#include "../SideEuler.h"
#include "../global.h"
#include "../SourceEuler.h"
#include "../quantities.h"
#include <cassert>

namespace art_visc{

void update_with_artificial_viscosity(t_data &data, const double time, const double dt){
	if (parameters::artificial_viscosity ==
	 parameters::artificial_viscosity_WT) {
	recalculate_viscosity(data, time);
	art_visc::update_with_artificial_viscosity_TW(data, dt);
	} else { // SN or TW (just for Vazimuthal ghost cells)
	art_visc::update_with_artificial_viscosity_SN(data, dt);
	recalculate_viscosity(data, time);
	}

	if(ECC_GROWTH_MONITOR){
		quantities::calculate_disk_delta_ecc_peri(data, delta_ecc_art_visc, delta_peri_art_visc);
	}
}


void update_with_artificial_viscosity_TW(t_data &data, const double dt)
{

	/// Do not Apply sub keplerian boundary for boundary conditions that set
	/// Vphi themselves
	const bool add_kep_inner =
	(parameters::boundary_inner !=
	 parameters::boundary_condition_evanescent) &&
	(parameters::boundary_inner !=
	 parameters::boundary_condition_boundary_layer) &&
	(parameters::boundary_inner !=
	 parameters::boundary_condition_precribed_time_variable) &&
	(!parameters::domegadr_zero);

	if (add_kep_inner) {
	ApplySubKeplerianBoundaryInner(data[t_data::V_AZIMUTHAL]);
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

	t_polargrid &DIV_V = data[t_data::DIV_V];
	t_polargrid &SIGMA_ART_NU = data[t_data::SIGMA_ART_VISC];

	//t_polargrid &Trr = data[t_data::TAU_R_R];
	//t_polargrid &Tpp = data[t_data::TAU_PHI_PHI];

	t_polargrid &Tau_art = data[t_data::TAU_R_R];

	t_polargrid &Qp = data[t_data::QPLUS];
	t_polargrid &visc = data[t_data::VISCOSITY];

	const unsigned int Nr = DIV_V.get_size_radial();
	const unsigned int Nphi = DIV_V.get_size_azimuthal();

	// calculate div(v)
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		// div(v) = 1/r d(r*v_r)/dr + 1/r d(v_phi)/dphi
		const double naz_next = naz == vphi.get_max_azimuthal() ? 0 : naz + 1;

		DIV_V(nr, naz) = (vr(nr + 1, naz) * Ra[nr + 1] -
		 vr(nr, naz) * Ra[nr]) * InvDiffRsup[nr] * InvRb[nr] +
		(vphi(nr, naz_next) - vphi(nr, naz)) * invdphi * InvRb[nr];

		const double Dr = Ra[nr+1] - Ra[nr];
		const double rDphi = Rmed[nr] * dphi;
		const double dx_2 = std::pow(std::min(Dr, rDphi), 2);

		double sigma_nu_art;
		if (DIV_V(nr, naz) < 0) {
			sigma_nu_art = parameters::artificial_viscosity_factor *
				 density(nr, naz) * dx_2 * (-DIV_V(nr, naz));
		} else {
			sigma_nu_art = 0;
		}

		Tau_art(nr, naz) = sigma_nu_art * DIV_V(nr, naz);
		/*Trr(nr, naz) =
			nu_art * DIV_V(nr, naz);
		Tpp(nr, naz) =
			nu_art * DIV_V(nr, naz);*/
		SIGMA_ART_NU(nr, naz) = sigma_nu_art;

		if (parameters::Adiabatic && parameters::artificial_viscosity_dissipation) {
		if(nr > 1 && nr < density.get_max_radial()){
			if(parameters::ALPHAVISCOSITY > 0){
				if(parameters::ALPHAVISCOSITY > 0.0 && (nr > 0)){

					double qplus =
							1.0 / (2.0 * visc(nr, naz) *
								   density(nr, naz)) * 2.0 * std::pow(Tau_art(nr, naz), 2);
					// Qplus ~ 1 / energy
					// dE/dt = const / E
					// solve analytically:
					// -> E(t) = sqrt(2) * sqrt(c * t + k)
					const double c = qplus * energy(nr, naz); // qplus ~ 1 / E -> qplus * E = const
					const double energy_old = energy(nr, naz);
					const double k = 0.5 * std::pow(energy_old, 2);
					const double energy_new = std::sqrt(2.0) * std::sqrt(c * dt + k);

					//const double de_dt = (energy_new - energy_old)/dt;
					//Qp(nr, naz) = de_dt;
					energy(nr, naz) = energy_new;
				}

			} else {
			if (visc(nr, naz) != 0.0) {

				double qplus =
						1.0 / (2.0 * visc(nr, naz) *
							   density(nr, naz)) * 2.0 * std::pow(Tau_art(nr, naz), 2);
				/*double qplus =
					1.0 / (2.0 * visc(nr, naz) *
					 density(nr, naz)) *
					(std::pow(Trr(nr, naz), 2)
					 + std::pow(Tpp(nr, naz), 2));*/

				qplus *= parameters::heating_viscous_factor;
				Qp(nr, naz) = qplus;
				//energy(nr, naz) += dt*qplus;
			}
			}
		}
		}

	}
	}

	if (StabilizeArtViscosity > 0) {
	const unsigned int Nrv = vr.get_max_radial(); // = (Nr+1) - 1 (vr is a vector)
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nrv; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {

		unsigned int naz_prev = (naz == 0 ? Nphi - 1 : naz - 1);

		/// Calc V_phi correction factor
		/// //////////////////////////////////////

		const double SigNuArt = SIGMA_ART_NU(nr, naz);
		const double SigNuArt_jm = SIGMA_ART_NU(nr, naz_prev);

		const double cphi_pp = -(SigNuArt_jm + SigNuArt) *
				(invdphi * invdphi * InvRmed[nr]);

		const double sigma_phi_avg = 0.5 *
				(density(nr, naz) + density(nr, naz_prev));

		const double c1_phi =
			cphi_pp / (sigma_phi_avg * Rmed[nr]);
		data[t_data::ART_VISCOSITY_CORRECTION_FACTOR_PHI](nr, naz) = c1_phi;

		/// END Calc V_phi correction factor
		/// ////////////////////////////////

		/// Calc V_r correction factor
		/// //////////////////////////////////////
		const double sigma_r_avg =
			0.5 * (density(nr, naz) +
			   density(nr - 1, naz));


		const double SigNuArt_im = SIGMA_ART_NU(nr - 1, naz);

		const double cr_pp_1 = SigNuArt * Ra[nr] * InvDiffRsupRb[nr];
		const double cr_pp_2 = -SigNuArt_im * Ra[nr] * InvDiffRsupRb[nr - 1];

		const double cr_rr_1 = -SigNuArt * Ra[nr] * InvDiffRsup[nr];
		const double cr_rr_2 = -SigNuArt_im * Ra[nr] * InvDiffRsup[nr - 1];

		const double cr_pp = -0.5 * (cr_pp_1 + cr_pp_2);
		const double cr_rr = InvDiffRmed[nr] * (cr_rr_1 + cr_rr_2);

		const double Rmed_mid = 0.5 * (Rb[nr] + Rb[nr - 1]);
		const double c1_r = parameters::radial_viscosity_factor *
					(cr_rr + cr_pp) / (sigma_r_avg * Rmed_mid);

		assert(c1_r < 0.0);
		assert(c1_phi < 0.0);

		data[t_data::ART_VISCOSITY_CORRECTION_FACTOR_R](nr,	naz) = c1_r;
		/// END Calc V_r correction factor
		/// //////////////////////////////////////
		}
	}
	}


	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {

		const int naz_prev = (naz == 0 ? Nphi-1 : naz - 1);

		const double sigma_phi_avg = 0.5 * (density(nr, naz) + density(nr, naz_prev));

		// See D'Angelo et al. 2002 Nested-grid calculations of disk-planet
		// interaction It is important to use the conservative form here and
		// not the one from Fargo / Baruteau. a_phi = 1/(r*Sigma) ( 1/r
		// d(r^2 * tau_r_phi)/dr + d(tau_phi_phi)/dphi )
		/*double dVp =
		dt * InvRb[nr] / (sigma_phi_avg) *
		(Tpp(nr, naz) - Tpp(nr, naz_minus)) * invdphi;*/

		double dVp =
		dt * InvRb[nr] / (sigma_phi_avg) *
		(Tau_art(nr, naz) - Tau_art(nr, naz_prev)) * invdphi;


		if (StabilizeArtViscosity == 1) {
		const double cphi =
			data[t_data::ART_VISCOSITY_CORRECTION_FACTOR_PHI](nr, naz);
		const double corr =
			1.0 / (std::max(1.0 + dt * cphi, 0.0) - dt * cphi);
		dVp *= corr;
		}

		vphi(nr, naz) += dVp;

		// a_r = 1/(r*Sigma) ( d(r*tau_r_r)/dr + d(tau_r_phi)/dphi -
		// tau_phi_phi )
		const double sigma_r_avg = 0.5 * (density(nr, naz) + density(nr - 1, naz));

		/*double dVr =
		dt / (sigma_r_avg)*parameters::radial_viscosity_factor * 2.0 /
		(Rb[nr] + Rb[nr - 1]) *
		((Rb[nr] * Trr(nr, naz) - Rb[nr - 1] * Trr(nr - 1, naz)) *
			 InvDiffRmed[nr]
				- 0.5 * (Tpp(nr, naz) + Tpp(nr - 1, naz)));*/
		double dVr =
				dt / (sigma_r_avg)*parameters::radial_viscosity_factor * 2.0 /
				(Rb[nr] + Rb[nr - 1]) *
				((Rb[nr] * Tau_art(nr, naz) - Rb[nr - 1] * Tau_art(nr - 1, naz)) *
					 InvDiffRmed[nr]
						- 0.5 * (Tau_art(nr, naz) + Tau_art(nr - 1, naz)));

		if (StabilizeArtViscosity == 1) {
		const double cr =
			data[t_data::ART_VISCOSITY_CORRECTION_FACTOR_R](nr, naz);
		const double corr =
			1.0 / (std::max(1.0 + dt * cr, 0.0) - dt * cr);
		dVr *= corr;
		}

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
	 parameters::boundary_condition_evanescent) &&
	(parameters::boundary_inner !=
	 parameters::boundary_condition_boundary_layer) &&
	(parameters::boundary_inner !=
	 parameters::boundary_condition_precribed_time_variable) &&
	(!parameters::domegadr_zero);

	if (add_kep_inner) {
	ApplySubKeplerianBoundaryInner(data[t_data::V_AZIMUTHAL]);
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
		parameters::artificial_viscosity_SN &&
	parameters::EXPLICIT_VISCOSITY) {

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
		for (unsigned int nr = 0; nr < Nr; ++nr) {

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
	for (unsigned int nr = 1; nr < Nr; ++nr) { // Nr = Nr_vr-1, so were leaving out the ghost cells
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
	for (unsigned int nr = 0; nr < Nr; ++nr) {

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
