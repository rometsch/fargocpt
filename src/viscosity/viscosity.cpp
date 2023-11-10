/**
	\file viscosity.cpp

	Calculation of the viscous force.

	The function FViscosity() returns the (kinematic) viscosity as a
   function of the radius (it handles all case: alpha or uniform viscosity, and
   inner cavity with a different viscosity). The update of the velocity is done
   in ViscousTerm(), which properly evaluate the stress tensor in 2D cylindrical
   coordinates.
*/
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>

#include "../global.h"
#include "../parameters.h"
#include "../units.h"
#include "../util.h"
#include "../pvte_law.h"
#include "../quantities.h"
#include "../constants.h"
#include "viscosity.h"
#include <cassert>
#include <cmath>

namespace viscosity
{

double get_alpha(const int nr, const int naz, t_data &data)
{
	switch (parameters::AlphaMode){
	case CONST_ALPHA:
	{
	return parameters::viscous_alpha;
	}
        case SCURVE_ALPHA:
        {
            const double temperatureCGS = data[t_data::TEMPERATURE](nr, naz) * units::temperature;
            const double alpha_cool = parameters::alphaCold*std::pow(Rmed[nr]/0.4, 0.3);
            const double alpha_hot = parameters::alphaHot;
            const double alpha =
                std::pow(10.0, 0.5*(std::log10(alpha_hot)-std::log10(alpha_cool))*
                                       (1.0-std::tanh((4.0-std::log10(temperatureCGS))/0.4)) + std::log10(alpha_cool));
            data[t_data::ALPHA](nr, naz) = alpha;
            return alpha;
	}
	case ALPHA_STAR_DIST_DEPENDEND:
			{

			static const unsigned int N_planets =
			data.get_planetary_system().get_number_of_planets();

			const int cell = get_cell_id(nr, naz);
			const double x = CellCenterX->Field[cell];
			const double y = CellCenterY->Field[cell];

			double alpha = parameters::alphaHot;
			for (unsigned int k = 0; k < N_planets; k++) {

				const double dx = x - g_xpl[k];
				const double dy = y - g_ypl[k];
				const double dist_2 = std::pow(dx, 2) + std::pow(dy, 2);
				const double d = std::sqrt(dist_2);

				const double dist_start = 0.35; // 7.0au
				const double dist_end   = 0.55; // 11.0au
				const double scale_unbound = (d - dist_start)/(dist_end - dist_start);
				const double scale = std::max(0.0, std::min(scale_unbound, 1.0));
				const double alpha_new = parameters::alphaCold + (parameters::alphaHot-parameters::alphaCold)*scale;
				alpha = std::min(alpha, alpha_new);
			}

			return alpha;
			}		
	case SCURVE_IONFRACTION:
	{
			const double temperatureCGS = data[t_data::TEMPERATURE](nr, naz) * units::temperature;
			const double sigma = data[t_data::SIGMA](nr, naz);
			const double scale_height = data[t_data::SCALE_HEIGHT](nr, naz);
			const double densityCGS =
			    sigma / (parameters::density_factor * scale_height) * units::density;
			const double alpha = parameters::alphaCold +
					     (parameters::alphaHot - parameters::alphaCold) *
						 std::min( 1000.0 * pvte::H_ionization_fraction(densityCGS, temperatureCGS), 1.0);
			return alpha;
	}
	}

	return 0.0;
}

/**
	updates nu-grid. If ViscosityAlpha is enabled, soundspeed-grid is needed
   for this
*/
void update_viscosity(t_data &data)
{
    static bool calculated = false;
    // if alpha-viscosity
    if (parameters::viscous_alpha > 0) {
	const unsigned int Nr = data[t_data::VISCOSITY].get_size_radial();
	const unsigned int Nphi = data[t_data::VISCOSITY].get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {
		// H = c_s^iso / Omega_K = c_s_adb / Omega_K / sqrt(gamma)
		// c_s_adb^2 = gamma * c_s_iso^2

		// nu = alpha * c_s_adb * H = alpha * c_s_adb^2 / sqrt(gamma) /
		// Omega_K
		const double alpha = get_alpha(nr, naz, data);
		const double c_s_adb = data[t_data::SOUNDSPEED](nr, naz);
		const double H = data[t_data::SCALE_HEIGHT](nr, naz);
		const double nu = alpha * H * c_s_adb;

		data[t_data::VISCOSITY](nr, naz) = nu;
	    }
	}
    } else {
	if (!calculated) {
		const unsigned int Nr = data[t_data::VISCOSITY].get_size_radial();
		const unsigned int Nphi = data[t_data::VISCOSITY].get_size_azimuthal();

		#pragma omp parallel for collapse(2)
		for (unsigned int nr = 0;	 nr < Nr; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {
			data[t_data::VISCOSITY](nr, naz) = parameters::constant_viscosity;
		}
	    }
	}

	calculated = true;
    }
}

void compute_viscous_stress_tensor(t_data &data)
{

	const unsigned int Nr = data[t_data::DIV_V].get_size_radial();
	const unsigned int Nphi = data[t_data::DIV_V].get_size_azimuthal();

	#pragma omp parallel
	{
    // calculate div(v)
	#pragma omp for collapse(2)
	for (unsigned int nr = 0; nr < Nr ; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const unsigned int naz_next = (naz == Nphi-1 ? 0 : naz + 1);

		// div(v) = 1/r d(r*v_r)/dr + 1/r d(v_phi)/dphi
		data[t_data::DIV_V](nr, naz) =
		(data[t_data::V_RADIAL](nr+1, naz)*Ra[nr+1] -
		 data[t_data::V_RADIAL](nr, naz)*Ra[nr])
				* InvDiffRsupRb[nr] +
		(data[t_data::V_AZIMUTHAL](nr, naz_next) -
		 data[t_data::V_AZIMUTHAL](nr, naz)) * invdphi * InvRb[nr];
	}
    }

    // calculate tau_r_r
	#pragma omp for collapse(2) nowait
	for (unsigned int nr = 0; nr < Nr; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {
	    // d(v_r)/dr (cell centered)
	    const double drr =
		(data[t_data::V_RADIAL](nr + 1, naz) -
		 data[t_data::V_RADIAL](nr, naz)) *
		InvDiffRsup[nr];

		// tau_r_r = 2*nu*Sigma*( d(v_r)/dr - 1/3 div(v))
		data[t_data::TAU_R_R](nr, naz) =
		2.0 * data[t_data::VISCOSITY](nr, naz) *
		data[t_data::SIGMA](nr, naz) *
		(drr - 1.0 / 3.0 * data[t_data::DIV_V](nr, naz));
	}
    }

    // calculate tau_phi_phi
	#pragma omp for collapse(2) nowait
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const unsigned int naz_next = (naz == Nphi-1 ? 0 : naz + 1);

	    // 1/r d(v_phi)/dphi + v_r/r (cell centered)
	    const double dpp =
		(data[t_data::V_AZIMUTHAL](nr, naz_next) -
		 data[t_data::V_AZIMUTHAL](nr, naz)) *
			invdphi * InvRmed[nr] +
		0.5 *
			(data[t_data::V_RADIAL](nr + 1, naz) +
			 data[t_data::V_RADIAL](nr, naz)) *
			InvRmed[nr];

	    // tau_phi_phi = 2*nu*Sigma*( 1/r d(v_phi)/dphi + v_r/r - 1/3 div(v)
	    // )
		const double nu = data[t_data::VISCOSITY](nr, naz);
		const double sigma = data[t_data::SIGMA](nr, naz);
		data[t_data::TAU_PHI_PHI](nr, naz) =
		2.0 * nu * sigma *
		(dpp - 1.0 / 3.0 * data[t_data::DIV_V](nr, naz));

	    const double correction_helper_value = nu * sigma;
		data[t_data::VISCOSITY_SIGMA](nr, naz) =
		correction_helper_value;
	}
    }

    // calculate tau_r_phi
	#pragma omp for collapse(2) nowait
	for (unsigned int nr = 1; nr < Nr; ++nr) { // Nr_vec -1 = Nr
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const double naz_prev = (naz == 0 ? Nphi-1 : naz - 1);

	    // d(v_phi/r)/dr
	    double dvphirdr =
		(data[t_data::V_AZIMUTHAL](nr, naz) * InvRb[nr] -
		 data[t_data::V_AZIMUTHAL](nr - 1, naz) * InvRb[nr - 1]) *
		InvDiffRmed[nr];
	    // d(v_r)/dphi
	    double dvrdphi =
		(data[t_data::V_RADIAL](nr, naz) -
		 data[t_data::V_RADIAL](nr, naz_prev)) *
		invdphi;

	    // r*d(v_phi/r)/dr + 1/r d(v_r)/dphi (edge)
	    const double drp =
		Ra[nr] * dvphirdr + dvrdphi * InvRa[nr];

	    // averaged nu over 4 corresponding cells
	    const double nu =
		0.25 *
		(data[t_data::VISCOSITY](nr, naz) +
		 data[t_data::VISCOSITY](nr - 1, naz) +
		 data[t_data::VISCOSITY](nr, naz_prev) +
		 data[t_data::VISCOSITY](nr - 1, naz_prev));

	    // averaged sigma over 4 corresponding cells
	    const double sigma =
		0.25 * (data[t_data::SIGMA](nr, naz) +
			data[t_data::SIGMA](nr - 1, naz) +
			data[t_data::SIGMA](nr, naz_prev) +
			data[t_data::SIGMA](nr - 1, naz_prev));

	    // tau_r_phi = nu*Sigma*( r*d(v_phi/r)/dr + 1/r d(v_r)/dphi )
		data[t_data::TAU_R_PHI](nr, naz) = nu * sigma * drp;

	    const double correction_helper_value = nu * sigma;
		data[t_data::VISCOSITY_SIGMA_RP](nr, naz) =
		correction_helper_value;
	}
    }

    if (StabilizeViscosity) {
	#pragma omp for collapse(2)
	for (unsigned int nr = 1; nr < Nr; ++nr) { // Nr_vr - 1 = Nr
		for (unsigned int naz = 0; naz < Nphi; ++naz) {

		/// Load general data
		/// ////////////////////////////////////////////////////////////
		const unsigned int naz_prev = (naz == 0 ? Nphi-1 : naz-1);
		const unsigned int naz_next = (naz == Nphi-1 ? 0 : naz+1);

		const double NuSig_rp =
			data[t_data::VISCOSITY_SIGMA_RP](nr, naz);
		const double NuSig_rp_ip =
			data[t_data::VISCOSITY_SIGMA_RP](nr + 1, naz);
		const double NuSig_rp_jp = data[t_data::VISCOSITY_SIGMA_RP](
			nr, naz_next);

		const double NuSigma =
			data[t_data::VISCOSITY_SIGMA](nr, naz);
		const double NuSigma_jm =
			data[t_data::VISCOSITY_SIGMA](nr, naz_prev);
		const double NuSigma_im =
			data[t_data::VISCOSITY_SIGMA](nr - 1, naz);

		/// END Load general data
		/// ////////////////////////////////////////////////////////////

		/// Calc V_phi correction factor
		/// //////////////////////////////////////
		const double Ra3NuSigmaInvDiffRmed = NuSig_rp *
							 std::pow(Ra[nr], 3) *
							 InvDiffRmed[nr];
		const double Ra3NuSigmaInvDiffRmed_p =
			NuSig_rp_ip * std::pow(Ra[nr + 1], 3) *
			InvDiffRmed[nr + 1];

		const double cphi_rp =
			-InvRmed[nr] * TwoDiffRaSq[nr] *
		    (Ra3NuSigmaInvDiffRmed_p + Ra3NuSigmaInvDiffRmed);
		double cphi_pp =
			-FourThirdInvRbInvdphiSq[nr] * (NuSigma + NuSigma_jm);

		const double sigma_avg_phi =
			0.5 * (data[t_data::SIGMA](nr, naz) +
			   data[t_data::SIGMA](nr, naz_prev));

		const double c1_phi =
			(cphi_rp + cphi_pp) / (sigma_avg_phi * Rmed[nr]);
		data[t_data::VISCOSITY_CORRECTION_FACTOR_PHI](
			nr, naz) = c1_phi;

		/// END Calc V_phi correction factor
		/// ////////////////////////////////

		/// Calc V_r correction factor
		/// //////////////////////////////////////
		const double sigma_avg_r =
			0.5 * (data[t_data::SIGMA](nr, naz) +
			   data[t_data::SIGMA](nr - 1, naz));

		const double cr_rp =
			-(NuSig_rp_jp + NuSig_rp) / (dphi * dphi * Ra[nr]);

		double cr_pp_1 = 2.0 * NuSigma *
			(0.5 * InvRmed[nr] + 1.0 / 3.0 * Ra[nr] * InvDiffRsupRb[nr]);
		double cr_pp_2 = 2.0 * NuSigma_im *
			(0.5 * InvRmed[nr - 1] - 1.0 / 3.0 * Ra[nr] * InvDiffRsupRb[nr - 1]);

		double cr_rr_1 =
			Rmed[nr] * 2.0 * NuSigma *
			(-InvDiffRsup[nr] + 1.0 / 3.0 * Ra[nr] * InvDiffRsupRb[nr]);
		double cr_rr_2 =
			-1.0 * Rmed[nr - 1] * 2.0 * NuSigma_im *
			(InvDiffRsup[nr - 1] - 1.0 / 3.0 * Ra[nr] * InvDiffRsupRb[nr - 1]);

		const double cr_pp = -0.5 * (cr_pp_1 + cr_pp_2);
		const double cr_rr = InvDiffRmed[nr] * (cr_rr_1 + cr_rr_2);

		const double Rmed_mid = 0.5 * (Rb[nr] + Rb[nr - 1]);
		const double c1_r = parameters::radial_viscosity_factor *
				    (cr_rr + cr_rp + cr_pp) /
				    (sigma_avg_r * Rmed_mid);

		assert(c1_r < 0.0);
		assert(c1_phi < 0.0);

		data[t_data::VISCOSITY_CORRECTION_FACTOR_R](nr,
								naz) = c1_r;
		/// END Calc V_r correction factor
		/// //////////////////////////////////////
	    }
	}
    }
	}
}

/**
	Update velocities with viscous source term of Navier-Stokes equations
*/
void update_velocities_with_viscosity(t_data &data, const double dt)
{

	t_polargrid &v_azimuthal = data[t_data::V_AZIMUTHAL];
	t_polargrid &v_radial = data[t_data::V_RADIAL];
    const t_polargrid &Sigma = data[t_data::SIGMA];
    const t_polargrid &Trp = data[t_data::TAU_R_PHI];
    const t_polargrid &Trr = data[t_data::TAU_R_R];
    const t_polargrid &Tpp = data[t_data::TAU_PHI_PHI];

	const unsigned int Nr = v_azimuthal.get_size_radial();
	const unsigned int Nphi = v_azimuthal.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nr-1; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const int naz_prev = (naz == 0 ? Nphi-1 : naz - 1);

		const double sigma_avg = 0.5 * (Sigma(nr, naz) + Sigma(nr, naz_prev));

	    // See D'Angelo et al. 2002 Nested-grid calculations of disk-planet
	    // interaction It is important to use the conservative form here and
	    // not the one from Fargo / Baruteau. a_phi = 1/(r*Sigma) ( 1/r
	    // d(r^2 * tau_r_phi)/dr + d(tau_phi_phi)/dphi )
	    double dVp =
		dt * InvRb[nr] / (sigma_avg) *
		((2.0 / (std::pow(Ra[nr + 1], 2) - std::pow(Ra[nr], 2))) *
		     (std::pow(Ra[nr + 1], 2) * Trp(nr + 1, naz) -
		      std::pow(Ra[nr], 2) * Trp(nr, naz)) +
		 (Tpp(nr, naz) - Tpp(nr, naz_prev)) * invdphi);

	    if (StabilizeViscosity == 1) {
		const double cphi =
		    data[t_data::VISCOSITY_CORRECTION_FACTOR_PHI](nr, naz);
		const double corr = 1.0 / (std::max(1.0 + dt * cphi, 0.0) - dt * cphi);
		dVp *= corr;
	    }

	    v_azimuthal(nr, naz) += dVp;
	}}

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = One_no_ghost_vr; nr < MaxMo_no_ghost_vr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const int naz_next = (naz == Nphi-1 ? 0 : naz + 1);

	    // a_r = 1/(r*Sigma) ( d(r*tau_r_r)/dr + d(tau_r_phi)/dphi -
	    // tau_phi_phi )
		const double sigma_avg = 0.5 * (Sigma(nr, naz) + Sigma(nr - 1, naz));

	    double dVr =
		dt / (sigma_avg)*parameters::radial_viscosity_factor * 2.0 /
		(Rb[nr] + Rb[nr - 1]) *
		((Rb[nr] * Trr(nr, naz) - Rb[nr - 1] * Trr(nr - 1, naz)) *
		     InvDiffRmed[nr] +
		 (Trp(nr, naz_next) - Trp(nr, naz)) * invdphi -
		 0.5 * (Tpp(nr, naz) + Tpp(nr - 1, naz)));

	    if (StabilizeViscosity == 1) {
		const double cr = data[t_data::VISCOSITY_CORRECTION_FACTOR_R](nr, naz);
		const double corr = 1.0 / (std::max(1.0 + dt * cr, 0.0) - dt * cr);
		dVr *= corr;
	    }

	    v_radial(nr, naz) += dVr;
	}
    }

	if(ECC_GROWTH_MONITOR){
		quantities::calculate_disk_delta_ecc_peri(data, delta_ecc_visc, delta_peri_visc);
	}
}

} // namespace viscosity
