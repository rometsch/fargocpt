/**
	\file viscosity.cpp

	Calculation of the viscous force.

	The function FViscosity() returns the (kinematic) viscosity as a
   function of the radius (it handles all case: alpha or uniform viscosity, and
   inner cavity with a different viscosity). The update of the velocity is done
   in ViscousTerm(), which properly evaluate the stress tensor in 2D cylindrical
   coordinates.
*/

#include <math.h>

#include "LowTasks.h"
#include "Theo.h"
#include "axilib.h"
#include "constants.h"
#include "global.h"
#include "output.h"
#include "parameters.h"
#include "units.h"
#include "util.h"
#include "viscosity.h"
#include <cmath>
#include <cassert>

namespace viscosity
{

/**
	updates nu-grid. If ViscosityAlpha is enabled, soundspeed-grid is needed
   for this
*/
void update_viscosity(t_data &data)
{
    static bool calculated = false;
    // if alpha-viscosity
    if (ViscosityAlpha) {
	for (unsigned int n_rad = 0;
		 n_rad <= data[t_data::VISCOSITY].get_max_radial(); ++n_rad) {
	    const double inv_omega_kepler =
		1.0 / calculate_omega_kepler(Rb[n_rad]);

		for (unsigned int n_az = 0;
		 n_az <= data[t_data::VISCOSITY].get_max_azimuthal();
		 ++n_az) {
			// H = c_s^iso / Omega_K = c_s_adb / Omega_K / sqrt(gamma)
			// c_s_adb^2 = gamma * c_s_iso

			// nu = alpha * c_s_adb * H = alpha * c_s_adb^2 / sqrt(gamma) /
			// Omega_K
			const double alpha = ALPHAVISCOSITY;
			const double gamma = ADIABATICINDEX;
			const double c_s_adb = data[t_data::SOUNDSPEED](n_rad, n_az);
			const double nu = alpha * std::pow(c_s_adb, 2) *
					  inv_omega_kepler / std::sqrt(gamma);

			data[t_data::VISCOSITY](n_rad, n_az) = nu;
	    }
	}
    } else {
	if (!calculated) {
	    for (unsigned int n_radial = 0;
		 n_radial <= data[t_data::VISCOSITY].get_max_radial();
		 ++n_radial) {
		for (unsigned int n_azimuthal = 0;
		     n_azimuthal <= data[t_data::VISCOSITY].get_max_azimuthal();
		     ++n_azimuthal) {
		    data[t_data::VISCOSITY](n_radial, n_azimuthal) = VISCOSITY;
		}
	    }
	}

	calculated = true;
    }
}

void compute_viscous_terms(t_data &data, bool include_artifical_viscosity)
{

    // calculate div(v)
    for (unsigned int n_radial = 0;
	 n_radial <= data[t_data::DIV_V].get_max_radial(); ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::DIV_V].get_max_azimuthal();
	     ++n_azimuthal) {
	    // div(v) = 1/r d(r*v_r)/dr + 1/r d(v_phi)/dphi
	    data[t_data::DIV_V](n_radial, n_azimuthal) =
		(data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) *
		     Ra[n_radial + 1] -
		 data[t_data::V_RADIAL](n_radial, n_azimuthal) * Ra[n_radial]) *
		    InvDiffRsup[n_radial] * InvRb[n_radial] +
		(data[t_data::V_AZIMUTHAL](
		     n_radial,
		     n_azimuthal == data[t_data::DIV_V].get_max_azimuthal()
			 ? 0
			 : n_azimuthal + 1) -
		 data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal)) *
		    invdphi * InvRb[n_radial];
	}
    }

    // calculate tau_r_r
    for (unsigned int n_radial = 0;
	 n_radial <= data[t_data::TAU_R_R].get_max_radial(); ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::TAU_R_R].get_max_azimuthal();
	     ++n_azimuthal) {
	    // d(v_r)/dr (cell centered)
		const double drr = (data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) -
		   data[t_data::V_RADIAL](n_radial, n_azimuthal)) *
		  InvDiffRsup[n_radial];

	    // tau_r_r = 2*nu*Sigma*( d(v_r)/dr - 1/3 div(v) + eta2 div(v))
	    data[t_data::TAU_R_R](n_radial, n_azimuthal) =
		2.0 * data[t_data::VISCOSITY](n_radial, n_azimuthal) *
		data[t_data::DENSITY](n_radial, n_azimuthal) *
		(drr - 1.0 / 3.0 * data[t_data::DIV_V](n_radial, n_azimuthal));
	}
    }

    // calculate tau_phi_phi
    for (unsigned int n_radial = 0;
	 n_radial <= data[t_data::TAU_PHI_PHI].get_max_radial(); ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::TAU_PHI_PHI].get_max_azimuthal();
	     ++n_azimuthal) {
	    // 1/r d(v_phi)/dphi + v_r/r (cell centered)
		const double dpp = (data[t_data::V_AZIMUTHAL](
		       n_radial,
		       n_azimuthal ==
			       data[t_data::TAU_PHI_PHI].get_max_azimuthal()
			   ? 0
			   : n_azimuthal + 1) -
		   data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal)) *
		      invdphi * InvRmed[n_radial] +
		  0.5 *
		      (data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) +
		       data[t_data::V_RADIAL](n_radial, n_azimuthal)) *
		      InvRmed[n_radial];

	    // tau_phi_phi = 2*nu*Sigma*( 1/r d(v_phi)/dphi + v_r/r - 1/3 div(v)
	    // )
		const double nu = data[t_data::VISCOSITY](n_radial, n_azimuthal);
		const double sigma = data[t_data::DENSITY](n_radial, n_azimuthal);
	    data[t_data::TAU_PHI_PHI](n_radial, n_azimuthal) =
		2.0 * nu * sigma * (dpp - 1.0 / 3.0 * data[t_data::DIV_V](n_radial, n_azimuthal));

		const double correction_helper_value = nu * sigma;
		data[t_data::VISCOSITY_SIGMA](n_radial, n_azimuthal) = correction_helper_value;
	}
    }

    // calculate tau_r_phi
    for (unsigned int n_radial = 1;
	 n_radial <= data[t_data::TAU_R_PHI].get_max_radial() - 1; ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::TAU_R_PHI].get_max_azimuthal();
	     ++n_azimuthal) {
	    // d(v_phi/r)/dr
	    double dvphirdr =
		(data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) *
		     InvRb[n_radial] -
		 data[t_data::V_AZIMUTHAL](n_radial - 1, n_azimuthal) *
		     InvRb[n_radial - 1]) *
		InvDiffRmed[n_radial];
	    // d(v_r)/dphi
	    double dvrdphi =
		(data[t_data::V_RADIAL](n_radial, n_azimuthal) -
		 data[t_data::V_RADIAL](
		     n_radial, n_azimuthal == 0
				   ? data[t_data::V_RADIAL].get_max_azimuthal()
				   : n_azimuthal - 1)) *
		invdphi;

	    // r*d(v_phi/r)/dr + 1/r d(v_r)/dphi (edge)
		const double drp = Ra[n_radial] * dvphirdr + dvrdphi * InvRa[n_radial];

	    unsigned int n_azimuthal_minus =
		(n_azimuthal == 0 ? data[t_data::VISCOSITY].get_max_azimuthal()
				  : n_azimuthal - 1);

	    // averaged nu over 4 corresponding cells
		const double nu =
		0.25 *
		(data[t_data::VISCOSITY](n_radial, n_azimuthal) +
		 data[t_data::VISCOSITY](n_radial - 1, n_azimuthal) +
		 data[t_data::VISCOSITY](n_radial, n_azimuthal_minus) +
		 data[t_data::VISCOSITY](n_radial - 1, n_azimuthal_minus));

	    // averaged sigma over 4 corresponding cells
		const double sigma =
		0.25 * (data[t_data::DENSITY](n_radial, n_azimuthal) +
			data[t_data::DENSITY](n_radial - 1, n_azimuthal) +
			data[t_data::DENSITY](n_radial, n_azimuthal_minus) +
			data[t_data::DENSITY](n_radial - 1, n_azimuthal_minus));

	    // tau_r_phi = nu*Sigma*( r*d(v_phi/r)/dr + 1/r d(v_r)/dphi )
	    data[t_data::TAU_R_PHI](n_radial, n_azimuthal) = nu * sigma * drp;

		const double correction_helper_value = nu * sigma;
		data[t_data::VISCOSITY_SIGMA_RP](n_radial, n_azimuthal) = correction_helper_value;
	}
    }

    if (include_artifical_viscosity) {
	for (unsigned int n_radial = 0;
	     n_radial <= data[t_data::DIV_V].get_max_radial(); ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::DIV_V].get_max_azimuthal();
		 ++n_azimuthal) {
		double nu_art;

		const double dx_2 = std::pow(std::min(
									Rsup[n_radial] - Rinf[n_radial],
									Rmed[n_radial] * dphi), 2);
		if (data[t_data::DIV_V](n_radial, n_azimuthal) < 0) {
		    nu_art =
			parameters::artificial_viscosity_factor *
			data[t_data::DENSITY](n_radial, n_azimuthal) *
					dx_2 *
			(-data[t_data::DIV_V](n_radial, n_azimuthal));
		} else {
		    nu_art = 0;
		}

		data[t_data::TAU_R_R](n_radial, n_azimuthal) +=
		    nu_art * data[t_data::DIV_V](n_radial, n_azimuthal);
		data[t_data::TAU_PHI_PHI](n_radial, n_azimuthal) +=
			nu_art * data[t_data::DIV_V](n_radial, n_azimuthal);
		data[t_data::ARTIFICIAL_VISCOSITY](n_radial, n_azimuthal) = nu_art;
	    }
	}
	}

	if (StabilizeViscosity) {
	for (unsigned int n_radial = 1; n_radial < data[t_data::V_RADIAL].get_max_radial();
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::V_RADIAL].get_size_azimuthal(); ++n_azimuthal){


		/// Load general data	////////////////////////////////////////////////////////////
		unsigned int n_azimuthal_minus =
			(n_azimuthal == 0
			 ? data[t_data::VISCOSITY].get_max_azimuthal()
			 : n_azimuthal - 1);
		const int n_azimuthal_plus =
			n_azimuthal == data[t_data::TAU_PHI_PHI].get_max_azimuthal()
			? 0
			: n_azimuthal + 1;

		const double NuSig_rp = data[t_data::VISCOSITY_SIGMA_RP](n_radial, n_azimuthal);
		const double NuSig_rp_ip = data[t_data::VISCOSITY_SIGMA_RP](n_radial+1, n_azimuthal);
		const double NuSig_rp_jp = data[t_data::VISCOSITY_SIGMA_RP](n_radial, n_azimuthal_plus);

		const double NuSigma = data[t_data::VISCOSITY_SIGMA](n_radial, n_azimuthal);
		const double NuSigma_jm = data[t_data::VISCOSITY_SIGMA](n_radial, n_azimuthal_minus);
		const double NuSigma_im = data[t_data::VISCOSITY_SIGMA](n_radial-1, n_azimuthal);

		/// END Load general data	////////////////////////////////////////////////////////////

		/// Calc V_phi correction factor	//////////////////////////////////////
		const double Ra3NuSigmaInvDiffRmed = NuSig_rp * std::pow(Ra[n_radial], 3) * InvDiffRmed[n_radial];
		const double Ra3NuSigmaInvDiffRmed_p = NuSig_rp_ip * std::pow(Ra[n_radial+1], 3) * InvDiffRmed[n_radial+1];

		const double cphi_rp = - InvRmed[n_radial] * TwoDiffRaSq[n_radial] * (Ra3NuSigmaInvDiffRmed_p + Ra3NuSigmaInvDiffRmed);
		double cphi_pp = - FourThirdInvRbInvdphiSq[n_radial] * (NuSigma + NuSigma_jm);

		if(include_artifical_viscosity){
			const double NuArt = data[t_data::ARTIFICIAL_VISCOSITY](n_radial, n_azimuthal);
			const double NuArt_jm = data[t_data::ARTIFICIAL_VISCOSITY](n_radial, n_azimuthal_minus);
			cphi_pp += (NuArt + NuArt_jm) * (-invdphi * invdphi * InvRmed[n_radial]);
		}

		const double sigma_avg_phi =
		0.5 * (data[t_data::DENSITY](n_radial, n_azimuthal) +
			   data[t_data::DENSITY](n_radial, n_azimuthal_minus));

		const double c1_phi = (cphi_rp + cphi_pp) /(sigma_avg_phi * Rmed[n_radial]);
		data[t_data::VISCOSITY_CORRECTION_FACTOR_PHI](n_radial, n_azimuthal) = c1_phi;

		/// END Calc V_phi correction factor	////////////////////////////////


		/// Calc V_r correction factor	//////////////////////////////////////
		const double sigma_avg_r =
				0.5 * (data[t_data::DENSITY](n_radial, n_azimuthal) +
				data[t_data::DENSITY](n_radial - 1, n_azimuthal));

		const double cr_rp = -(NuSig_rp_jp + NuSig_rp) / (dphi*dphi * Ra[n_radial]);

		double cr_pp_1 = 2.0*NuSigma * (0.5*InvRmed[n_radial] + 1.0/3.0 * Ra[n_radial]*InvDiffRsupRb[n_radial]);
		double cr_pp_2 = 2.0*NuSigma_im * (0.5*InvRmed[n_radial-1] - 1.0/3.0 * Ra[n_radial]*InvDiffRsupRb[n_radial-1]);

		double cr_rr_1 = Rmed[n_radial] * 2.0 * NuSigma * (-InvDiffRsup[n_radial] + 1.0/3.0 * Ra[n_radial]*InvDiffRsupRb[n_radial]);
		double cr_rr_2 = -1.0 * Rmed[n_radial-1] * 2.0 * NuSigma_im * (InvDiffRsup[n_radial-1] - 1.0/3.0 * Ra[n_radial]*InvDiffRsupRb[n_radial-1]);

		if(include_artifical_viscosity){
		const double NuArt = data[t_data::ARTIFICIAL_VISCOSITY](n_radial, n_azimuthal);
		const double NuArt_im = data[t_data::ARTIFICIAL_VISCOSITY](n_radial-1, n_azimuthal);

		cr_pp_1 -= NuArt * Ra[n_radial]*InvDiffRsupRb[n_radial];
		cr_pp_2 += NuArt_im * Ra[n_radial]*InvDiffRsupRb[n_radial-1];

		cr_rr_1 -= NuArt*Ra[n_radial]*InvDiffRsup[n_radial];
		cr_rr_2 -= NuArt_im*Ra[n_radial]*InvDiffRsup[n_radial-1];
		}

		const double cr_pp = -0.5*(cr_pp_1 + cr_pp_2);
		const double cr_rr = InvDiffRmed[n_radial] * (cr_rr_1 + cr_rr_2);

		const double c1_r = parameters::radial_viscosity_factor * (cr_rr + cr_rp + cr_pp) / (sigma_avg_r * Rinf[n_radial]);

		assert(c1_r < 0.0);
		assert(c1_phi < 0.0);

		data[t_data::VISCOSITY_CORRECTION_FACTOR_R](n_radial, n_azimuthal) = c1_r;
		/// END Calc V_r correction factor	//////////////////////////////////////
	}}}
}

/**
	Update velocities with viscous source term of Navier-Stokes equations
*/
void update_velocities_with_viscosity(t_data &data, t_polargrid &v_radial,
				      t_polargrid &v_azimuthal, double dt)
{

    for (unsigned int n_radial = 1; n_radial <= v_radial.get_max_radial() - 1;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= v_radial.get_max_azimuthal(); ++n_azimuthal) {
		const int n_azimuthal_plus =
		(n_azimuthal == data[t_data::DENSITY].get_max_azimuthal()
		     ? 0
		     : n_azimuthal + 1);
		const int n_azimuthal_minus =
		(n_azimuthal == 0 ? data[t_data::DENSITY].get_max_azimuthal()
				  : n_azimuthal - 1);

	    double sigma_avg =
		0.5 * (data[t_data::DENSITY](n_radial, n_azimuthal) +
		       data[t_data::DENSITY](n_radial, n_azimuthal_minus));

		// a_phi = 1/(r*Sigma) ( d(r*tau_r_phi)/dr + d(tau_phi_phi)/dphi + tau_r_phi )
		double dVp =
		dt * InvRb[n_radial] / (sigma_avg) *
		((2.0 / (std::pow(Ra[n_radial + 1], 2) - std::pow(Ra[n_radial], 2))) *
			 (std::pow(Ra[n_radial + 1], 2) *
			  data[t_data::TAU_R_PHI](n_radial + 1, n_azimuthal) -
			  std::pow(Ra[n_radial], 2) *
			  data[t_data::TAU_R_PHI](n_radial, n_azimuthal)) +
		 (data[t_data::TAU_PHI_PHI](n_radial, n_azimuthal) -
		  data[t_data::TAU_PHI_PHI](n_radial, n_azimuthal_minus)) *
			 invdphi);


		if(StabilizeViscosity == 1){
			const double cphi = data[t_data::VISCOSITY_CORRECTION_FACTOR_PHI](n_radial, n_azimuthal);
			const double corr = 1.0/(std::max(1.0 + dt*cphi, 0.0) -dt*cphi);
			dVp *= corr;
		}

		v_azimuthal(n_radial, n_azimuthal) += dVp;



	    // a_r = 1/(r*Sigma) ( d(r*tau_r_r)/dr + d(tau_r_phi)/dphi -
	    // tau_phi_phi )
	    sigma_avg =
		0.5 * (data[t_data::DENSITY](n_radial, n_azimuthal) +
		       data[t_data::DENSITY](n_radial - 1, n_azimuthal));

		double dVr = dt * InvRinf[n_radial] /
				(sigma_avg)*parameters::radial_viscosity_factor *
				((Rmed[n_radial] * data[t_data::TAU_R_R](n_radial, n_azimuthal) -
				  Rmed[n_radial - 1] *
					  data[t_data::TAU_R_R](n_radial - 1, n_azimuthal)) *
					 InvDiffRmed[n_radial] +
				 (data[t_data::TAU_R_PHI](n_radial, n_azimuthal_plus) -
				  data[t_data::TAU_R_PHI](n_radial, n_azimuthal)) *
					 invdphi -
				 0.5 * (data[t_data::TAU_PHI_PHI](n_radial, n_azimuthal) +
					data[t_data::TAU_PHI_PHI](n_radial - 1, n_azimuthal)));

		if(StabilizeViscosity == 1){
			const double cr = data[t_data::VISCOSITY_CORRECTION_FACTOR_R](n_radial, n_azimuthal);
			const double corr = 1.0/(std::max(1.0 + dt*cr, 0.0) -dt*cr);
			dVr *= corr;
		}

		v_radial(n_radial, n_azimuthal) += dVr;
	}
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
/// Only debug functionality below
///////////////////////////////////////////////////////////////////////////////////////////////////////

///
/// \brief get_tau_rp
/// Split up the Tau_r_phi in two components such that Tau_r_phi = Trp1*vphi + Trp2
///
static void get_tau_rp(t_data &data, double &tau_rp_1, double &tau_rp_2, const int n_radial, const int n_azimuthal, bool r_id_p)
{

	if(n_radial == data[t_data::TAU_R_PHI].get_max_radial()){
		tau_rp_1 = 0.0;
		tau_rp_2 = 0.0;
		return;
	}

	int n_azimuthal_minus =
	(n_azimuthal == 0 ? data[t_data::DENSITY].get_max_azimuthal()
			  : n_azimuthal - 1);

	double dvphirdr =
	(data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) *
		 InvRb[n_radial] -
	 data[t_data::V_AZIMUTHAL](n_radial - 1, n_azimuthal) *
		 InvRb[n_radial - 1]) *
	InvDiffRmed[n_radial];
	// d(v_r)/dphi
	double dvrdphi =
	(data[t_data::V_RADIAL](n_radial, n_azimuthal) -
	 data[t_data::V_RADIAL](
		 n_radial, n_azimuthal == 0
			   ? data[t_data::V_RADIAL].get_max_azimuthal()
			   : n_azimuthal - 1)) *
	invdphi;

	// r*d(v_phi/r)/dr + 1/r d(v_r)/dphi (edge)
	const double drp_org = Ra[n_radial] * dvphirdr + dvrdphi * InvRa[n_radial];

	double drp_1 = Ra[n_radial] * (data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) *
			InvRb[n_radial]) *
	   InvDiffRmed[n_radial];

	double drp_2 = Ra[n_radial] * (-
		data[t_data::V_AZIMUTHAL](n_radial - 1, n_azimuthal) *
			InvRb[n_radial - 1]) *
	   InvDiffRmed[n_radial];

	double drp_r = dvrdphi * InvRa[n_radial];

	double drp = drp_1 + drp_2 + drp_r;

	if(drp != 0.0){
	if(std::fabs((drp - drp_org)/drp) > 1e-5){
		printf("drp (%d %d) failed with upd = %.5e	alt = %.5e	rel = %.5e\n", n_radial, n_azimuthal, drp_org, drp, std::fabs(drp_org - drp)/drp_org);
	}}


	// averaged nu over 4 corresponding cells
	const double nu =
	0.25 *
	(data[t_data::VISCOSITY](n_radial, n_azimuthal) +
	 data[t_data::VISCOSITY](n_radial - 1, n_azimuthal) +
	 data[t_data::VISCOSITY](n_radial, n_azimuthal_minus) +
	 data[t_data::VISCOSITY](n_radial - 1, n_azimuthal_minus));

	// averaged sigma over 4 corresponding cells
	const double sigma =
	0.25 * (data[t_data::DENSITY](n_radial, n_azimuthal) +
		data[t_data::DENSITY](n_radial - 1, n_azimuthal) +
		data[t_data::DENSITY](n_radial, n_azimuthal_minus) +
		data[t_data::DENSITY](n_radial - 1, n_azimuthal_minus));

	if(!r_id_p){

	drp_1 = Ra[n_radial] * InvRb[n_radial] * InvDiffRmed[n_radial];

	tau_rp_1 = nu * sigma * drp_1;
	tau_rp_2 = nu * sigma * (drp_2 + drp_r);
	} else {

	drp_2 = Ra[n_radial] * (- InvRb[n_radial - 1]) *
		   InvDiffRmed[n_radial];
	tau_rp_1 = nu * sigma * drp_2;
	tau_rp_2 = nu * sigma * (drp_1 + drp_r);
	}
}

///
/// \brief get_tau_pp
/// Split up the Tau_phi_phi in two components such that Tau_phi_phi = Tpp1*vphi + Tpp2
///
static void get_phi_pp(t_data &data, double &tau_pp_1, double &tau_pp_2, const int n_radial, const int n_azimuthal, bool include_artifical_viscosity, bool p_id_m)
{


	// div(v) = 1/r d(r*v_r)/dr + 1/r d(v_phi)/dphi

	const double DIV_V_org = data[t_data::DIV_V](n_radial, n_azimuthal) =
			(data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) *
				 Ra[n_radial + 1] -
			 data[t_data::V_RADIAL](n_radial, n_azimuthal) * Ra[n_radial]) *
				InvDiffRsup[n_radial] * InvRb[n_radial]
			+
			(data[t_data::V_AZIMUTHAL](
				 n_radial,
				 n_azimuthal == data[t_data::DIV_V].get_max_azimuthal()
				 ? 0
				 : n_azimuthal + 1) -
			 data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal)) *
				invdphi * InvRb[n_radial];

	const double DIV_V_R = (data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) *
			Ra[n_radial + 1] -
		data[t_data::V_RADIAL](n_radial, n_azimuthal) * Ra[n_radial]) *
		   InvDiffRsup[n_radial] * InvRb[n_radial];


	double DIV_V_P_1 = -data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) *
			   invdphi * InvRb[n_radial];


	double DIV_V_P_2 = data[t_data::V_AZIMUTHAL](
				n_radial,
				n_azimuthal == data[t_data::DIV_V].get_max_azimuthal()
				? 0
				: n_azimuthal + 1) *
			   invdphi * InvRb[n_radial];

	const double DIV_V =
	DIV_V_R + (DIV_V_P_2 + DIV_V_P_1);

	if(DIV_V != 0.0){
	if(std::fabs((DIV_V - DIV_V_org)/DIV_V) > 1e-4){
		printf("DIV_V (%d %d) failed with upd = %.5e	alt = %.5e	rel = %.5e\n", n_radial, n_azimuthal, DIV_V_org, DIV_V, std::fabs(DIV_V_org - DIV_V)/DIV_V_org);
	}}


	const double dpp_org = (data[t_data::V_AZIMUTHAL](
		   n_radial,
		   n_azimuthal ==
			   data[t_data::TAU_PHI_PHI].get_max_azimuthal()
		   ? 0
		   : n_azimuthal + 1) -
	   data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal)) *
		  invdphi * InvRmed[n_radial] +
	  0.5 *
		  (data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) +
		   data[t_data::V_RADIAL](n_radial, n_azimuthal)) *
		  InvRmed[n_radial];

	// 1/r d(v_phi)/dphi + v_r/r (cell centered)
	double dpp_1 = -data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) *
		  invdphi * InvRmed[n_radial];

	double dpp_2 = data[t_data::V_AZIMUTHAL](
		   n_radial,
		   n_azimuthal ==
			   data[t_data::TAU_PHI_PHI].get_max_azimuthal()
		   ? 0
		   : n_azimuthal + 1) *
		  invdphi * InvRmed[n_radial];

	const double dpp_r = 0.5 * (data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) +
		   data[t_data::V_RADIAL](n_radial, n_azimuthal)) *
		  InvRmed[n_radial];

	const double dpp = dpp_1 + dpp_2 + dpp_r;

	if(dpp != 0.0){
	if(std::fabs((dpp - dpp_org)/dpp) > 1e-9){
		printf("dpp failed with upd = %.5e	alt = %.5e	rel = %.5e\n", dpp_org, dpp, std::fabs(dpp_org - dpp)/dpp_org);
	}}

	// tau_phi_phi = 2*nu*Sigma*( 1/r d(v_phi)/dphi + v_r/r - 1/3 div(v)
	// )
	const double nu = data[t_data::VISCOSITY](n_radial, n_azimuthal);
	const double sigma = data[t_data::DENSITY](n_radial, n_azimuthal);

	double tpp = 2.0 * nu * sigma * (dpp_org - 1.0 / 3.0 * DIV_V_org);

	const double ttau_pp_1 = 2.0 * nu * sigma * (dpp_1 - 1.0 / 3.0 * DIV_V_P_1);
	const double ttau_pp_2 = 2.0 * nu * sigma * ((dpp_2 + dpp_r) - 1.0 / 3.0 * (DIV_V_P_2 + DIV_V_R));
	const double vp = data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal);
	const double vp_jm = data[t_data::V_AZIMUTHAL](
				n_radial,
				n_azimuthal ==
					data[t_data::TAU_PHI_PHI].get_max_azimuthal()
				? 0
				: n_azimuthal + 1);

	const double tppt = data[t_data::TAU_PHI_PHI](n_radial, n_azimuthal);


	if(tpp != 0.0){
	if(std::fabs((tpp - (ttau_pp_1 + ttau_pp_2))/tpp) > 1e-4){
		printf("tpp1 (%d	%d) failed with	%.5e	alt = %.5e	rel = %.5e\n", n_radial, n_azimuthal, tpp, (ttau_pp_1 + ttau_pp_2), std::fabs((tpp - (ttau_pp_1 + ttau_pp_2))/tpp));
	}}


	if(!p_id_m){

	dpp_1 = - invdphi * InvRmed[n_radial];
	DIV_V_P_1 = - invdphi * InvRb[n_radial];
	tau_pp_1 = 2.0 * nu * sigma * (dpp_1 - 1.0 / 3.0 * DIV_V_P_1);
	tau_pp_2 = 2.0 * nu * sigma * ((dpp_2 + dpp_r) - 1.0 / 3.0 * (DIV_V_P_2 + DIV_V_R));
	} else {

	DIV_V_P_2 = invdphi * InvRb[n_radial];
	dpp_2 = invdphi * InvRmed[n_radial];
	tau_pp_1 = 2.0 * nu * sigma * (dpp_2 - 1.0 / 3.0 * DIV_V_P_2);
	tau_pp_2 = 2.0 * nu * sigma * ((dpp_1 + dpp_r) - 1.0 / 3.0 * (DIV_V_P_1 + DIV_V_R));
	}


	if(include_artifical_viscosity){
	double nu_art;
	if (DIV_V < 0) {
		nu_art =
		parameters::artificial_viscosity_factor *
		data[t_data::DENSITY](n_radial, n_azimuthal) *
		std::pow(std::min(
			Rsup[n_radial] - Rinf[n_radial],
			Rmed[n_radial] * 2 * M_PI /
			data[t_data::DENSITY].get_size_azimuthal()), 2) *
		(-DIV_V);
	} else {
		nu_art = 0;
	}

	if(!p_id_m){
	tau_pp_1 += nu_art * DIV_V_P_1;
	tau_pp_2 += nu_art * (DIV_V_P_2 + DIV_V_R);

	const double tau_pp = vp*tau_pp_1 + tau_pp_2;
	if(tppt != 0.0){
	if(std::fabs((tppt - tau_pp)/tppt) > 1e-4){
		printf("tppt true (%d	%d) failed with	%.5e	alt = %.5e	rel = %.5e\n", n_radial, n_azimuthal, tppt, tau_pp, std::fabs((tppt - tau_pp)/tppt));
	}}

	} else {
	tau_pp_1 += nu_art * DIV_V_P_2;
	tau_pp_2 += nu_art * (DIV_V_P_1 + DIV_V_R);

	const double tau_pp = vp_jm*tau_pp_1 + tau_pp_2;
	if(tppt != 0.0){
	if(std::fabs((tppt - tau_pp)/tppt) > 1e-4){
		printf("tppt else (%d	%d) failed with	%.5e	alt = %.5e	rel = %.5e\n", n_radial, n_azimuthal, tppt, tau_pp, std::fabs((tppt - tau_pp)/tppt));
	}}
	}
	}
}

///
/// \brief get_tau_pp
/// Split up the Tau_phi_phi in two componentss such that Tau_phi_phi = Tpp1*vr + Tpp2
///
static void get_r_pp(t_data &data, double &tau_pp_1, double &tau_pp_2, const int n_radial, const int n_azimuthal, bool include_artifical_viscosity, bool r_id_m)
{
	// div(v) = 1/r d(r*v_r)/dr + 1/r d(v_phi)/dphi
	const double DIV_V_org = data[t_data::DIV_V](n_radial, n_azimuthal) =
			(data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) *
				 Ra[n_radial + 1] -
			 data[t_data::V_RADIAL](n_radial, n_azimuthal) * Ra[n_radial]) *
				InvDiffRsup[n_radial] * InvRb[n_radial] +
			(data[t_data::V_AZIMUTHAL](
				 n_radial,
				 n_azimuthal == data[t_data::DIV_V].get_max_azimuthal()
				 ? 0
				 : n_azimuthal + 1) -
			 data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal)) *
				invdphi * InvRb[n_radial];

	double DIV_V_R_1 = (- data[t_data::V_RADIAL](n_radial, n_azimuthal) * Ra[n_radial]) *
		   InvDiffRsup[n_radial] * InvRb[n_radial];

	double DIV_V_R_2 = (data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) *
			Ra[n_radial + 1]) * InvDiffRsup[n_radial] * InvRb[n_radial];

	const double DIV_V_P = (data[t_data::V_AZIMUTHAL](
				n_radial,
				n_azimuthal == data[t_data::DIV_V].get_max_azimuthal()
				? 0
				: n_azimuthal + 1) -
			data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal)) *
			   invdphi * InvRb[n_radial];

	const double DIV_V =
	(DIV_V_R_2 + DIV_V_R_1) + DIV_V_P;

	if(DIV_V != 0.0){
	if(std::fabs((DIV_V - DIV_V_org)/DIV_V) > 1e-4){
		printf("DIV_V (%d %d) failed with upd = %.5e	alt = %.5e	rel = %.5e\n", n_radial, n_azimuthal, DIV_V_org, DIV_V, std::fabs(DIV_V_org - DIV_V)/DIV_V_org);
	}}


	const double dpp_org = (data[t_data::V_AZIMUTHAL](
		   n_radial,
		   n_azimuthal ==
			   data[t_data::TAU_PHI_PHI].get_max_azimuthal()
		   ? 0
		   : n_azimuthal + 1) -
	   data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal)) *
		  invdphi * InvRmed[n_radial] +
	  0.5 *
		  (data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) +
		   data[t_data::V_RADIAL](n_radial, n_azimuthal)) *
		  InvRmed[n_radial];

	// 1/r d(v_phi)/dphi + v_r/r (cell centered)
	const double dpp_p = (data[t_data::V_AZIMUTHAL](
				n_radial,
				n_azimuthal ==
					data[t_data::TAU_PHI_PHI].get_max_azimuthal()
				? 0
				: n_azimuthal + 1) -
			data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal)) *
			   invdphi * InvRmed[n_radial];

	double dpp_1 = 0.5 * data[t_data::V_RADIAL](n_radial, n_azimuthal) *
		  InvRmed[n_radial];

	double dpp_2 = 0.5 * data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) *
		  InvRmed[n_radial];

	const double dpp = dpp_p + dpp_1 + dpp_2;

	if(dpp != 0.0){
	if(std::fabs((dpp - dpp_org)/dpp) > 1e-9){
		printf("dpp failed with upd = %.5e	alt = %.5e	rel = %.5e\n", dpp_org, dpp, std::fabs(dpp_org - dpp)/dpp_org);
	}}

	//tau_phi_phi = 2*nu*Sigma*( 1/r d(v_phi)/dphi + v_r/r - 1/3 div(v))
	const double nu = data[t_data::VISCOSITY](n_radial, n_azimuthal);
	const double sigma = data[t_data::DENSITY](n_radial, n_azimuthal);

	if(!r_id_m){

	DIV_V_R_1 = - Ra[n_radial] * InvDiffRsup[n_radial] * InvRb[n_radial];
	dpp_1 = 0.5 * InvRmed[n_radial];

	/*if(n_radial == 62 && n_azimuthal == 5){
		//printf("calc pp1 = %.5e	%.5e	%.5e", -nu * sigma, dpp_1 - 1.0 / 3.0 * DIV_V_R_1, DIV_V_R_1);
		printf("calc pp1\n");
		printf("NuSigma1 = %.12e\n", -nu * sigma);
		printf("dpp1 = %.12e\n", dpp_1);
		printf("DIV_V1 = %.12e\n", - 1.0 / 3.0 * DIV_V_R_1);
		printf("Nu1 = %.12e\n", -2.0 * nu * sigma * (dpp_1 - 1.0 / 3.0 * DIV_V_R_1));
	}*/

	tau_pp_1 = 2.0 * nu * sigma * (dpp_1 - 1.0 / 3.0 * DIV_V_R_1);
	tau_pp_2 = 2.0 * nu * sigma * ((dpp_p + dpp_2) - 1.0 / 3.0 * (DIV_V_R_2 + DIV_V_P));
	} else {

	DIV_V_R_2 = Ra[n_radial + 1] * InvDiffRsup[n_radial] * InvRb[n_radial];
	dpp_2 = 0.5 * InvRmed[n_radial];

	/*if(n_radial == 61 && n_azimuthal == 5){
		//printf("calc pp2 = %.5e	%.5e	%.5e",- nu * sigma, dpp_2 - 1.0 / 3.0 * DIV_V_R_2, DIV_V_R_2);
		printf("calc pp2\n");
		printf("NuSigma2 = %.12e\n", -nu * sigma);
		printf("dpp2 = %.12e\n", dpp_2);
		printf("DIV_V2 = %.12e\n", - 1.0 / 3.0 * DIV_V_R_2);
		printf("Nu2 = %.12e\n", -2.0 * nu * sigma * (dpp_2 - 1.0 / 3.0 * DIV_V_R_2));
	}*/

	tau_pp_1 = 2.0 * nu * sigma * (dpp_2 - 1.0 / 3.0 * DIV_V_R_2);
	tau_pp_2 = 2.0 * nu * sigma * ((dpp_p + dpp_1) - 1.0 / 3.0 * (DIV_V_R_1 + DIV_V_P));
	}


	if(include_artifical_viscosity){
	double nu_art;
	if (DIV_V < 0) {
		nu_art =
		parameters::artificial_viscosity_factor *
		data[t_data::DENSITY](n_radial, n_azimuthal) *
		std::pow(std::min(
			Rsup[n_radial] - Rinf[n_radial],
			Rmed[n_radial] * 2 * M_PI /
			data[t_data::DENSITY].get_size_azimuthal()), 2) *
		(-DIV_V);
	} else {
		nu_art = 0;
	}

	if(!r_id_m){
	tau_pp_1 += nu_art * DIV_V_R_1;
	//if(n_radial == 62 && n_azimuthal == 5)
	//	printf("NuArt1 = %.12e	%.5e	%.5e\n\n", -0.5*nu_art * DIV_V_R_1, nu_art, DIV_V_R_1);

	tau_pp_2 += nu_art * (DIV_V_R_2 + DIV_V_P);
	} else {
	tau_pp_1 += nu_art * DIV_V_R_2;
	//if(n_radial == 61 && n_azimuthal == 5)
	//	printf("NuArt2 = %.12e\n\n", -0.5*nu_art * DIV_V_R_2);

	tau_pp_2 += nu_art * (DIV_V_R_1 + DIV_V_P);
	}
	}

}

///
/// \brief get_tau_r_rp
/// Split up the Tau_r_phi in two componentss such that Tau_r_phi = Trp1*vr + Trp2
///
static void get_tau_r_rp(t_data &data, double &tau_rp_1, double &tau_rp_2, const int n_radial, const int n_azimuthal, bool p_id_p)
{

	const int n_azimuthal_minus =
	(n_azimuthal == 0 ? data[t_data::DENSITY].get_max_azimuthal()
			  : n_azimuthal - 1);

	double dvphirdr =
	(data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) *
		 InvRb[n_radial] -
	 data[t_data::V_AZIMUTHAL](n_radial - 1, n_azimuthal) *
		 InvRb[n_radial - 1]) *
	InvDiffRmed[n_radial];
	// d(v_r)/dphi
	double dvrdphi =
	(data[t_data::V_RADIAL](n_radial, n_azimuthal) -
	 data[t_data::V_RADIAL](
		 n_radial, n_azimuthal == 0
			   ? data[t_data::V_RADIAL].get_max_azimuthal()
			   : n_azimuthal - 1)) *
	invdphi;

	// r*d(v_phi/r)/dr + 1/r d(v_r)/dphi (edge)
	double drp_org = Ra[n_radial] * dvphirdr + dvrdphi * InvRa[n_radial];

	double drp_p = Ra[n_radial] * dvphirdr;

	double drp_1 = data[t_data::V_RADIAL](n_radial, n_azimuthal) *
		   invdphi * InvRa[n_radial];

	double drp_2 = -data[t_data::V_RADIAL](
				n_radial, n_azimuthal == 0
					  ? data[t_data::V_RADIAL].get_max_azimuthal()
					  : n_azimuthal - 1) *
		   invdphi * InvRa[n_radial];

	double drp = drp_p + (drp_1 + drp_2);

	if(drp != 0.0){
	if(std::fabs((drp - drp_org)/drp) > 1e-10){
		printf("drp (%d %d) failed with upd = %.5e	alt = %.5e	rel = %.5e\n", n_radial, n_azimuthal, drp_org, drp, std::fabs(drp_org - drp)/drp_org);
	}}


	// averaged nu over 4 corresponding cells
	const double nu =
	0.25 *
	(data[t_data::VISCOSITY](n_radial, n_azimuthal) +
	 data[t_data::VISCOSITY](n_radial - 1, n_azimuthal) +
	 data[t_data::VISCOSITY](n_radial, n_azimuthal_minus) +
	 data[t_data::VISCOSITY](n_radial - 1, n_azimuthal_minus));

	// averaged sigma over 4 corresponding cells
	const double sigma =
	0.25 * (data[t_data::DENSITY](n_radial, n_azimuthal) +
		data[t_data::DENSITY](n_radial - 1, n_azimuthal) +
		data[t_data::DENSITY](n_radial, n_azimuthal_minus) +
		data[t_data::DENSITY](n_radial - 1, n_azimuthal_minus));


	const double vr = data[t_data::V_RADIAL](n_radial, n_azimuthal);
	const double vr_im = data[t_data::V_RADIAL](
				n_radial, n_azimuthal == 0
					  ? data[t_data::V_RADIAL].get_max_azimuthal()
					  : n_azimuthal - 1);
	const double trp = data[t_data::TAU_R_PHI](n_radial, n_azimuthal);
	const double ttau_rp_1 = nu * sigma * drp_1;
	const double ttau_rp_2 = nu * sigma * (drp_2 + drp_p);
	if(trp != 0.0){
	if(std::fabs((trp - (ttau_rp_1 + ttau_rp_2))/trp) > 2e-10){
		printf("trp1 failed with	%.5e	alt = %.5e	rel = %.5e\n", trp, (ttau_rp_1 + ttau_rp_2), std::fabs((trp - (ttau_rp_1 + ttau_rp_2))/trp));
	}}


	const double ettau_rp_1 = nu * sigma * drp_2;
	const double ettau_rp_2 = nu * sigma * (drp_1 + drp_p);
	if(trp != 0.0){
	if(std::fabs((trp - (ettau_rp_1 + ettau_rp_2))/trp) > 1e-10){
		printf("etrp1 failed with	%.5e	alt = %.5e	rel = %.5e\n", trp, (ettau_rp_1 + ettau_rp_2), std::fabs((trp - (ettau_rp_1 + ettau_rp_2))/trp));
	}}

	if(!p_id_p){
	drp_1 = invdphi * InvRa[n_radial];
	tau_rp_1 = nu * sigma * drp_1;
	tau_rp_2 = nu * sigma * (drp_2 + drp_p);

	if(trp != 0.0){
	if(std::fabs((trp - (vr*tau_rp_1 + tau_rp_2))/trp) > 2e-10){
	printf("\ntrue\n");
	printf("trp2 (%d %d) %.3e failed with	%.5e	alt = %.5e	rel = %.5e\n", n_radial, n_azimuthal, vr, trp, (vr*tau_rp_1 + tau_rp_2), std::fabs((trp - (vr*tau_rp_1 + tau_rp_2))/trp));
	printf("trp test	(%.5e	%.5e	%.5e)\n", data[t_data::TAU_R_PHI](n_radial, n_azimuthal), ttau_rp_1+ttau_rp_2, vr*tau_rp_1 + tau_rp_2);
	printf("trp12 test	(%.5e	%.5e)	(%.5e	%.5e)\n", ttau_rp_1, vr*tau_rp_1, ttau_rp_2, tau_rp_2);
	}}

	} else {
	drp_2 = - invdphi * InvRa[n_radial];
	tau_rp_1 = nu * sigma * drp_2;
	tau_rp_2 = nu * sigma * (drp_1 + drp_p);

	if(trp != 0.0){
	if(std::fabs((trp - (vr_im*tau_rp_1 + tau_rp_2))/trp) > 1e-10){
		printf("\nelse\n");
		printf("trp2 (%d %d) %.3e failed with	%.5e	alt = %.5e	rel = %.5e\n", n_radial, n_azimuthal, vr_im, trp, (vr_im*tau_rp_1 + tau_rp_2), std::fabs((trp - (vr_im*tau_rp_1 + tau_rp_2))/trp));
		printf("trp test	(%.5e	%.5e	%.5e)\n", data[t_data::TAU_R_PHI](n_radial, n_azimuthal), ettau_rp_1+ettau_rp_2, vr_im*tau_rp_1 + tau_rp_2);
		printf("trp12 test	(%.5e	%.5e)	(%.5e	%.5e)\n", ettau_rp_1, vr_im*tau_rp_1, ettau_rp_2, tau_rp_2);
	}}
	}
}

///
/// \brief get_tau_rr
/// Split up the Tau_r_r in two componentss such that Tau_r_r = Trr1*vr + Trr2
///
static void get_tau_rr(t_data &data, double &tau_rr_1, double &tau_rr_2, const int n_radial, const int n_azimuthal, bool include_artifical_viscosity, bool r_id_m)
{

	// div(v) = 1/r d(r*v_r)/dr + 1/r d(v_phi)/dphi
	const double DIV_V_org =
			(data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) *
				 Ra[n_radial + 1] -
			 data[t_data::V_RADIAL](n_radial, n_azimuthal) * Ra[n_radial]) *
				InvDiffRsup[n_radial] * InvRb[n_radial] +
			(data[t_data::V_AZIMUTHAL](
				 n_radial,
				 n_azimuthal == data[t_data::DIV_V].get_max_azimuthal()
				 ? 0
				 : n_azimuthal + 1) -
			 data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal)) *
				invdphi * InvRb[n_radial];

	double DIV_V_R_1 = (- data[t_data::V_RADIAL](n_radial, n_azimuthal) * Ra[n_radial]) *
		   InvDiffRsup[n_radial] * InvRb[n_radial];

	double DIV_V_R_2 = (data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) *
			Ra[n_radial + 1]) * InvDiffRsup[n_radial] * InvRb[n_radial];

	const double DIV_V_P = (data[t_data::V_AZIMUTHAL](
				n_radial,
				n_azimuthal == data[t_data::DIV_V].get_max_azimuthal()
				? 0
				: n_azimuthal + 1) -
			data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal)) *
			   invdphi * InvRb[n_radial];

	const double DIV_V =
	(DIV_V_R_2 + DIV_V_R_1) + DIV_V_P;

	if(DIV_V != 0.0){
	if(std::fabs((DIV_V - DIV_V_org)/DIV_V) > 1e-5){
		printf("DIV_V (%d %d) failed with upd = %.5e	alt = %.5e	rel = %.5e\n", n_radial, n_azimuthal, DIV_V_org, DIV_V, std::fabs(DIV_V_org - DIV_V)/DIV_V_org);
	}}

	// calculate tau_r_r
	// d(v_r)/dr (cell centered)
	double drr_1 = -data[t_data::V_RADIAL](n_radial, n_azimuthal) *
			InvDiffRsup[n_radial];
	double drr_2 = data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) *
			InvDiffRsup[n_radial];

	// tau_r_r = 2*nu*Sigma*( d(v_r)/dr - 1/3 div(v) + eta2 div(v))
	if(!r_id_m){

		drr_1 = - InvDiffRsup[n_radial];
		DIV_V_R_1 = - Ra[n_radial] * InvDiffRsup[n_radial] * InvRb[n_radial];
	/*if(n_radial == 50 && n_azimuthal == 0){
		printf("\ncalc crr1\n");
		printf("NuSigma1 = %.12e\n", 2.0 * data[t_data::VISCOSITY](n_radial, n_azimuthal) *
				data[t_data::DENSITY](n_radial, n_azimuthal));
		printf("drr1 = %.12e\n", drr_1);
		printf("DIV_V1 = %.12e\n", - 1.0 / 3.0 * DIV_V_R_1);
		printf("Nu1 = %.12e\n", 2.0 * data[t_data::VISCOSITY](n_radial, n_azimuthal) *
				data[t_data::DENSITY](n_radial, n_azimuthal) *
				(drr_1 - 1.0 / 3.0 * DIV_V_R_1));
	}*/

	tau_rr_1 = 2.0 * data[t_data::VISCOSITY](n_radial, n_azimuthal) *
			data[t_data::DENSITY](n_radial, n_azimuthal) *
			(drr_1 - 1.0 / 3.0 * DIV_V_R_1);
	tau_rr_2 = 2.0 * data[t_data::VISCOSITY](n_radial, n_azimuthal) *
			data[t_data::DENSITY](n_radial, n_azimuthal) *
			(drr_2 - 1.0 / 3.0 * (DIV_V_R_2 + DIV_V_P));
	} else {

		drr_2 = InvDiffRsup[n_radial];
		DIV_V_R_2 =	Ra[n_radial + 1] * InvDiffRsup[n_radial] * InvRb[n_radial];

		/*if(n_radial == 49 && n_azimuthal == 0){
			printf("\ncalc crr2\n");
			printf("NuSigma2 = %.12e\n", 2.0 * data[t_data::VISCOSITY](n_radial, n_azimuthal) *
					data[t_data::DENSITY](n_radial, n_azimuthal));
			printf("drr2 = %.12e\n", drr_2);
			printf("DIV_V2 = %.12e\n",  - 1.0 / 3.0 * DIV_V_R_2);
			printf("Nu2 = %.12e\n", 2.0 * data[t_data::VISCOSITY](n_radial, n_azimuthal) *
					data[t_data::DENSITY](n_radial, n_azimuthal) *
					(drr_2 - 1.0 / 3.0 * DIV_V_R_2));
		}*/

	tau_rr_1 = 2.0 * data[t_data::VISCOSITY](n_radial, n_azimuthal) *
			data[t_data::DENSITY](n_radial, n_azimuthal) *
			(drr_2 - 1.0 / 3.0 * DIV_V_R_2);
	tau_rr_2 = 2.0 * data[t_data::VISCOSITY](n_radial, n_azimuthal) *
			data[t_data::DENSITY](n_radial, n_azimuthal) *
			(drr_1 - 1.0 / 3.0 * (DIV_V_R_1 + DIV_V_P));
	}

	if(include_artifical_viscosity){
	double nu_art;
	if (DIV_V < 0) {
		nu_art =
		parameters::artificial_viscosity_factor *
		data[t_data::DENSITY](n_radial, n_azimuthal) *
		std::pow(std::min(
			Rsup[n_radial] - Rinf[n_radial],
			Rmed[n_radial] * 2 * M_PI /
			data[t_data::DENSITY].get_size_azimuthal()), 2) *
		(-DIV_V);
	} else {
		nu_art = 0;
	}
	if(!r_id_m){

	/*if(n_radial == 50 && n_azimuthal == 0){
		printf("NuArt1 = %.12e\n\n", nu_art * DIV_V_R_1);
	}*/

	tau_rr_1 += nu_art * DIV_V_R_1;
	tau_rr_2 += nu_art * (DIV_V_R_2 + DIV_V_P);
	} else {

		/*if(n_radial == 49 && n_azimuthal == 0){
			printf("NuArt2 = %.12e\n\n", nu_art * DIV_V_R_2);
		}*/

	tau_rr_1 += nu_art * DIV_V_R_2;
	tau_rr_2 += nu_art * (DIV_V_R_1 + DIV_V_P);
	}
	}
}

///
/// \brief debug_function_viscous_terms
///	took the standart viscosity velocity updates, split them into their components
/// then sort them by the velocity of the cell such that delta v = c1*v + c2
/// compared the resulting velocity update with the original computation
/// and compared the resulting constants with the one computed in compute_viscous_terms
///
void debug_function_viscous_terms(t_data &data, bool include_artifical_viscosity, const double dt)
{

	int n_azimuthal_plus, n_azimuthal_minus;

	for (unsigned int n_radial = 1; n_radial <= data[t_data::V_RADIAL].get_max_radial() - 1;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::V_RADIAL].get_max_azimuthal(); ++n_azimuthal) {
		n_azimuthal_plus =
		(n_azimuthal == data[t_data::DENSITY].get_max_azimuthal()
			 ? 0
			 : n_azimuthal + 1);
		n_azimuthal_minus =
		(n_azimuthal == 0 ? data[t_data::DENSITY].get_max_azimuthal()
				  : n_azimuthal - 1);

		// a_phi = 1/(r*Sigma) ( d(r*tau_r_phi)/dr + d(tau_phi_phi)/dphi +
		const double sigma_avg_phi =
		0.5 * (data[t_data::DENSITY](n_radial, n_azimuthal) +
			   data[t_data::DENSITY](n_radial, n_azimuthal_minus));
		const double v_phi_upd = dt * InvRb[n_radial] / (sigma_avg_phi) *
				(TwoDiffRaSq[n_radial] * (std::pow(Ra[n_radial + 1], 2) *
					  data[t_data::TAU_R_PHI](n_radial + 1,
		n_azimuthal) - std::pow(Ra[n_radial], 2) *
					  data[t_data::TAU_R_PHI](n_radial,
		n_azimuthal)) + (data[t_data::TAU_PHI_PHI](n_radial,
		n_azimuthal) - data[t_data::TAU_PHI_PHI](n_radial,
		n_azimuthal_minus)) * invdphi);

		const double v_phi_upd_rp = dt * InvRb[n_radial] / (sigma_avg_phi) *
				(TwoDiffRaSq[n_radial] * (std::pow(Ra[n_radial + 1], 2) *
					  data[t_data::TAU_R_PHI](n_radial + 1,
		n_azimuthal) - std::pow(Ra[n_radial], 2) *
					  data[t_data::TAU_R_PHI](n_radial,
		n_azimuthal)));

		const double v_phi_upd_pp = dt * InvRb[n_radial] / (sigma_avg_phi) *
				((data[t_data::TAU_PHI_PHI](n_radial,
		n_azimuthal) - data[t_data::TAU_PHI_PHI](n_radial,
		n_azimuthal_minus)) * invdphi);

		//// PHI UPD

		const double vp = data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal);
		double TAU_RP_1;
		double TAU_RP_2;
		double TAU_RP_ip_1;
		double TAU_RP_ip_2;
		get_tau_rp(data, TAU_RP_1, TAU_RP_2, n_radial, n_azimuthal, false);
		get_tau_rp(data, TAU_RP_ip_1, TAU_RP_ip_2, n_radial+1, n_azimuthal, true);

		double crp_org = dt * InvRb[n_radial] / (sigma_avg_phi) * (TwoDiffRaSq[n_radial] * (std::pow(Ra[n_radial + 1], 2) * (vp*TAU_RP_ip_1 + TAU_RP_ip_2) - std::pow(Ra[n_radial], 2) * (vp*TAU_RP_1 + TAU_RP_2)));

		double crp_1 = InvRb[n_radial] / (sigma_avg_phi) * (TwoDiffRaSq[n_radial] * (std::pow(Ra[n_radial + 1], 2) * TAU_RP_ip_1 - std::pow(Ra[n_radial], 2) * TAU_RP_1));
		double crp_2 = InvRb[n_radial] / (sigma_avg_phi) * (TwoDiffRaSq[n_radial] * (std::pow(Ra[n_radial + 1], 2) * TAU_RP_ip_2 - std::pow(Ra[n_radial], 2) * TAU_RP_2));

		const double crp = dt*(vp*crp_1 + crp_2);

		if(crp != 0.0){
		if(std::fabs((crp - crp_org)/crp) > 1e-4 && crp / (std::fabs(vp*crp_1) + std::fabs(crp_2)) > 5e-14){
			printf("crp failed with upd = %.5e	alt = %.5e	rel = %.5e	c12 = (%.3e	%.3e)\n", crp_org, crp, std::fabs(crp_org - crp)/crp_org, vp*crp_1, crp_2);
		}}

		double TAU_PP_1;
		double TAU_PP_2;
		double TAU_PP_jm_1;
		double TAU_PP_jm_2;
		get_phi_pp(data, TAU_PP_1, TAU_PP_2, n_radial, n_azimuthal, include_artifical_viscosity, false);
		get_phi_pp(data, TAU_PP_jm_1, TAU_PP_jm_2, n_radial, n_azimuthal_minus, include_artifical_viscosity, true);
		double cpp_org = dt * InvRb[n_radial] / (sigma_avg_phi) * ((vp*TAU_PP_1 + TAU_PP_2) - (vp*TAU_PP_jm_1 + TAU_PP_jm_2)) * invdphi;

		double cpp_1 = InvRb[n_radial] / (sigma_avg_phi) * (TAU_PP_1 - TAU_PP_jm_1) * invdphi;
		double cpp_2 = InvRb[n_radial] / (sigma_avg_phi) * (TAU_PP_2 - TAU_PP_jm_2) * invdphi;

		double cpp = dt*(vp*cpp_1 + cpp_2);

		if(cpp != 0.0){
		if(std::fabs((cpp - cpp_org)/cpp) > 1e-4 && cpp / (std::fabs(vp*cpp_1) + std::fabs(cpp_2)) > 1e-13){
			printf("cpp failed with upd = %.5e	alt = %.5e	rel = %.5e\n", cpp_org, cpp, std::fabs(cpp_org - cpp)/cpp_org);
		}}

		if(v_phi_upd != 0.0){
		if(std::fabs((v_phi_upd - v_phi_upd_rp - v_phi_upd_pp)/v_phi_upd) > 1e-9){
			printf("v_phi consts (%d %d) failed with upd = %.5e	alt = %.5e	rel = %.5e\n", n_radial, n_azimuthal, v_phi_upd, v_phi_upd_rp+v_phi_upd_pp, std::fabs(v_phi_upd - v_phi_upd_rp - v_phi_upd_pp)/v_phi_upd);
		}}

		if(v_phi_upd_rp != 0.0){
		if(std::fabs((v_phi_upd_rp - crp)/v_phi_upd_rp) > 1e-5 && crp / (std::fabs(vp*crp_1) + std::fabs(crp_2)) > 5e-13){
			printf("v_phi_rp (%d %d) failed with upd = %.5e	alt = %.5e	rel = %.5e\n", n_radial, n_azimuthal, v_phi_upd_rp, crp, std::fabs(v_phi_upd_rp - crp)/v_phi_upd_rp);
		}}

		if(v_phi_upd_pp != 0.0){
		if(std::fabs((v_phi_upd_pp - cpp)/v_phi_upd_pp) > 1e-2  && cpp / (std::fabs(vp*cpp_1) + std::fabs(cpp_2)) > 1e-14){
			printf("v_phi_pp (%d %d) failed with upd = %.5e	alt = %.5e	rel = %.5e\n", n_radial, n_azimuthal, v_phi_upd_pp, cpp, std::fabs(v_phi_upd_pp - cpp)/v_phi_upd_pp);
		}}

		const double add_crp_cpp = (crp + cpp);
		if(v_phi_upd != 0.0){
		if(std::fabs((v_phi_upd - add_crp_cpp)/v_phi_upd) > 1e-3  && std::fabs(add_crp_cpp/vp) > 1e-13){
			printf("v_phi (%d %d) failed with upd = %.5e	alt = %.5e	vp = %.5e	rel = %.5e\n", n_radial, n_azimuthal, v_phi_upd, cpp+crp, vp, std::fabs(v_phi_upd - add_crp_cpp)/v_phi_upd);
		}}

		/// test passed
		//if(n_radial == 50 && n_azimuthal == 0)
		//	printf("calc cpp = %.12e	%.12e\n", cpp_1, crp_1);


		const double c1_alt = (cpp_1 + crp_1);
		const double c1 = data[t_data::VISCOSITY_CORRECTION_FACTOR_PHI](n_radial, n_azimuthal);
		if(c1 != 0.0){
		if(std::fabs((c1 - c1_alt)/c1) > 1e-9){
			printf("c1 (%d %d) failed with upd = %.5e	alt = %.5e	rel = %.5e\n", n_radial, n_azimuthal, c1, c1_alt, std::fabs(c1 - c1_alt)/c1);
		}}


		//if(n_radial == 50 && n_azimuthal == 0)
		//	printf("cpp = %.5e	%.5e\n", c1, c1_alt);

		//// END PHI UPD

		/// R UPD
		const double vr = data[t_data::V_RADIAL](n_radial, n_azimuthal);
		// a_r = 1/(r*Sigma) ( d(r*tau_r_r)/dr + d(tau_r_phi)/dphi -
		// tau_phi_phi )
		const double sigma_avg_r =
		0.5 * (data[t_data::DENSITY](n_radial, n_azimuthal) +
			   data[t_data::DENSITY](n_radial - 1, n_azimuthal));

		const double v_r_upd_org =
		dt * InvRinf[n_radial] /
		(sigma_avg_r)*parameters::radial_viscosity_factor *
		((Rmed[n_radial] * data[t_data::TAU_R_R](n_radial, n_azimuthal) -
		  Rmed[n_radial - 1] *
			  data[t_data::TAU_R_R](n_radial - 1, n_azimuthal)) *
			 InvDiffRmed[n_radial] +
		 (data[t_data::TAU_R_PHI](n_radial, n_azimuthal_plus) -
		  data[t_data::TAU_R_PHI](n_radial, n_azimuthal)) *
			 invdphi -
		 0.5 * (data[t_data::TAU_PHI_PHI](n_radial, n_azimuthal) +
			data[t_data::TAU_PHI_PHI](n_radial - 1, n_azimuthal)));

		const double v_r_upd_r_org =
		dt * InvRinf[n_radial] /
		(sigma_avg_r)*parameters::radial_viscosity_factor *
		((Rmed[n_radial] * data[t_data::TAU_R_R](n_radial, n_azimuthal) -
		  Rmed[n_radial - 1] *
			  data[t_data::TAU_R_R](n_radial - 1, n_azimuthal)) *
			 InvDiffRmed[n_radial]);

		double TAU_RR_1;
		double TAU_RR_2;
		double TAU_RR_im_1;
		double TAU_RR_im_2;
		get_tau_rr(data, TAU_RR_1, TAU_RR_2, n_radial, n_azimuthal, include_artifical_viscosity, false);
		get_tau_rr(data, TAU_RR_im_1, TAU_RR_im_2, n_radial - 1, n_azimuthal, include_artifical_viscosity, true);

		/// test passed
		/*
		if(n_radial == 50 && n_azimuthal == 0)
		printf("calc cr_rr = %.12e	%.12e	%.3e\n\n", Rmed[n_radial]*TAU_RR_1, -Rmed[n_radial - 1]*TAU_RR_im_1,
				(Rmed[n_radial] * TAU_RR_1 - Rmed[n_radial - 1] * TAU_RR_im_1) *
				   InvDiffRmed[n_radial]);*/

		const double crr_r1 =
		InvRinf[n_radial] /
		(sigma_avg_r)*parameters::radial_viscosity_factor *
		(Rmed[n_radial] * TAU_RR_1 - Rmed[n_radial - 1] * TAU_RR_im_1) *
			 InvDiffRmed[n_radial];

		const double crr_r2 =
		InvRinf[n_radial] /
		(sigma_avg_r)*parameters::radial_viscosity_factor *
		(Rmed[n_radial] * TAU_RR_2 - Rmed[n_radial - 1] * TAU_RR_im_2) *
			 InvDiffRmed[n_radial];

		/// test passed
		//if(n_radial == 50 && n_azimuthal == 0)
		//printf("calc cr_rr = %.12e	%.12e	%.12e\n", Rmed[n_radial] * TAU_RR_1 * InvDiffRmed[n_radial], - Rmed[n_radial - 1] * TAU_RR_im_1 * InvDiffRmed[n_radial], (Rmed[n_radial] * TAU_RR_1 - Rmed[n_radial - 1] * TAU_RR_im_1) *
		//		InvDiffRmed[n_radial]);

		const double crr = dt*(vr*crr_r1 + crr_r2);

		if(crr != 0.0){
		if(std::fabs((crr - v_r_upd_r_org)/crr) > 1e-7){
			printf("v_r_upd_r failed with upd = %.5e	alt = %.5e	rel = %.5e\n", v_r_upd_r_org, crr, std::fabs(crr - v_r_upd_r_org)/crr);
		}}

		const double v_r_upd_p_org =
		dt * InvRinf[n_radial] /
		(sigma_avg_r)*parameters::radial_viscosity_factor *
		(-0.5 * (data[t_data::TAU_PHI_PHI](n_radial, n_azimuthal) +
			data[t_data::TAU_PHI_PHI](n_radial - 1, n_azimuthal)));

		double TAU_PP_im_1;
		double TAU_PP_im_2;
		get_r_pp(data, TAU_PP_1, TAU_PP_2, n_radial, n_azimuthal, include_artifical_viscosity, false);
		get_r_pp(data, TAU_PP_im_1, TAU_PP_im_2, n_radial-1, n_azimuthal, include_artifical_viscosity, true);

		const double cpp_r1 =
				InvRinf[n_radial] /
				(sigma_avg_r)*parameters::radial_viscosity_factor *
				(-0.5 * (TAU_PP_1 + TAU_PP_im_1));

		/// test passed
		//if(n_radial == 62 && n_azimuthal == 5)
		//printf("calc cr_pp = %.12e	%.12e\n\n", -0.5 * TAU_PP_1,-0.5 * TAU_PP_im_1);

		const double cpp_r2 =
				InvRinf[n_radial] /
				(sigma_avg_r)*parameters::radial_viscosity_factor *
				(-0.5 * (TAU_PP_2 + TAU_PP_im_2));

		const double cppr = dt*(vr*cpp_r1 + cpp_r2);


		if(cppr != 0.0){
		if(std::fabs((cppr - v_r_upd_p_org)/cppr) > 1e-7){
			printf("v_r_upd_pp failed with upd = %.5e	alt = %.5e	rel = %.5e\n", v_r_upd_p_org, cppr, std::fabs(cppr - v_r_upd_p_org)/cppr);
		}}

		const double v_r_upd_rp_org =
		dt * InvRinf[n_radial] /
		(sigma_avg_r)*parameters::radial_viscosity_factor *
		((data[t_data::TAU_R_PHI](n_radial, n_azimuthal_plus) -
		  data[t_data::TAU_R_PHI](n_radial, n_azimuthal)) *
			 invdphi);

		double TAU_RP_jp_1;
		double TAU_RP_jp_2;
		get_tau_r_rp(data, TAU_RP_1, TAU_RP_2, n_radial, n_azimuthal, false);
		get_tau_r_rp(data, TAU_RP_jp_1, TAU_RP_jp_2, n_radial, n_azimuthal_plus, true);

		const double crpr_1 =
		InvRinf[n_radial] /
		(sigma_avg_r)*parameters::radial_viscosity_factor *
		(TAU_RP_jp_1 - TAU_RP_1) * invdphi;

		/// test passed
		//if(n_radial == 50 && n_azimuthal == 0)
		//printf("calc cr_rp = %.12e	%.12e	%.5e\n\n", (TAU_RP_jp_1) * invdphi, (-TAU_RP_1) * invdphi, (TAU_RP_jp_1 - TAU_RP_1) * invdphi);

		const double v_r_upd_rp_2 =
				InvRinf[n_radial] /
				(sigma_avg_r)*parameters::radial_viscosity_factor *
				(TAU_RP_jp_2 - TAU_RP_2) * invdphi;

		const double crpr = dt*(vr * crpr_1 + v_r_upd_rp_2);

		if(crpr != 0.0){
		if(std::fabs((crpr - v_r_upd_rp_org)/crpr) > 1e-4){
			printf("v_r_upd_rp failed with upd = %.5e	alt = %.5e	rel = %.5e	compare (rr=%.5e	pp=%.5e)\n", v_r_upd_rp_org, crpr, std::fabs(crpr - v_r_upd_rp_org)/crpr, v_r_upd_r_org, v_r_upd_p_org);
		}}

		if(v_phi_upd != 0.0){
		if(std::fabs((v_r_upd_org - crr - cppr - crpr)/v_r_upd_org) > 1e-5){
			printf("v_r consts (%d %d) failed with upd = %.5e	alt = %.5e	rel = %.5e\n", n_radial, n_azimuthal, v_r_upd_org, crr+cppr+crpr, std::fabs(v_r_upd_org - crr - cppr - crpr)/v_r_upd_org);
		}}

		/// test passed
		//if(n_radial == 50 && n_azimuthal == 0)
		//	printf("calc crr = %.12e	%.12e	%.12e		%.12e\n", crr_r1, crpr_1, cpp_r1, (data[t_data::VISCOSITY_CORRECTION_FACTOR_R](n_radial, n_azimuthal) - (crr_r1 + crpr_1 + cpp_r1))/data[t_data::VISCOSITY_CORRECTION_FACTOR_R](n_radial, n_azimuthal));

		//if(n_radial == 62 && n_azimuthal == 5)
		//	printf("calc crr = %.12e	%.12e	%.12e		%.12e\n", crr_r1, crpr_1, cpp_r1, (data[t_data::VISCOSITY_CORRECTION_FACTOR_R](n_radial, n_azimuthal) - (crr_r1 + crpr_1 + cpp_r1))/data[t_data::VISCOSITY_CORRECTION_FACTOR_R](n_radial, n_azimuthal));

		const double c1r_alt = cpp_r1 + crpr_1 + crr_r1;
		const double c1r = data[t_data::VISCOSITY_CORRECTION_FACTOR_R](n_radial, n_azimuthal);
		if(c1r != 0.0){
		if(std::fabs((c1r - c1r_alt)/c1r) > 1e-12){
			printf("c1r (%d %d) failed with upd = %.5e	alt = %.5e	rel = %.5e\n", n_radial, n_azimuthal, c1r, c1r_alt, std::fabs(c1r - c1r_alt)/c1r);
		}}


		//if(n_radial == 50 && n_azimuthal == 0)
		//	printf("crr = %.5e	%.5e\n", c1r, c1r_alt);


		/// END R UPD


	}
	}
}

} // namespace viscosity
