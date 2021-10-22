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
	for (unsigned int n_radial = 0;
	     n_radial <= data[t_data::VISCOSITY].get_max_radial(); ++n_radial) {
	    const double inv_omega_kepler =
		1.0 / calculate_omega_kepler(Rb[n_radial]);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::VISCOSITY].get_max_azimuthal();
		 ++n_azimuthal) {
		// ν = α * c_s^2 / Ω_K
		data[t_data::VISCOSITY](n_radial, n_azimuthal) =
		    ALPHAVISCOSITY *
		    pow2(data[t_data::SOUNDSPEED](n_radial, n_azimuthal)) *
		    inv_omega_kepler;
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
    double drr, dpp, drp;

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
	    drr = (data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) -
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
	    dpp = (data[t_data::V_AZIMUTHAL](
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
	    data[t_data::TAU_PHI_PHI](n_radial, n_azimuthal) =
		2.0 * data[t_data::VISCOSITY](n_radial, n_azimuthal) *
		data[t_data::DENSITY](n_radial, n_azimuthal) *
		(dpp - 1.0 / 3.0 * data[t_data::DIV_V](n_radial, n_azimuthal));
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
	    drp = Ra[n_radial] * dvphirdr + dvrdphi * InvRa[n_radial];

	    unsigned int n_azimuthal_minus =
		(n_azimuthal == 0 ? data[t_data::VISCOSITY].get_max_azimuthal()
				  : n_azimuthal - 1);

	    // averaged nu over 4 corresponding cells
	    double nu =
		0.25 *
		(data[t_data::VISCOSITY](n_radial, n_azimuthal) +
		 data[t_data::VISCOSITY](n_radial - 1, n_azimuthal) +
		 data[t_data::VISCOSITY](n_radial, n_azimuthal_minus) +
		 data[t_data::VISCOSITY](n_radial - 1, n_azimuthal_minus));

	    // averaged sigma over 4 corresponding cells
	    double sigma =
		0.25 * (data[t_data::DENSITY](n_radial, n_azimuthal) +
			data[t_data::DENSITY](n_radial - 1, n_azimuthal) +
			data[t_data::DENSITY](n_radial, n_azimuthal_minus) +
			data[t_data::DENSITY](n_radial - 1, n_azimuthal_minus));

	    // tau_r_phi = nu*Sigma*( r*d(v_phi/r)/dr + 1/r d(v_r)/dphi )
	    data[t_data::TAU_R_PHI](n_radial, n_azimuthal) = nu * sigma * drp;
	}
    }

    if (include_artifical_viscosity) {
	for (unsigned int n_radial = 0;
	     n_radial <= data[t_data::DIV_V].get_max_radial(); ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::DIV_V].get_max_azimuthal();
		 ++n_azimuthal) {
		double nu_art;

		if (data[t_data::DIV_V](n_radial, n_azimuthal) < 0) {
		    nu_art =
			parameters::artificial_viscosity_factor *
			data[t_data::DENSITY](n_radial, n_azimuthal) *
			pow2(min(
			    Rsup[n_radial] - Rinf[n_radial],
				Rmed[n_radial] * 2 * M_PI /
				data[t_data::DENSITY].get_size_azimuthal())) *
			(-data[t_data::DIV_V](n_radial, n_azimuthal));
		} else {
		    nu_art = 0;
		}

		data[t_data::TAU_R_R](n_radial, n_azimuthal) +=
		    nu_art * data[t_data::DIV_V](n_radial, n_azimuthal);
		data[t_data::TAU_PHI_PHI](n_radial, n_azimuthal) +=
		    nu_art * data[t_data::DIV_V](n_radial, n_azimuthal);
	    }
	}
    }
}

/**
	Update velocities with viscous source term of Navier-Stokes equations
*/
void update_velocities_with_viscosity(t_data &data, t_polargrid &v_radial,
				      t_polargrid &v_azimuthal, double dt)
{
    double invdphi;

    invdphi = 1.0 / (2.0 * M_PI / (double)data[t_data::DENSITY].Nsec);

    double n_azimuthal_plus, n_azimuthal_minus;

    for (unsigned int n_radial = 1; n_radial <= v_radial.get_max_radial() - 1;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= v_radial.get_max_azimuthal(); ++n_azimuthal) {
	    n_azimuthal_plus =
		(n_azimuthal == data[t_data::DENSITY].get_max_azimuthal()
		     ? 0
		     : n_azimuthal + 1);
	    n_azimuthal_minus =
		(n_azimuthal == 0 ? data[t_data::DENSITY].get_max_azimuthal()
				  : n_azimuthal - 1);

	    // a_phi = 1/(r*Sigma) ( d(r*tau_r_phi)/dr + d(tau_phi_phi)/dphi +
	    // tau_r_phi )
	    /*			v_azimuthal(n_radial,n_azimuthal) +=
	       dt*InvRmed[n_radial]/((data[t_data::DENSITY](n_radial,n_azimuthal)+data[t_data::DENSITY](n_radial,n_azimuthal_minus))/2.0)
								       *((Rsup[n_radial]*data[t_data::TAU_R_PHI](n_radial+1,n_azimuthal)-Rinf[n_radial]*data[t_data::TAU_R_PHI](n_radial,n_azimuthal))*InvDiffRsup[n_radial]
								       +
	       (data[t_data::TAU_PHI_PHI](n_radial,n_azimuthal)-data[t_data::TAU_PHI_PHI](n_radial,n_azimuthal_minus))*invdphi
								       +
	       0.5*(data[t_data::TAU_R_PHI](n_radial,n_azimuthal)+data[t_data::TAU_R_PHI](n_radial+1,n_azimuthal)));
	     */

	    // a_phi = 1/(r*Sigma) ( 1/r d(r^2*tau_r_phi)/dr +
	    // d(tau_phi_phi)/dphi )
	    /*			v_azimuthal(n_radial,n_azimuthal) +=
	       dt*InvRmed[n_radial]/(0.5*(data[t_data::DENSITY](n_radial,n_azimuthal)+data[t_data::DENSITY](n_radial,n_azimuthal_minus)))
								       *(InvRmed[n_radial]*(pow2(Rsup[n_radial])*data[t_data::TAU_R_PHI](n_radial+1,n_azimuthal)-pow2(Rinf[n_radial])*data[t_data::TAU_R_PHI](n_radial,n_azimuthal))*InvDiffRsup[n_radial]
								       +
	       (data[t_data::TAU_PHI_PHI](n_radial,n_azimuthal)-data[t_data::TAU_PHI_PHI](n_radial,n_azimuthal_minus))*invdphi);
	    */
	    double sigma_avg =
		0.5 * (data[t_data::DENSITY](n_radial, n_azimuthal) +
		       data[t_data::DENSITY](n_radial, n_azimuthal_minus));
	    v_azimuthal(n_radial, n_azimuthal) +=
		dt * InvRb[n_radial] / (sigma_avg) *
		((2.0 / (pow2(Ra[n_radial + 1]) - pow2(Ra[n_radial]))) *
		     (pow2(Ra[n_radial + 1]) *
			  data[t_data::TAU_R_PHI](n_radial + 1, n_azimuthal) -
		      pow2(Ra[n_radial]) *
			  data[t_data::TAU_R_PHI](n_radial, n_azimuthal)) +
		 (data[t_data::TAU_PHI_PHI](n_radial, n_azimuthal) -
		  data[t_data::TAU_PHI_PHI](n_radial, n_azimuthal_minus)) *
		     invdphi);

	    // a_r = 1/(r*Sigma) ( d(r*tau_r_r)/dr + d(tau_r_phi)/dphi -
	    // tau_phi_phi )
	    sigma_avg =
		0.5 * (data[t_data::DENSITY](n_radial, n_azimuthal) +
		       data[t_data::DENSITY](n_radial - 1, n_azimuthal));
	    v_radial(n_radial, n_azimuthal) +=
		dt * InvRinf[n_radial] /
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
	}
    }
}

} // namespace viscosity
