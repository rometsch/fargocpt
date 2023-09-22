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

namespace boundary_conditions
{

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

    if (CPU_Rank != 0) {
	return;
    }

    t_polargrid &Sigma = data[t_data::SIGMA];
    t_polargrid &Energy = data[t_data::ENERGY];
    t_polargrid &vrad = data[t_data::V_RADIAL];
    t_polargrid &vaz = data[t_data::V_AZIMUTHAL];

    const unsigned int Naz = Sigma.get_max_azimuthal();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {
	// zero gradient
	Sigma(0, naz) = Sigma(1, naz);
	Energy(0, naz) = Energy(1, naz);

	// set vrad to fraction of Keplerian velocity
	{
	    const double c = parameters::vrad_fraction_of_kepler;
	    const double vK = compute_v_kepler(hydro_center_mass, Rinf[1]);
	    vrad(1, naz) = -c * vK;
	    vrad(0, naz) = vrad(1, naz);
	}

	// set vphi to stellar rotation rate
	{
	    const double c = parameters::stellar_rotation_rate;
	    const double vK = compute_v_kepler(hydro_center_mass, Rmed[0]);
	    vaz(0, naz) = c * vK;
	}
    }
}

/**
	Outer boundary: floating boundary conditions & pressure correction for
   Omega
*/

void boundary_layer_outer_boundary(t_data &data)
{

    if (CPU_Rank != CPU_Highest) {
	return;
    }

    t_polargrid &Sigma = data[t_data::SIGMA];
    t_polargrid &Energy = data[t_data::ENERGY];
    t_polargrid &vrad = data[t_data::V_RADIAL];
    t_polargrid &vaz = data[t_data::V_AZIMUTHAL];

    const unsigned int Naz = Sigma.get_max_azimuthal();
    const unsigned int Nrad = Sigma.get_max_radial();
    const unsigned int Nrad_v = vrad.get_max_radial();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {
		// floating BCs
		Sigma(Nrad, naz) = Sigma(Nrad - 1, naz) * std::sqrt(Ra[Nrad - 1] / Ra[Nrad]);
		Energy(Nrad, naz) = Energy(Nrad - 1, naz) * std::pow(Ra[Nrad - 1] / Ra[Nrad], 1.25);

		vrad(Nrad_v - 1, naz) = -1. * fabs(vrad(Nrad_v - 2, naz)) * std::sqrt(Ra[Nrad_v - 2] / Ra[Nrad_v - 1]);
		vrad(Nrad_v, naz) = -1. * fabs(vrad(Nrad_v - 2, naz)) * std::sqrt(Ra[Nrad_v - 2] / Ra[Nrad_v]);

		// Omega at outer boundary equals calculate_omega_kepler (plus leading order pressure correction)
		vaz(Nrad, naz) = compute_v_kepler(hydro_center_mass, Rmed[Nrad]);
		// TODO: Include pressure correction, like in uphi[*jN]
		// = 1./sqrt(Rb[*jN]) +
		// 0.5/Sigma[*jN]*sqrt(pow3(Rb[*jN])*pow2(Rb[*jN]))*.5/Rb[*jN]*(P[*jN+1]-P[*jN-1])/DeltaRa[*jN+1];
	}
}

} // namespace boundary_conditions
