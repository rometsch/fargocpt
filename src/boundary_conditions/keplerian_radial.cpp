/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"

#include "../Theo.h"
#include "../frame_of_reference.h"
#include "../global.h"
#include "../parameters.h"

namespace boundary_conditions
{

void keplerian_radial_inner(t_polargrid &vaz,
			    [[maybe_unused]] t_polargrid &dummy,
			    [[maybe_unused]] t_data &ddummy)
{

    if (CPU_Rank != 0) {
	return;
    }

    const unsigned int Iaz = vaz.get_max_azimuthal();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {

	for (unsigned int k = 0; k <= 1; k++) {
	    // set first two rings
	    const double vKep = compute_v_kepler(Rmed[k], hydro_center_mass);
	    const double OmegaF = refframe::OmegaFrame;
	    const double r = Rmed[k];
	    const double c = keplerian_radial_inner_factor;
	    vaz(k, naz) = c * vKep;
	}
    }
}

void keplerian_radial_outer(t_polargrid &vaz,
			    [[maybe_unused]] t_polargrid &dummy,
			    [[maybe_unused]] t_data &ddummy)
{

    if (CPU_Rank != CPU_Highest) {
	return;
    }

    const unsigned int Iaz = vaz.get_max_azimuthal();
    const unsigned int Irad = vaz.get_max_radial();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {

	for (unsigned int k = Irad; k >= Irad - 1; k--) {
	    // set last two rings
	    const double vKep = compute_v_kepler(Rmed[k], hydro_center_mass);
	    const double OmegaF = refframe::OmegaFrame;
	    const double r = Rmed[k];
	    const double c = keplerian_radial_outer_factor;
	    vaz(k, naz) = c * vKep;
	}
    }
}

} // namespace boundary_conditions
