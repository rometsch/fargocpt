/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"

#include "../Theo.h"
#include "../global.h"
#include "../parameters.h"
#include "../frame_of_reference.h"

namespace boundary_conditions
{

void keplerian_azimuthal_inner(t_polargrid &vaz, [[maybe_unused]] t_polargrid &dummy, [[maybe_unused]] t_data &ddummy) {


	if (CPU_Rank != 0) {
		return;
	}
		
	const unsigned int Iaz = vaz.get_max_azimuthal();
	const double vKep = compute_v_kepler(Rmed[0], hydro_center_mass);
	const double OmegaF = refframe::OmegaFrame;
	const double r = Rmed[0];
	const double c = keplerian_azimuthal_inner_factor;

	#pragma omp parallel for
	for (unsigned int naz = 0; naz <= Iaz; ++naz) {
		vaz(0, naz) = c*vKep - r*OmegaF;
	}

	
}


void keplerian_azimuthal_outer(t_polargrid &vaz, [[maybe_unused]] t_polargrid &dummy, [[maybe_unused]] t_data &ddummy) {


	if (CPU_Rank != CPU_Highest) {
		return;
	}
		
	const unsigned int Iaz = vaz.get_max_azimuthal();
	const unsigned int Irad = vaz.get_max_radial();
	const double vKep = compute_v_kepler(Rmed[Irad], hydro_center_mass);
	const double OmegaF = refframe::OmegaFrame;
	const double r = Rmed[Irad];
	const double c = keplerian_azimuthal_outer_factor;

	#pragma omp parallel for
	for (unsigned int naz = 0; naz <= Iaz; ++naz) {
		vaz(Irad, naz) = c*vKep - r*OmegaF;
	}

	
}


} // namespace boundary_conditions
