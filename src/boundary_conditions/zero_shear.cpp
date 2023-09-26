/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"

#include "../global.h"


namespace boundary_conditions
{

/****************************************
/ Enforce zero shear at the boundaries
****************************************/
void zero_shear_inner(t_polargrid &vaz, [[maybe_unused]] t_polargrid &dummy, [[maybe_unused]] t_data &ddummy) {
	/// Vphi_i / r_i - Vphi_(i-1) / r_(i-1) = 0
	const unsigned int Iaz = vaz.get_max_azimuthal();

	if (CPU_Rank == 0) {
		#pragma omp parallel for
	    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
			const double v_active = vaz(1, naz);
			const double r_active = Rmed[1];
			const double Omega_active = v_active / r_active;
			const double r = Rmed[0];
			vaz(0, naz) = r * Omega_active;
	    }
	}

}

void zero_shear_outer(t_polargrid &vaz, [[maybe_unused]] t_polargrid &dummy, [[maybe_unused]] t_data &ddummy) {
	const unsigned int Iaz = vaz.get_max_azimuthal();
	const unsigned int Irad = vaz.get_max_radial();

	if (CPU_Rank == CPU_Highest) {
		#pragma omp parallel for
	    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
			const double v_active = vaz(Irad - 1, naz);
			const double r_active = Rmed[Irad - 1];
			const double Omega_active = v_active / r_active;
			const double r = Rmed[Irad];
			vaz(Irad, naz) = r * Omega_active;
	    }
	}

}
} // namespace boundary_conditions
