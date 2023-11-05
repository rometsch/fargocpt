/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"

namespace boundary_conditions
{

void reflecting_inner(t_polargrid &vrad, [[maybe_unused]] t_polargrid &dummy,
		      [[maybe_unused]] t_data &ddummy)
{

    const unsigned int Iaz = vrad.get_max_azimuthal();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
	vrad(0, naz) = -vrad(2, naz);
	vrad(1, naz) = 0;
    }
}

void reflecting_outer(t_polargrid &vrad, [[maybe_unused]] t_polargrid &dummy,
		      [[maybe_unused]] t_data &ddummy)
{

    const unsigned int Iaz = vrad.get_max_azimuthal();
    const unsigned int Irad = vrad.get_max_radial();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
	vrad(Irad, naz) = -vrad(Irad - 2, naz);
	vrad(Irad - 1, naz) = 0;
    }
}

} // namespace boundary_conditions
