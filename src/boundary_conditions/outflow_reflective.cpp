/**
	\file outflow_reflective.cpp

	Special conditions for the radial direction for the spreading ring test.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "../global.h"
#include "boundary_conditions.h"

namespace boundary_conditions
{

void outflow_reflective_inner(t_polargrid &vr, [[maybe_unused]] t_polargrid &dummy,
		   [[maybe_unused]] t_data &ddummy)
{
    if (CPU_Rank != 0) {
	return;
    }

    const unsigned int Iaz = vr.get_max_azimuthal();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
	if (vr(2, naz) > 0.0) {
	    vr(1, naz) = -vr(2, naz);
	    vr(0, naz) = -vr(2, naz);
	} else {
	    vr(1, naz) = vr(2, naz);
	    vr(0, naz) = vr(2, naz);
	}
    }
}

void outflow_reflective_outer(t_polargrid &vr, [[maybe_unused]] t_polargrid &dummy,
		   [[maybe_unused]] t_data &ddummy)
{
    if (CPU_Rank != CPU_Highest) {
	return;
    }

    const unsigned int Iaz = vr.get_max_azimuthal();
    const unsigned int Irad = vr.get_max_radial();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
	if (vr(Irad - 2, naz) < 0.0) {
	    vr(Irad - 1, naz) = -vr(Irad - 2, naz);
	    vr(Irad, naz) = -vr(Irad - 2, naz);
	} else {
	    vr(Irad - 1, naz) = vr(Irad - 2, naz);
	    vr(Irad, naz) = vr(Irad - 2, naz);
	}
    }
}

} // namespace boundary_conditions
