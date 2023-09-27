/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "../global.h"
#include "boundary_conditions.h"

namespace boundary_conditions
{

static void zero_gradient_inner_scalar(t_polargrid &x)
{

    const unsigned int Iaz = x.get_max_azimuthal();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
	x(0, naz) = x(1, naz);
    }
}

static void zero_gradient_inner_vector(t_polargrid &x)
{

    const unsigned int Iaz = x.get_max_azimuthal();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
	x(0, naz) = x(2, naz);
	x(1, naz) = x(2, naz);
    }
}

void zero_gradient_inner(t_polargrid &x, [[maybe_unused]] t_polargrid &dummy,
			 [[maybe_unused]] t_data &ddummy)
{
    if (CPU_Rank != 0) {
	return;
    }
    if (x.is_scalar()) {
	zero_gradient_inner_scalar(x);
    } else {
	zero_gradient_inner_vector(x);
    }
}

static void zero_gradient_outer_scalar(t_polargrid &x)
{

    const unsigned int Iaz = x.get_max_azimuthal();
    const unsigned int Irad = x.get_max_radial();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
	x(Irad, naz) = x(Irad - 1, naz);
    }
}

static void zero_gradient_outer_vector(t_polargrid &x)
{

    const unsigned int Iaz = x.get_max_azimuthal();
    const unsigned int Irad = x.get_max_radial();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
	x(Irad, naz) = x(Irad - 2, naz);
	x(Irad - 1, naz) = x(Irad - 2, naz);
    }
}

void zero_gradient_outer(t_polargrid &x, [[maybe_unused]] t_polargrid &dummy,
			 [[maybe_unused]] t_data &ddummy)
{
    if (CPU_Rank != CPU_Highest) {
	return;
    }
    if (x.is_scalar()) {
	zero_gradient_outer_scalar(x);
    } else {
	zero_gradient_outer_vector(x);
    }
}

} // namespace boundary_conditions
