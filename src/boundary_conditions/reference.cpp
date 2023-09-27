/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "../global.h"
#include "../simulation.h"
#include "boundary_conditions.h"

namespace boundary_conditions
{

static void reference_inner_scalar(t_polargrid &x, t_polargrid &x0)
{

    const unsigned int Iaz = x.get_max_azimuthal();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
	x(0, naz) = x0(0, naz);
    }
}

static void reference_inner_vector(t_polargrid &x, t_polargrid &x0)
{

    const unsigned int Iaz = x.get_max_azimuthal();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
	x(0, naz) = x0(0, naz);
	x(1, naz) = x0(1, naz);
    }
}

void reference_inner(t_polargrid &x, [[maybe_unused]] t_polargrid &x0,
		     [[maybe_unused]] t_data &ddummy)
{
    if (CPU_Rank != 0) {
	return;
    }
    if (x.is_scalar()) {
	reference_inner_scalar(x, x0);
    } else {
	reference_inner_vector(x, x0);
    }
}

static void reference_outer_scalar(t_polargrid &x, t_polargrid &x0)
{

    const unsigned int Iaz = x.get_max_azimuthal();
    const unsigned int Irad = x.get_max_radial();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
	x(Irad, naz) = x0(Irad, naz);
    }
}

static void reference_outer_vector(t_polargrid &x, t_polargrid &x0)
{

    const unsigned int Iaz = x.get_max_azimuthal();
    const unsigned int Irad = x.get_max_radial();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
	x(Irad, naz) = x0(Irad, naz);
	x(Irad - 1, naz) = x0(Irad - 1, naz);
    }
}

void reference_outer(t_polargrid &x, [[maybe_unused]] t_polargrid &x0,
		     [[maybe_unused]] t_data &ddummy)
{
    if (CPU_Rank != CPU_Highest) {
	return;
    }
    if (x.is_scalar()) {
	reference_outer_scalar(x, x0);
    } else {
	reference_outer_vector(x, x0);
    }
}

} // namespace boundary_conditions
