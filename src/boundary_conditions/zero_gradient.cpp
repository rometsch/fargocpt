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

static void zero_gradient_inner_scalar(t_polargrid &x) {

	const unsigned int Iaz = x.get_max_azimuthal();

	#pragma omp parallel for
	for (unsigned int naz = 0; naz <= Iaz; ++naz) {
		x(0, naz) = x(1, naz);
	}
}

static void zero_gradient_inner_vector(t_polargrid &x) {
	
	const unsigned int Iaz = x.get_max_azimuthal();

	#pragma omp parallel for
	for (unsigned int naz = 0; naz <= Iaz; ++naz) {
		x(0, naz) = x(2, naz);
		x(1, naz) = x(2, naz);
	}
}


void zero_gradient_inner(t_polargrid &x, [[maybe_unused]] t_polargrid &dummy, [[maybe_unused]] t_data &ddummy) {
	if (CPU_Rank != 0) {
		return;
	}
	if (x.is_scalar()) {
		zero_gradient_inner_scalar(x);
	} else {
		zero_gradient_inner_vector(x);
	}
}

static void zero_gradient_outer_scalar(t_polargrid &x) {

    const unsigned int Iaz = x.get_max_azimuthal();
    const unsigned int Irad = x.get_max_radial();

    #pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
        x(Irad, naz) = x(Irad-1, naz);
    }
}

static void zero_gradient_outer_vector(t_polargrid &x) {

	const unsigned int Iaz = x.get_max_azimuthal();
	const unsigned int Irad = x.get_max_radial();

	#pragma omp parallel for
	for (unsigned int naz = 0; naz <= Iaz; ++naz) {
		x(Irad, naz) = x(Irad-2, naz);
		x(Irad-1, naz) = x(Irad-2, naz);
	}
}


void zero_gradient_outer(t_polargrid &x, [[maybe_unused]] t_polargrid &dummy, [[maybe_unused]] t_data &ddummy) {
    if (CPU_Rank != CPU_Highest) {
        return;
    }
    if (x.is_scalar()) {
        zero_gradient_outer_scalar(x);
    } else {
        zero_gradient_outer_vector(x);
    }
}

void zero_gradient_boundary_inner(t_data &data)
{
    if (CPU_Rank != 0){
		return;
	}

	t_polargrid &sig = data[t_data::SIGMA];
	t_polargrid &eng = data[t_data::ENERGY];
	t_polargrid &vrad = data[t_data::V_RADIAL];

	const unsigned int Naz = sig.get_max_azimuthal();

	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {
		// copy first ring into ghost ring
		sig(0, naz) = sig(1, naz);
		eng(0, naz) = eng(1, naz);

		vrad(1, naz) = vrad(2, naz);
		vrad(0, naz) = vrad(2, naz);
    }
}

void zero_gradient_boundary_outer(t_data &data)
{
    if (CPU_Rank != CPU_Highest) {
		return;
	}

	t_polargrid &sig = data[t_data::SIGMA];
	t_polargrid &eng = data[t_data::ENERGY];
	t_polargrid &vrad = data[t_data::V_RADIAL];

	const unsigned int Nrad = sig.get_max_radial();
	const unsigned int Nrad_v = vrad.get_max_radial();
	const unsigned int Naz = sig.get_max_azimuthal();


	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {
	// copy last ring into ghost ring
	sig(Nrad, naz) = sig(Nrad - 1, naz);
	eng(Nrad, naz) = eng(Nrad - 1, naz);

	vrad(Nrad_v - 1, naz) = vrad(Nrad_v - 2, naz);
	vrad(Nrad_v, naz) = vrad(Nrad_v - 2, naz);
    }
}

} // namespace boundary_conditions
