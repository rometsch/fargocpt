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

void zero_gradient_boundary_inner(t_data &data)
{
    if (CPU_Rank != 0){
		return;
	}

	t_polargrid &sig = data[t_data::SIGMA];
	t_polargrid &e = data[t_data::ENERGY];
	t_polargrid &vr = data[t_data::V_RADIAL];

	const unsigned int Naz = sig.get_max_azimuthal();

	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {
		// copy first ring into ghost ring
		sig(0, naz) = sig(1, naz);
		e(0, naz) = e(1, naz);

		vr(1, naz) = vr(2, naz);
		vr(0, naz) = vr(2, naz);
    }
}

void zero_gradient_boundary_outer(t_data &data)
{
    if (CPU_Rank != CPU_Highest) {
		return;
	}

	t_polargrid &sig = data[t_data::SIGMA];
	t_polargrid &e = data[t_data::ENERGY];
	t_polargrid &vr = data[t_data::V_RADIAL];

	const unsigned int Nrad = sig.get_max_radial();
	const unsigned int Nrad_v = vr.get_max_radial();
	const unsigned int Naz = sig.get_max_azimuthal();


	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {
	// copy last ring into ghost ring
	sig(Nrad, naz) = sig(Nrad - 1, naz);
	e(Nrad, naz) = e(Nrad - 1, naz);

	vr(Nrad_v - 1, naz) =
	    vr(Nrad_v - 2, naz);
	vr(Nrad_v, naz) = vr(Nrad_v - 2, naz);
    }
}

} // namespace boundary_conditions
