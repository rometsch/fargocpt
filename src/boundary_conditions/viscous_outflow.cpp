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

namespace boundary_conditions
{

void viscous_outflow_boundary_inner(t_data &data)
{
    if (CPU_Rank != 0) {
		return;
	}

	t_polargrid &sig = data[t_data::SIGMA];
	t_polargrid &e = data[t_data::ENERGY];
	t_polargrid &vr = data[t_data::V_RADIAL];
	t_polargrid &visc = data[t_data::VISCOSITY];

	const unsigned int Naz = sig.get_max_azimuthal();

	for (unsigned int naz = 0; naz <= Naz; ++naz) {
	    sig(0, naz) = sig(1, naz);
	    e(0, naz) = e(1, naz);
	}

	const double s = parameters::viscous_outflow_speed;

	#pragma omp parallel for
	for (unsigned int naz = 0; naz <= Naz; ++naz) {

	    const double Nu0 = visc(0, naz);
	    const double Nu1 = visc(1, naz);
	    const double Nu = 0.5 * (Nu0 + Nu1);

	    // V_rad =  - 1.5 / r * Nu (Kley, Papaloizou and Ogilvie, 2008)
	    vr(1, naz) = -1.5 * s / Rinf[1] * Nu;
	    vr(0, naz) = -1.5 * s / Rinf[0] * Nu;
	}
}

} // namespace boundary_conditions
