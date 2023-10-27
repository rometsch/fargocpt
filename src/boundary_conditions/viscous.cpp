/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "../global.h"
#include "../parameters.h"
#include "boundary_conditions.h"

namespace boundary_conditions
{

void viscous_outflow_inner(t_polargrid &vr, [[maybe_unused]] t_polargrid &dummy,
			   t_data &data)
{
    if (CPU_Rank != 0) {
	return;
    }

    const unsigned int Iaz = vr.get_max_azimuthal();
    const double s = viscous_outflow_speed;

    const t_polargrid &visc = data[t_data::VISCOSITY];

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {

	const double Nu0 = visc(0, naz);
	const double Nu1 = visc(1, naz);
	const double Nu = 0.5 * (Nu0 + Nu1);

	// V_rad =  - 1.5 / r * Nu (Kley, Papaloizou and Ogilvie, 2008)
	vr(1, naz) = -1.5 * s / Rinf[1] * Nu;
	vr(0, naz) = -1.5 * s / Rinf[0] * Nu;
    }
}

void viscous_inflow_outer(t_polargrid &vr, [[maybe_unused]] t_polargrid &dummy,
			  t_data &data)
{
    if (CPU_Rank != 0) {
	return;
    }

    const unsigned int Iaz = vr.get_max_azimuthal();
    const unsigned int Irad = vr.get_max_radial();
    const double s = viscous_outflow_speed;

    const t_polargrid &visc = data[t_data::VISCOSITY];

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {

	const double Nu0 = visc(Irad, naz);
	const double Nu1 = visc(Irad - 1, naz);
	const double Nu = 0.5 * (Nu0 + Nu1);

	// V_rad =  - 1.5 / r * Nu (Kley, Papaloizou and Ogilvie, 2008)
	vr(Irad, naz) = -1.5 * s / Rinf[Irad] * Nu;
	vr(Irad - 1, naz) = -1.5 * s / Rinf[Irad - 1] * Nu;
    }
}

} // namespace boundary_conditions
