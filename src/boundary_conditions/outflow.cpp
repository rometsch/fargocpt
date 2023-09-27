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

void outflow_inner(t_polargrid &vr, [[maybe_unused]] t_polargrid &dummy,
		   [[maybe_unused]] t_data &ddummy)
{
    if (CPU_Rank != 0) {
	return;
    }

    const unsigned int Iaz = vr.get_max_azimuthal();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
	if (vr(2, naz) > 0.0) {
	    vr(1, naz) = 0.0;
	    vr(0, naz) = 0.0;
	} else {
	    vr(1, naz) = vr(2, naz);
	    vr(0, naz) = vr(2, naz);
	}
    }
}

void outflow_outer(t_polargrid &vr, [[maybe_unused]] t_polargrid &dummy,
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
	    vr(Irad - 1, naz) = 0.0;
	    vr(Irad, naz) = 0.0;
	} else {
	    vr(Irad - 1, naz) = vr(Irad - 2, naz);
	    vr(Irad, naz) = vr(Irad - 2, naz);
	}
    }
}

/**
	inner open boundary condition
*/
void open_boundary_inner(t_data &data)
{
    if (CPU_Rank != 0) {
	return;
    }

    t_polargrid &Sigma = data[t_data::SIGMA];
    t_polargrid &Energy = data[t_data::ENERGY];
    t_polargrid &vrad = data[t_data::V_RADIAL];

    const unsigned int Naz = Sigma.get_max_azimuthal();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {
	// copy first ring into ghost ring
	Sigma(0, naz) = Sigma(1, naz);
	Energy(0, naz) = Energy(1, naz);

	// set velocity to min(v[+1],0) (allow only outflow)
	if (vrad(2, naz) > 0.0) {
	    vrad(1, naz) = 0.0;
	    vrad(0, naz) = 0.0;
	} else {
	    vrad(1, naz) = vrad(2, naz);
	    vrad(0, naz) = vrad(2, naz);
	}
    }
}

/**
	outer open boundary condition
 */
void open_boundary_outer(t_data &data)
{
    if (CPU_Rank != CPU_Highest) {
	return;
    }

    t_polargrid &Sigma = data[t_data::SIGMA];
    t_polargrid &Energy = data[t_data::ENERGY];
    t_polargrid &vrad = data[t_data::V_RADIAL];

    const unsigned int Nrad = Sigma.get_max_radial();
    const unsigned int Nrad_v = vrad.get_max_radial();
    const unsigned int Naz = Sigma.get_max_azimuthal();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {
	// copy last ring into ghost ring
	Sigma(Nrad, naz) = Sigma(Nrad - 1, naz);
	Energy(Nrad, naz) = Energy(Nrad - 1, naz);

	// set velocity to min(v[+1],0) (allow only outflow)
	if (vrad(Nrad_v - 2, naz) < 0.0) {
	    vrad(Nrad_v - 1, naz) = 0.0;
	    vrad(Nrad_v, naz) = 0.0;
	} else {
	    vrad(Nrad_v - 1, naz) = vrad(Nrad_v - 2, naz);
	    vrad(Nrad_v, naz) = vrad(Nrad_v - 2, naz);
	}
    }
}

} // namespace boundary_conditions
