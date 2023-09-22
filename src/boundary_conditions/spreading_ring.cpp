/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"
#include "../global.h"
#include "../Theo.h"
#include "../frame_of_reference.h"

namespace boundary_conditions
{

/**
	for viscous spreading ring comparison simulations for Jibin Joseph
*/
void spreading_ring_inner(t_data &data)
{
    if (CPU_Rank != 0) {
		return;
	}

    const double R = Rmed[0];

	const double vaz_new = initial_locally_isothermal_smoothed_v_az(R, 1.0) - R * refframe::OmegaFrame;

	auto &vrad = data[t_data::V_RADIAL];
	auto &vaz = data[t_data::V_AZIMUTHAL];
	auto &sig = data[t_data::SIGMA];

	const unsigned int Naz = sig.get_max_azimuthal();

	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {

	// copy first ring into ghost ring
	sig(0, naz) = sig(1, naz);
	vaz(0, naz) = vaz_new;

	if (vrad(2, naz) <= 0.0) { // outflow
	    vrad(1, naz) = vrad(2, naz);
	    vrad(0, naz) = vrad(2, naz);
	} else { // reflective
	    vrad(1, naz) = -vrad(2, naz);
	    vrad(0, naz) = -vrad(2, naz);
	}
    }
}

/**
	for viscous spreading ring comparison simulations for Jibin Joseph
*/
void spreading_ring_outer(t_data &data)
{
    if (CPU_Rank != CPU_Highest) {
		return;
	}

	auto &vaz = data[t_data::V_AZIMUTHAL];

	const double R = Rmed[data[t_data::V_AZIMUTHAL].get_max_radial()];

	const double vaz_new = initial_locally_isothermal_smoothed_v_az(R, 1.0) - R * refframe::OmegaFrame;

	const unsigned int Nrad = vaz.get_max_radial();
	const unsigned int Naz = vaz.get_max_azimuthal();

	#pragma omp parallel for
	for (unsigned int naz = 0; naz <= Naz; ++naz) {
		vaz(Nrad, naz) = vaz_new;
	}
}

} // namespace boundary_conditions
