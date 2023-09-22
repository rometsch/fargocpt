/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"
#include "../global.h"
#include "../simulation.h"


namespace boundary_conditions
{

/**
	set inner boundary values to initial values
 */
void reference_value_boundary_inner(t_data &data)
{
    if (CPU_Rank != 0){
		return;
	}
	if (sim::N_hydro_iter == 0) {
		return;
	}

    auto &sig = data[t_data::SIGMA];
    auto &eng = data[t_data::ENERGY];
	auto &vr = data[t_data::V_RADIAL];

    auto &sig0 = data[t_data::SIGMA0];
    auto &eng0 = data[t_data::ENERGY0];
	auto &vr0 = data[t_data::V_RADIAL0];

	const unsigned int Naz = sig.get_max_azimuthal();

	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {
		sig(0, naz) = sig0(0,naz);
		eng(0, naz) = eng0(0,naz);
		vr(0, naz) = vr0(0,naz);
		vr(1, naz) = vr0(1,naz);
	}
}

/**
	set outer boundary values to initial values
 */
void reference_value_boundary_outer(t_data &data)
{
    if (CPU_Rank != CPU_Highest){
		return;
	}
	if (sim::N_hydro_iter == 0) {
		return;
	}

    auto &sig = data[t_data::SIGMA];
    auto &eng = data[t_data::ENERGY];
	auto &vr = data[t_data::V_RADIAL];

    auto &sig0 = data[t_data::SIGMA0];
    auto &eng0 = data[t_data::ENERGY0];
	auto &vr0 = data[t_data::V_RADIAL0];

	const unsigned int Naz = sig.get_max_azimuthal();

	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {
		const unsigned int nr = sig.get_max_radial();
		sig(nr, naz) = sig0(nr,naz);
		eng(nr, naz) = eng0(nr,naz);
		const unsigned int nrv = vr.get_max_radial();
		vr(nrv, naz) = vr0(nrv,naz);
		vr(nrv-1, naz) = vr0(nrv-1,naz);
	}
}


} // namespace boundary_conditions
