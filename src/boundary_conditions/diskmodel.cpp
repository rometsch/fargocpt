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

static double calc_sig(const double R) {
	const double sig0 = parameters::sigma0;
	const double sigexp = -parameters::SIGMASLOPE;
	return sig0 * std::pow(R, sigexp);
}

static double calc_eng(const double R) {
	const double gamma = parameters::ADIABATICINDEX;
	const double sig0 = parameters::sigma0;
	const double sigexp = -parameters::SIGMASLOPE;
	const double f = parameters::FLARINGINDEX;
	const double h0 = parameters::ASPECTRATIO_REF;
	return 1.0 / (gamma - 1.0) * sig0 * std::pow(h0, 2) * std::pow(R, sigexp - 1.0 + 2.0 * f);
}

void diskmodel_inner_sigma(t_polargrid &sig, [[maybe_unused]] t_polargrid &dummy,
			 [[maybe_unused]] t_data &ddummy) {

	const unsigned int Iaz = sig.get_max_azimuthal();

	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
		sig(0, naz) = calc_sig(Rmed[0]);
    }
}

void diskmodel_inner_energy(t_polargrid &eng, [[maybe_unused]] t_polargrid &dummy,
			 [[maybe_unused]] t_data &ddummy) {

	const unsigned int Iaz = eng.get_max_azimuthal();

	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
		eng(0, naz) = calc_eng(Rmed[0]);
    }
}

void diskmodel_outer_sigma(t_polargrid &sig, [[maybe_unused]] t_polargrid &dummy,
			 [[maybe_unused]] t_data &ddummy) {

	const unsigned int Iaz = sig.get_max_azimuthal();
	const unsigned int Irad = sig.get_max_radial();

	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
		sig(Irad, naz) = calc_sig(Rmed[Irad]);
    }
}

void diskmodel_outer_energy(t_polargrid &eng, [[maybe_unused]] t_polargrid &dummy,
			 [[maybe_unused]] t_data &ddummy) {

	const unsigned int Iaz = eng.get_max_azimuthal();
	const unsigned int Irad = eng.get_max_radial();

	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Iaz; ++naz) {
		eng(Irad, naz) = calc_eng(Rmed[Irad]);
    }
}

} // namespace boundary_conditions
