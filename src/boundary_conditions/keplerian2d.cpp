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
#include "../constants.h"


namespace boundary_conditions
{

static double calc_sig(const double R) {
	const double sig0 = parameters::sigma0;
	const double sigexp = -parameters::SIGMASLOPE;
	return sig0 * std::pow(R, sigexp);
}

static double calc_e(const double R) {
	const double gamma = parameters::ADIABATICINDEX;
	const double sig0 = parameters::sigma0;
	const double sigexp = -parameters::SIGMASLOPE;
	const double f = parameters::FLARINGINDEX;
	const double h0 = parameters::ASPECTRATIO_REF;
	return 1.0 / (gamma - 1.0) * sig0 * std::pow(h0, 2) * std::pow(R, sigexp - 1.0 + 2.0 * f);
}

static double calc_T(const double e, const double sig) {
	const double gamma = parameters::ADIABATICINDEX;
	return e / sig * (gamma - 1.0) * parameters::MU * constants::R;
}

static double calc_vK(const double R) {
	return std::sqrt(constants::G * hydro_center_mass / R);
}

void keplerian2d_boundary_inner(t_data &data)
{

	auto &sig = data[t_data::SIGMA];
	auto &e = data[t_data::ENERGY];
	auto &T = data[t_data::TEMPERATURE];
	auto &vrad = data[t_data::V_RADIAL];
	auto &vaz = data[t_data::V_AZIMUTHAL];

	const unsigned int Naz = sig.get_max_azimuthal();

	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {
		sig(0, naz) = calc_sig(Rmed[0]);

		e(0, naz) = calc_e(Rmed[0]);

		T(0, naz) = calc_T(e(0, naz), sig(0, naz));

		vrad(1, naz) = 0.0;
		vrad(0, naz) = -vrad(2, naz);

		vaz(0, naz) = calc_vK(Rmed[0]);
    }
}

void keplerian2d_boundary_outer(t_data &data)
{
	auto &sig = data[t_data::SIGMA];
	auto &e = data[t_data::ENERGY];
	auto &T = data[t_data::TEMPERATURE];
	auto &vrad = data[t_data::V_RADIAL];
	auto &vaz = data[t_data::V_AZIMUTHAL];

	const unsigned int Naz = sig.get_max_azimuthal();
	const unsigned int Nrad = sig.get_max_radial();
	const unsigned int Nrad_v = vrad.get_max_radial();

	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {
		sig(Nrad, naz) = calc_sig(Rmed[Nrad]);

		e(Nrad, naz) = calc_e(Rmed[Nrad]);

		T(Nrad, naz) = calc_T(e(Nrad, naz), sig(Nrad, naz));

		vrad(Nrad_v, naz) = -vrad(Nrad_v - 2, naz);
		vrad(Nrad_v - 1, naz) = 0.0;

		vaz(Nrad, naz) = calc_vK(Rmed[Nrad]);
	}
}


} // namespace boundary_conditions
