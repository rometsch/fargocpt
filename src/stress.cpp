#include "stress.h"
#include "constants.h"
#include "global.h"
#include "selfgravity.h"
#include "parameters.h"
#include "logging.h"

namespace stress
{

void calculate_gravitational_stress(t_data &data)
{
	#ifdef DISABLE_FFTW
		logging::print_master(LOG_ERROR
			"Self-gravity is not compiled in. Please recompile with FFTW enabled.\n");
		PersonalExit(1);
	#endif
	const unsigned int Nr = data[t_data::T_GRAVITATIONAL].get_size_radial();
	const unsigned int Nphi = data[t_data::T_GRAVITATIONAL].get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
	    // TODO: Factor sqrt(2*PI) instead of 2 should be more accurate
		data[t_data::T_GRAVITATIONAL](nr, naz) =
		1.0 / (4.0 * M_PI * constants::G) *
		selfgravity::g_radial[nr * NAzimuthal + naz] *
		selfgravity::g_azimuthal[nr * NAzimuthal + naz] *
		(2 * parameters::aspectratio_ref * Rmed[nr]);
	}
    }
}

void calculate_Reynolds_stress(t_data &data)
{

	const unsigned int Nr = data[t_data::T_REYNOLDS].get_size_radial();
	const unsigned int Nphi = data[t_data::V_AZIMUTHAL].get_size_azimuthal();

    // T_Reynolds = Sigma*(v_radial - v_radial_mean)*(v_azimuthal -
    // v_azimuthal_mean)
	#pragma omp parallel for
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	// calculate v_radial_mean & v_azimuthal_mean
	double v_radial_mean = 0.0;
	double v_azimuthal_mean = 0.0;

	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const unsigned int naz_next = (naz == Nphi-1 ? 0 : naz+1);
	    v_radial_mean +=
		0.5 * (data[t_data::V_RADIAL](nr, naz) +
			   data[t_data::V_RADIAL](nr + 1, naz));

		v_azimuthal_mean += 0.5*(data[t_data::V_AZIMUTHAL](nr, naz)
				+ data[t_data::V_AZIMUTHAL](nr, naz_next));
	}
	v_azimuthal_mean /= (double)Nphi;
	v_radial_mean /= (double)Nphi;

	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const unsigned int naz_next = (naz == Nphi-1 ? 0 : naz+1);
		data[t_data::T_REYNOLDS](nr, naz) =
		data[t_data::SIGMA](nr, naz) *
		(0.5 * (data[t_data::V_RADIAL](nr, naz) +
			data[t_data::V_RADIAL](nr + 1, naz)) - v_radial_mean) *
		(0.5 * (data[t_data::V_AZIMUTHAL](nr, naz) +
		 data[t_data::V_AZIMUTHAL](nr, naz_next)) - v_azimuthal_mean);
	}
    }
}

} // namespace stress
