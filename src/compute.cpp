#include "data.h"
#include "parameters.h"
#include "compute.h"
#include "SourceEuler.h"
#include "opacity.h"

namespace compute {


/**
	computes density rho
*/
void midplane_density(t_data &data, const double current_time)
{
	compute_scale_height(data, current_time);

	auto &rho = data[t_data::RHO];
	auto &Sig = data[t_data::SIGMA];
	auto &H = data[t_data::SCALE_HEIGHT];

	const unsigned int Nr = rho.get_size_radial();
	const unsigned int Naz = rho.get_size_azimuthal();
	const double factor = parameters::density_factor;

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
		for (unsigned int naz = 0; naz < Naz; ++naz) {
			rho(nr, naz) = Sig(nr, naz) / (factor * H(nr, naz));
		}
    }
}

/**
 * Compute opacity values from current temperature and midplane density.
 * Store the values in data[KAPPA]
*/
void opacity(t_data &data) {
	auto &rho = data[t_data::RHO];
	auto &T = data[t_data::TEMPERATURE];
	auto &kappa = data[t_data::KAPPA];

	const unsigned int Nr = rho.get_size_radial();
	const unsigned int Naz = rho.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
		for (unsigned int naz = 0; naz < Naz; ++naz) {
			kappa(nr, naz) = opacity::opacity(rho(nr, naz), T(nr, naz));
		}
    }
}

} // namespace compute