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
void kappa_eff(t_data &data) {
    auto &rho = data[t_data::RHO];
    auto &T = data[t_data::TEMPERATURE];
    auto &kappa = data[t_data::KAPPA];
    auto &tau = data[t_data::TAU];
    auto &tau_eff = data[t_data::TAU_EFF];
    auto &sigma = data[t_data::SIGMA];

    const unsigned int Nr = rho.get_size_radial();
    const unsigned int Naz = rho.get_size_azimuthal();

    #pragma omp parallel for collapse(2)
    for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Naz; ++naz) {
	    kappa(nr, naz) = opacity::opacity(rho(nr, naz), T(nr, naz));


	    // mean vertical optical depth: tau = 1/2 kappa Sigma
	    tau(nr, naz) =
		parameters::tau_factor *
		(1.0 / parameters::density_factor) *
		kappa(nr, naz) * sigma(nr, naz);

	    if(parameters::heating_star_enabled){
		//  irradiated disk tau_eff = 3/8 tau + 1/2 + 1/(4*tau+tau_min)
		//  compare D'Angelo & Marzari 2012
		tau_eff(nr, naz) =
		    3.0 / 8.0 * tau(nr, naz) +
		    0.5 +
		    1.0 / (4.0 * tau(nr, naz) + parameters::tau_min);
	    } else {
		//  non irradiated disk tau_eff = 3/8 tau + sqrt(3)/4 + 1/(4*tau+tau_min)
		tau_eff(nr, naz) =
		    3.0 / 8.0 * tau(nr, naz) +
		    std::sqrt(3.0) / 4.0 +
		    1.0 / (4.0 * tau(nr, naz) + parameters::tau_min);
	    }

	    if (parameters::opacity ==
		parameters::opacity_simple) { // Compare D'Angelo et. al
					      // 2003 eq.(28)
		tau_eff(nr, naz) = 3.0 / 8.0 * tau(nr, naz);
	    }

    }
    }
}

} // namespace compute