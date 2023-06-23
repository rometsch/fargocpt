#include "data.h"
#include "parameters.h"
#include "compute.h"
#include "SourceEuler.h"

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
	const unsigned int Nphi = rho.get_size_azimuthal();
	const double factor = parameters::density_factor;

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {
			rho(nr, naz) = Sig(nr, naz) / (factor * H(nr, naz));
		}
    }
}



}