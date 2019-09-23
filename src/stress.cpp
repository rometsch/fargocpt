#include "stress.h"
#include "units.h"
#include "global.h"
#include "constants.h"
#include "selfgravity.h"
#include "output.h"

namespace stress {

void calculate_gravitational_stress(t_data &data)
{
	for (unsigned int n_radial = 0; n_radial <= data[t_data::T_GRAVITATIONAL].get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::T_GRAVITATIONAL].get_max_azimuthal(); ++n_azimuthal) {
			// TODO: Factor sqrt(2*PI) instead of 2 should be more accurate
			data[t_data::T_GRAVITATIONAL](n_radial, n_azimuthal) = 1.0/(4.0*PI*constants::G) * selfgravity::g_radial[n_radial*NAzimuthal+n_azimuthal] * selfgravity::g_azimuthal[n_radial*NAzimuthal+n_azimuthal] * (2*ASPECTRATIO_REF*Rmed[n_radial]);
		}
	}
}

void calculate_Reynolds_stress(t_data &data)
{
	// T_Reynolds = Sigma*(v_radial - v_radial_mean)*(v_azimuthal - v_azimuthal_mean)
	double v_radial_mean, v_azimuthal_mean;
	for (unsigned int n_radial = 0; n_radial <= data[t_data::T_REYNOLDS].get_max_radial(); ++n_radial) {
		// calculate v_radial_mean & v_azimuthal_mean
		v_radial_mean = 0.0;
		v_azimuthal_mean = 0.0;

		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::V_AZIMUTHAL].get_max_azimuthal(); ++n_azimuthal) {
			v_radial_mean += 0.5*(data[t_data::V_RADIAL](n_radial, n_azimuthal) + data[t_data::V_RADIAL](n_radial+1, n_azimuthal));
			v_azimuthal_mean += data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal);
		}
		v_azimuthal_mean /= (double)data[t_data::V_AZIMUTHAL].get_size_azimuthal();
		v_radial_mean /= (double)data[t_data::V_RADIAL].get_size_azimuthal();

		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::T_REYNOLDS].get_max_azimuthal(); ++n_azimuthal) {
			data[t_data::T_REYNOLDS](n_radial, n_azimuthal) = data[t_data::DENSITY](n_radial,n_azimuthal) * (0.5*(data[t_data::V_RADIAL](n_radial, n_azimuthal)+data[t_data::V_RADIAL](n_radial+1, n_azimuthal))-v_radial_mean) * (0.5*(data[t_data::V_AZIMUTHAL](n_radial,n_azimuthal)+data[t_data::V_AZIMUTHAL](n_radial,n_azimuthal+1 == data[t_data::V_AZIMUTHAL].get_max_azimuthal() ? 0 : n_azimuthal+1))-v_azimuthal_mean);
		}
	}
}

}
