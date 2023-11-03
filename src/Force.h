#pragma once

#include "types.h"

Pair ComputeDiskOnPlanetAccel(t_data &data, const unsigned nb);
double compute_smoothing(t_data &data, const int n_radial,
			 const int n_azimuthal, const unsigned int nb);
double compute_smoothing_r(t_data &data, const int n_radial,
			    const int n_azimuthal, const unsigned int nb);
double compute_smoothing_az(t_data &data, const int n_radial,
			    const int n_azimuthal, const unsigned int nb);
