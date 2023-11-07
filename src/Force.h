#pragma once

#include "types.h"
#include "data.h"

Pair ComputeDiskOnPlanetAccel(t_data &data, const unsigned nb);
double compute_smoothing(t_data &data, const int n_radial,
			 const int n_azimuthal, const unsigned nb);
double compute_smoothing_iso_planet(const double Rp);

