#pragma once

#include "types.h"

Pair ComputeDiskOnPlanetAccel(t_data &data, const double x, const double y,
			      const double l1);
double compute_smoothing(t_data &data, const int n_radial,
			 const int n_azimuthal);
double compute_smoothing_iso_planet(const double Rp);
double compute_smoothing_r(t_data &data, const int n_radial,
			   const int n_azimuthal);
double compute_smoothing_az(t_data &data, const int n_radial,
			    const int n_azimuthal);
