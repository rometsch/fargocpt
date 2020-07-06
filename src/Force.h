#ifndef FORCE_H
#define FORCE_H

#include "types.h"

Pair ComputeAccel(t_data &data, double x, double y);
double compute_smoothing(double r, t_data &data, const int n_radial,
			 const int n_azimuthal);
#endif // FORCE_H
