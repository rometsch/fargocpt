#ifndef FORCE_H
#define FORCE_H

#include "types.h"

Force *AllocateForce(int dimfxy);
void FreeForce(Force *force);
void ComputeForce(t_data &data, Force *force, double x, double y, double mass);
double compute_smoothing_isothermal(double r);
double compute_smoothing(double r, t_data &data, const int n_radial,
			 const int n_azimuthal);
void UpdateLog(t_data &data, Force *fc, int outputnb, double time);
#endif // FORCE_H
