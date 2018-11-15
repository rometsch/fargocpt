#ifndef FORCE_H
#define FORCE_H

#include "types.h"

Force *AllocateForce(int dimfxy);
void FreeForce(Force* force);
void ComputeForce(t_data &data, Force* force, double x, double y, double rsmoothing, double mass);
double compute_smoothing(double r);
void UpdateLog(t_data &data, Force* fc, int outputnb, double time);
#endif // FORCE_H
