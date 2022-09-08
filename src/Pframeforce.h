#ifndef PFRAMEFORCE_H
#define PFRAMEFORCE_H

#include "data.h"

void ComputeDiskOnNbodyAccel(t_data &data);
void ComputeNbodyOnNbodyAccel(t_planetary_system &planetary_system);
void ComputeNbodyOnNbodyAccelRK5(t_data &data, const double dt);
void CalculateNbodyPotential(t_data &data);
void CalculateAccelOnGas(t_data &data);
void UpdatePlanetVelocitiesWithDiskForce(t_data &data, const double dt);
double ConstructSequence(double *u, double *v, int n);

#endif // PFRAMEFORCE_H
