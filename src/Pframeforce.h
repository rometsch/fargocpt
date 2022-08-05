#ifndef PFRAMEFORCE_H
#define PFRAMEFORCE_H

#include "data.h"

void ComputeIndirectTermDisk(t_data &data);
void ComputeIndirectTermNbody(t_data &data);
void ComputeIndirectTermNbodyAndFixVelocities(t_data &data, const double dt);
void compute_minimum_timestep_size(t_data &data);
void ComputeDiskOnNbodyAccel(t_data &data);
void ComputeNbodyOnNbodyAccel(t_planetary_system &planetary_system);
void ComputeNbodyOnNbodyAccelRK5(t_data &data, double dt);
void ComputeNbodyOnNbodyAccelRebound(t_planetary_system &planetary_system);
void CalculateNbodyPotential(t_data &data);
void CalculateAccelOnGas(t_data &data);
void UpdatePlanetVelocitiesWithDiskForce(t_data &data, double dt);
void AdvanceSystemRK5(t_data &data, double dt);
double ConstructSequence(double *u, double *v, int n);

#endif // PFRAMEFORCE_H
