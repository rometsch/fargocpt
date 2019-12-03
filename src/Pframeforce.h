#ifndef PFRAMEFORCE_H
#define PFRAMEFORCE_H

#include "data.h"
#include "types.h"

void ComputeIndirectTerm(t_data &data);
void ComputeDiskOnNbodyAccel(Force *force, t_data &data);
void ComputeNbodyOnNbodyAccel(t_planetary_system &planetary_system);
void CalculatePotential(t_data &data);
void UpdatePlanetVelocitiesWithDiskForce(t_data &data, double dt);
void AdvanceSystemRK5(t_data &data, double dt);
double ConstructSequence(double *u, double *v, int n);
void AccreteOntoPlanets(t_data &data, double dt);

#endif // PFRAMEFORCE_H
