#ifndef PFRAMEFORCE_H
#define PFRAMEFORCE_H

#include "types.h"
#include "data.h"

void ComputeIndirectTerm(Force* force,t_data &data);
void ComputeDiskOnNbodyAccel(Force* force,t_data &data);
void ComputeNbodyOnNbodyAccel(t_planetary_system &planetary_system);
void CalculatePotential(t_data &data);
void AdvanceSystemFromDisk(Force* force, t_data &data, double dt);
void AdvanceSystemRK5(t_data &data, double dt);
double ConstructSequence(double* u, double* v, int n);
void AccreteOntoPlanets(t_data &data, double dt);

#endif // PFRAMEFORCE_H
