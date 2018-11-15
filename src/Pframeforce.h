#ifndef PFRAMEFORCE_H
#define PFRAMEFORCE_H

#include "types.h"
#include "data.h"

Pair ComputeIndirectTerm();
void FillForcesArrays(t_data &data);
void AdvanceSystemFromDisk(Force* force, t_data &data, double dt);
void AdvanceSystemRK5(t_data &data, double dt);
void SolveOrbits(t_data &data);
double ConstructSequence(double* u, double* v, int n);
void AccreteOntoPlanets(t_data &data, double dt);
void FindOrbitalElements (double x, double y, double vx, double vy, double m, int n);

#endif // PFRAMEFORCE_H
