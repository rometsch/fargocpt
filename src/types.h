#ifndef TYPES_H
#define TYPES_H

/**
	\file types.h

	Definition of the structures used in the FARGO code.
*/

#include "polargrid.h"
#include <sys/times.h>

/**
	The boolean type will be used mainly for the variables corresponding to
   the command line switches
*/
typedef int boolean;

/**
	Set of two reals. It is used whenever a set of two reals is needed, and
   it usually represents a vector in the (x,y) plane (e.g., a %force), but not
   only: it is for instance used to store the mass inside and outside the orbit
   in MassInOut().
*/
struct pair {
    double x;
    double y;
};

typedef struct pair Pair;

class BoundaryFlow
{
  public:
    // Mass flow over inner and outer boundary and mass change by wave damping
    // boundaries and floor density
    double InnerBoundaryInflow;
    double InnerBoundaryOutflow;
    double OuterBoundaryInflow;
    double OuterBoundaryOutflow;
	double InnerWaveDampingMassCreation;
	double InnerWaveDampingMassRemoval;
	double OuterWaveDampingMassCreation;
	double OuterWaveDampingMassRemoval;
    double FloorMassCreation;

	BoundaryFlow() { reset(); }

    void reset()
    {
	InnerBoundaryInflow = 0.0;
	InnerBoundaryOutflow = 0.0;
	OuterBoundaryOutflow = 0.0;
	OuterBoundaryInflow = 0.0;
	InnerWaveDampingMassCreation = 0.0;
	InnerWaveDampingMassRemoval = 0.0;
	OuterWaveDampingMassCreation = 0.0;
	OuterWaveDampingMassRemoval = 0.0;
	FloorMassCreation = 0.0;
    }
};

#define YES true
#define NO false

#define EULER_INTEGRATOR 0
#define LEAPFROG_DRIFT_KICK_DRIFT 1
#define LEAPFROG_KICK_DRIFT_KICK 2

#define INDIRECT_TERM_REBOUND 0
#define INDIRECT_TERM_EULER 1

#define CONST_ALPHA 0
#define SCURVE_ALPHA 1
#define ALPHA_STAR_DIST_DEPENDEND 2

// TODO: This could be problematic on large NR
#define MAX1D 16384

#endif // TYPES_H
