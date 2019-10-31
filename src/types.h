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

struct torque_set {
    double InnerDisk;
    double OuterDisk;
    double ExcludeHill;
    double Total;
};

typedef struct torque_set TorqueSet;

struct triplet {
    double x;
    double y;
    double z;
};

typedef struct triplet Triplet;

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

/**
	An 8-component vector that contains some details on the force acting on
   a planet. The eight components arise from the subdivision of the disk into
   inner and outer disk (the %force arising from each part is evaluated
   separately) and from the fact that we contemplate the force arising from all
   the disk material and the force that excludes part of the Hill sphere.
*/
struct force {
    /// x-component of the force arising from the inner disk, without Hill
    /// sphere avoidance
    double fx_inner;
    /// y-component of the force arising from the inner disk, without Hill
    /// sphere avoidance
    double fy_inner;
    /// x-component of the force arising from the inner disk, with Hill
    /// sphere avoidance
    double fx_ex_inner;
    /// y-component of the force arising from the inner disk, with Hill
    /// sphere avoidance
    double fy_ex_inner;
    /// x-component of the force arising from the outer disk, without Hill
    /// sphere avoidance
    double fx_outer;
    /// y-component of the force arising from the outer disk, without Hill
    /// sphere avoidance
    double fy_outer;
    /// x-component of the force arising from the outer disk, with Hill
    /// sphere avoidance
    double fx_ex_outer;
    /// x-component of the force arising from the outer disk, with Hill
    /// sphere avoidance
    double fy_ex_outer;

    double *GlobalForce;
};

typedef struct force Force;

class BoundaryFlow
{
  public:
    // Mass flow over inner and outer boundary and mass change by wave
    // damping boundaries and floor density
    double InnerPositive;
    double InnerNegative;
    double OuterPositive;
    double OuterNegative;
    double WaveDampingPositive;
    double WaveDampingNegative;
    double FloorPositive;

    BoundaryFlow() { reset(); };

    void reset()
    {
	InnerPositive = 0.0;
	InnerNegative = 0.0;
	OuterPositive = 0.0;
	OuterNegative = 0.0;
	WaveDampingPositive = 0.0;
	WaveDampingNegative = 0.0;
	FloorPositive = 0.0;
    }
};

#define YES 1
#define NO 0

// TODO: This could be problematic on large NR
#define MAX1D 16384

#endif // TYPES_H
