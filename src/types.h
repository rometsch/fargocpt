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
    double InnerPositive;
    double InnerNegative;
    double OuterPositive;
    double OuterNegative;
	double InnerWaveDampingPositive;
	double InnerWaveDampingNegative;
	double OuterWaveDampingPositive;
	double OuterWaveDampingNegative;
    double FloorPositive;

	BoundaryFlow() { reset(); }

    void reset()
    {
	InnerPositive = 0.0;
	InnerNegative = 0.0;
	OuterPositive = 0.0;
	OuterNegative = 0.0;
	InnerWaveDampingPositive = 0.0;
	InnerWaveDampingNegative = 0.0;
	OuterWaveDampingPositive = 0.0;
	OuterWaveDampingNegative = 0.0;
	FloorPositive = 0.0;
    }
};

#define YES true
#define NO false


#define INDIRECT_TERM_REBOUND 0
#define INDIRECT_TERM_EULER 1

// TODO: This could be problematic on large NR
#define MAX1D 16384

#endif // TYPES_H
