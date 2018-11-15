#ifndef THEO_H
#define THEO_H

#include "constants.h"
#include "types.h"
#include <math.h>

/**
	Calculates Kepler angular velocity

	\param r radius
	\returns kepler angular velocity
**/
inline double omega_kepler(double r)
{
	// omega_kepler = sqrt(G*M/r^3)
	// as G=1 and M=1: omega_kepler = r^-1.5
	return pow(r,-1.5);
}

void RefillSigma(t_polargrid* Density);
void RefillEnergy(t_polargrid* Energy);

#endif // THEO_H
