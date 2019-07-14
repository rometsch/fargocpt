#include "util.h"
#include <math.h>
#include <stdint.h>

/**
	cutoff function
		x < point-width -> 1
		x > point+width -> 0
		x > point-width & x < point+width -> smooth transition from 1 to 0
*/
double cutoff(double point, double width, double x)
{
/*	if (x < point-width) {
		return 1.0;
	} else if (x > point+width) {
		return 0.0;
	} else {
		return 0.5+0.5*cos((x-(point-width))*M_PI/2.0/width);
	} */

	//return 1.0/(exp(8.0/width*(x-point))+1.0);
	return 1.0/(1.0+exp((x-point)/width));
}

bool is_big_endian(void)
{
    union {
        uint32_t i;
         char c[4];
    } bint = {0x01020304};

    return bint.c[0] == 1;
}

bool is_distance_zero(double x) {
	// in cgs units 1e-7 equals 1 nm
	// this is zero for astonomical scales
	return abs(x) < 1e-13;
}
