#include "util.h"
#include <math.h>

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

