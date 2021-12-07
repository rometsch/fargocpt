#include "util.h"
#include "global.h"
#include <cmath>
#include <cstring>
#include <stdint.h>
#include <string>

bool is_number(std::string s)
{
    for (std::uint32_t i = 0; i < s.length(); i++)
	if (isdigit(s[i]) == false)
	    return false;

    return true;
}

unsigned int get_next_azimuthal_id(const unsigned int id)
{
    unsigned int id_next = id + 1;
    if (id_next == NAzimuthal) {
	return 0;
    } else {
	return id_next;
    }
}

unsigned int get_prev_azimuthal_id(const unsigned int id)
{
    if (id == 0) {
	return NAzimuthal - 1;
    } else {
	return id - 1;
    }
}

unsigned int get_cell_id(const int nRadial, const int nAzimuthal)
{
    return nAzimuthal + (nRadial * NAzimuthal);
}

/** helper function to sum up a quantity inside the processes domain without
 * ghost cells */
void sum_without_ghost_cells(double &accumulator, const double &addend,
			     const unsigned int &n_radial)
{
    if (One_or_active <= n_radial && n_radial < Max_or_active) {
	accumulator += addend;
    }
}

/**
	cutoff function
		x < point-width -> 1
		x > point+width -> 0
		x > point-width & x < point+width -> smooth transition from 1 to
   0
*/
double cutoff_outer(double point, double width, double x)
{
    /*	if (x < point-width) {
		    return 1.0;
	    } else if (x > point+width) {
		    return 0.0;
	    } else {
		    return 0.5+0.5*cos((x-(point-width))*M_PI/2.0/width);
	    } */

    // return 1.0/(exp(8.0/width*(x-point))+1.0);
    return 1.0 / (1.0 + exp((x - point) / width));
}

/**
	cutoff function
		x > point-width -> 1
		x < point+width -> 0
		x < point-width & x < point+width -> smooth transition from 1 to
   0
*/
double cutoff_inner(double point, double width, double x)
{
    return 1.0 / (1.0 + exp((point - x) / width));
}

bool is_big_endian(void)
{
    union {
	uint32_t i;
	char c[4];
    } bint = {0x01020304};

    return bint.c[0] == 1;
}

bool is_distance_zero(double x)
{
    // in cgs units 1e-7 equals 1 nm
    // this is zero for astonomical scales
    return x * x < 1e-26;
}
