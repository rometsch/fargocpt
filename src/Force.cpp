/**
	\file Force.cpp

	Contains the function used to evaluate the %force due to disk, and the
function that writes the 'tqwk' log files. Although the planet mass is given as
an argument to ComputeForce(), this mass is used only to specify the distance
cutoff in the case of the Hill sphere avoidance. The force returned is a
specific force. It has therefore the dimension of an acceleration (LT^-2).
**/

#include <float.h>
#include <math.h>
#include <stdio.h>

#include "Force.h"
#include "LowTasks.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "parameters.h"
#include "viscosity.h"
#include "util.h"


/**
	Computes the force due to the disk on an object at position (x,y)
*/
Pair ComputeAccel(t_data &data, double x, double y, double mass)
{
	Pair acceleration;
	double localaccel[4] = {0., 0., 0., 0.},
	   globalaccel[4] = {0., 0., 0., 0.};
    double r_sm = 0.0;

    const unsigned int N_sec = data[t_data::DENSITY].Nsec;
    const double* cell_center_x = CellCenterX->Field;
    const double* cell_center_y = CellCenterY->Field;
	double fxi, fyi, fxo, fyo;
	fxi = fyi = fxo = fyo = 0.0;
    const double a = sqrt(pow2(x) + pow2(y));

    bool SmoothingEnabled = (a != 0.0);

	const unsigned int N_az_max = data[t_data::POTENTIAL].get_max_azimuthal();

	for (unsigned int n_rad = One_or_active; n_rad < MaxMO_or_active; ++n_rad) {
	for (unsigned int n_az = 0; n_az <= N_az_max; ++n_az) {
	    // calculate smoothing length if dependend on radius
	    // i.e. for thickness smoothing with scale height at cell location
	    if (SmoothingEnabled) {
		r_sm = compute_smoothing(Rmed[n_rad], data, n_rad, n_az);
		}
	    const unsigned int l = n_az + n_rad * N_sec;
	    const double xc = cell_center_x[l];
	    const double yc = cell_center_y[l];
	    const double cellmass =
		Surf[n_rad] * data[t_data::DENSITY](n_rad, n_az);
	    const double dx = xc - x;
	    const double dy = yc - y;
	    const double dist2 = pow2(dx) + pow2(dy) + pow2(r_sm);
	    const double distance = sqrt(dist2);
	    const double InvDist3 = 1.0 / dist2 / distance;
	    if (Rmed[n_rad] < a) {
		fxi += constants::G * cellmass * dx * InvDist3;
		fyi += constants::G * cellmass * dy * InvDist3;
	    } else {
		fxo += constants::G * cellmass * dx * InvDist3;
		fyo += constants::G * cellmass * dy * InvDist3;
	    }
	}
    }

	localaccel[0] = fxi;
	localaccel[1] = fyi;
	localaccel[2] = fxo;
	localaccel[3] = fyo;
	MPI_Allreduce(&localaccel, &globalaccel, 4, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);

	acceleration.x = globalaccel[0] + globalaccel[2];
	acceleration.y = globalaccel[1] + globalaccel[3];

	return acceleration;
}

double compute_smoothing_isothermal(double r)
{
    double smooth;
    const double scale_height =
	ASPECTRATIO_REF * pow(r, 1.0 + FLARINGINDEX); // = H
    smooth = parameters::thickness_smoothing * scale_height;
    return smooth;
}

double compute_smoothing(double r, t_data &data, const int n_radial,
			 const int n_azimuthal)
{
    double smooth;
    const double scale_height =
	data[t_data::ASPECTRATIO](n_radial, n_azimuthal) * r;
    smooth = parameters::thickness_smoothing * scale_height;
    return smooth;
}
