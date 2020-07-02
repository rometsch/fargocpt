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


/**
	Computes the force due to the disk on an object at position (x,y)
*/
Pair ComputeAccel(t_data &data, double x, double y, double mass)
{
	Pair acceleration;
    int l, ns;
	double localaccel[4] = {0., 0., 0., 0.},
	   globalaccel[4] = {0., 0., 0., 0.};
	double xc, yc, cellmass, dx, dy, distance, dist2, a;
	double InvDist3, fxi, fyi, fxo, fyo;
    double rsmoothing = 0.0;

    ns = data[t_data::DENSITY].Nsec;
    const double* cell_center_x = CellCenterX->Field;
    const double* cell_center_y = CellCenterY->Field;
	fxi = fyi = fxo = fyo = 0.0;
    a = sqrt(x * x + y * y);

    bool SmoothingEnabled = (a != 0.0);

	for (unsigned int n_radial = One_or_active; n_radial < MaxMO_or_active;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
	     ++n_azimuthal) {
	    // calculate smoothing length if dependend on radius
	    // i.e. for thickness smoothing with scale height at cell location
	    if (SmoothingEnabled) {
		rsmoothing = compute_smoothing(Rmed[n_radial], data, n_radial,
					       n_azimuthal);
		}
	    l = n_azimuthal + n_radial * ns;
	    xc = cell_center_x[l];
	    yc = cell_center_y[l];
	    cellmass =
		Surf[n_radial] * data[t_data::DENSITY](n_radial, n_azimuthal);
	    dx = xc - x;
	    dy = yc - y;
	    dist2 = dx * dx + dy * dy;
	    dist2 += rsmoothing * rsmoothing;
	    distance = sqrt(dist2);
	    InvDist3 = 1.0 / dist2 / distance;
	    if (Rmed[n_radial] < a) {
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
