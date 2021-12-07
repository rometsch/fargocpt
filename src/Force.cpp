/**
	\file Force.cpp

	Contains the function used to evaluate the %force due to disk, and the
function that writes the 'tqwk' log files. Although the planet mass is given as
an argument to ComputeForce(), this mass is used only to specify the distance
cutoff in the case of the Hill sphere avoidance. The force returned is a
specific force. It has therefore the dimension of an acceleration (LT^-2).
**/

#include <cmath>
#include <float.h>
#include <stdio.h>

#include "Force.h"
#include "LowTasks.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "parameters.h"
#include "util.h"
#include "viscosity.h"

/**
	Computes the acceleration due to the disk on an object at position (x,y)
*/
Pair ComputeDiskOnNbodyAccel(t_data &data, double x, double y)
{
    Pair acceleration;
    double localaccel[4] = {0., 0., 0., 0.};
    double globalaccel[4] = {0., 0., 0., 0.};
    double axi, ayi, axo, ayo;

    const unsigned int ns = data[t_data::DENSITY].Nsec;
    const double *cell_center_x = CellCenterX->Field;
    const double *cell_center_y = CellCenterY->Field;
    axi = ayi = axo = ayo = 0.0;
    const double a = sqrt(x * x + y * y);

    for (unsigned int n_rad = One_or_active; n_rad < MaxMO_or_active; ++n_rad) {
	for (unsigned int n_az = 0; n_az < ns; ++n_az) {
	    // calculate smoothing length if dependend on radius
	    // i.e. for thickness smoothing with scale height at cell location
	    const double smooth =
		compute_smoothing(Rmed[n_rad], data, n_rad, n_az);
	    const int cell_id = n_az + n_rad * ns;
	    const double xc = cell_center_x[cell_id];
	    const double yc = cell_center_y[cell_id];
	    const double cellmass =
		Surf[n_rad] * data[t_data::DENSITY](n_rad, n_az);

	    const double dx = xc - x;
	    const double dy = yc - y;
	    const double dist_2 = std::pow(dx, 2) + std::pow(dy, 2);
	    const double dist_sm_2 = dist_2 + std::pow(smooth, 2);
	    const double dist_sm_3 = dist_sm_2 * std::sqrt(dist_2);
	    const double inv_dist_sm_3 = 1.0 / dist_sm_3;

	    if (Rmed[n_rad] < a) {
		axi += constants::G * cellmass * dx * inv_dist_sm_3;
		ayi += constants::G * cellmass * dy * inv_dist_sm_3;
	    } else {
		axo += constants::G * cellmass * dx * inv_dist_sm_3;
		ayo += constants::G * cellmass * dy * inv_dist_sm_3;
	    }
	}
    }

    localaccel[0] = axi;
    localaccel[1] = ayi;
    localaccel[2] = axo;
    localaccel[3] = ayo;
    MPI_Allreduce(&localaccel, &globalaccel, 4, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);

    acceleration.x = globalaccel[0] + globalaccel[2];
    acceleration.y = globalaccel[1] + globalaccel[3];

    return acceleration;
}

double compute_smoothing(double r, t_data &data, const int n_radial,
			 const int n_azimuthal)
{
    double smooth;
    const double scale_height =
	data[t_data::SCALE_HEIGHT](n_radial, n_azimuthal);
    smooth = parameters::thickness_smoothing * scale_height;
    return smooth;
}

///
/// \brief compute_smoothing_r
/// compute smoothing length at the bottom interface of the cell where v_radial
/// is defined
///
double compute_smoothing_r(t_data &data, const int n_radial,
			   const int n_azimuthal)
{
    double smooth;
    const double scale_height =
	0.5 *
	(data[t_data::SCALE_HEIGHT](n_radial, n_azimuthal) +
	 data[t_data::SCALE_HEIGHT](n_radial - 1, n_azimuthal));
    smooth = parameters::thickness_smoothing * scale_height;
    return smooth;
}

///
/// \brief compute_smoothing_az
/// compute smoothing length at the left interface of the cell where v_azimuthal
/// is defined
///
double compute_smoothing_az(t_data &data, const int n_radial,
			    const int n_azimuthal)
{
    int prev_n_az = n_azimuthal == 0
			? data[t_data::SCALE_HEIGHT].get_max_azimuthal()
			: n_azimuthal - 1;
    double smooth;
    const double scale_height =
	0.5 *
	(data[t_data::SCALE_HEIGHT](n_radial, n_azimuthal) +
	 data[t_data::SCALE_HEIGHT](n_radial, prev_n_az));
    smooth = parameters::thickness_smoothing * scale_height;
    return smooth;
}
