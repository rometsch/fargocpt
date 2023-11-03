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
#include "viscosity/viscosity.h"

/**
	Computes the acceleration due to the disk on an object at position (x,y)
*/
Pair ComputeDiskOnPlanetAccel(t_data &data, const unsigned int nb)
{

	t_planetary_system & psys = data.get_planetary_system();
	t_planet & planet = psys.get_planet(nb);
	const double x = planet.get_x();
	const double y = planet.get_y();
	const double a = planet.get_r();


    Pair acceleration;
    double localaccel[4] = {0., 0., 0., 0.};
    double globalaccel[4] = {0., 0., 0., 0.};
    double axi, ayi, axo, ayo;

    const unsigned int ns = data[t_data::SIGMA].Nsec;
	const auto & sigma = data[t_data::SIGMA];
	const auto & sigma1d = data[t_data::SIGMA_1D];
    const double *cell_center_x = CellCenterX->Field;
    const double *cell_center_y = CellCenterY->Field;
    axi = ayi = axo = ayo = 0.0;

	#pragma omp parallel for collapse(2) reduction(+ : axi, ayi, axo, ayo)
    for (unsigned int n_rad = radial_first_active; n_rad < radial_active_size;
	 ++n_rad) {
	for (unsigned int n_az = 0; n_az < ns; ++n_az) {
	    // calculate smoothing length if dependend on radius
	    // i.e. for thickness smoothing with scale height at cell location
		double smooth = compute_smoothing(data, n_rad, n_az, nb);


		// Phi = GMm / r_sm
		// r_sm = sqrt(r**2 + (eps * H)**2)

	    const int cell_id = n_az + n_rad * ns;
	    const double xc = cell_center_x[cell_id];
	    const double yc = cell_center_y[cell_id];

		double cell_sigma = sigma(n_rad, n_az);
		if (parameters::correct_disk_selfgravity) {
			cell_sigma -= sigma1d(n_rad);
		}
	    const double cellmass = Surf[n_rad] * cell_sigma;

	    const double dx = xc - x;
	    const double dy = yc - y;
	    const double dist_2 = std::pow(dx, 2) + std::pow(dy, 2);
		const double dist_sm_2 = dist_2 + std::pow(smooth, 2);
		const double dist_sm = std::sqrt(dist_sm_2);
		const double dist_sm_3 = dist_sm_2 * dist_sm;
	    const double inv_dist_sm_3 = 1.0 / dist_sm_3;

	    // just to be consistent with the force the gas feels from the
	    // planets
	    double smooth_factor_klahr = 1.0;
	    if (parameters::ASPECTRATIO_MODE == 1) {
		/// scale height is reduced by the planets and can cause the
		/// epsilon smoothing be not sufficient for numerical stability.
		/// Thus we add the gravitational potential smoothing proposed
		/// by Klahr & Kley 2005; but the derivative of it, since we
		/// apply it directly on the force
		if (std::sqrt(x*x + y*y) > 10.e-10) {
			const double l1 = planet.get_dimensionless_roche_radius() *
				planet.get_distance_to_primary();
			const double r_sm = l1 * parameters::klahr_smoothing_radius;

			if (dist_sm < r_sm) {
			smooth_factor_klahr =
				-(3.0 * std::pow(dist_sm / r_sm, 4.0) -
				  4.0 * std::pow(dist_sm / r_sm, 3.0));
		    }
		}
		}

	    if (Rmed[n_rad] < a) {
		axi += constants::G * cellmass * dx * inv_dist_sm_3 *
			   smooth_factor_klahr;
		ayi += constants::G * cellmass * dy * inv_dist_sm_3 *
			   smooth_factor_klahr;
	    } else {
		axo += constants::G * cellmass * dx * inv_dist_sm_3 *
			   smooth_factor_klahr;
		ayo += constants::G * cellmass * dy * inv_dist_sm_3 *
			   smooth_factor_klahr;
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

inline static double compute_smoothing_scaleheight(t_data &data, const int n_radial,
			 const int n_azimuthal)
{
    const double scale_height =
	data[t_data::SCALE_HEIGHT](n_radial, n_azimuthal);
    const double smooth = parameters::thickness_smoothing * scale_height;
    return smooth;
}

inline static double compute_smoothing_iso_planet(t_data &data, const double nb)
{
	t_planetary_system & psys = data.get_planetary_system();
	t_planet & planet = psys.get_planet(nb);
	const double a = planet.get_r();
	const double h0 = parameters::ASPECTRATIO_REF;
	const double beta = parameters::FLARINGINDEX;
    const double scale_height = h0*std::pow(a, 1+beta);
    const double smooth = parameters::thickness_smoothing * scale_height;
    return smooth;
}

double compute_smoothing(t_data &data, const int n_radial,
			 const int n_azimuthal, const unsigned int nb)
{
	double rv;
	if (parameters::compatibility_no_star_smoothing && nb == 0) {
		return 0;
	}

	if (parameters::compatibility_smoothing_planetloc) {
		rv = compute_smoothing_iso_planet(data, nb);
	} else {
		rv = compute_smoothing_scaleheight(data, n_radial, n_azimuthal);
	}
	return rv;
}





inline static double compute_smoothing_r_scaleheight(t_data &data, const int n_radial,
			    const int n_azimuthal)
{
    const double scale_height =
	0.5 * (data[t_data::SCALE_HEIGHT](n_radial, n_azimuthal) +
	       data[t_data::SCALE_HEIGHT](n_radial - 1, n_azimuthal));
    const double smooth = parameters::thickness_smoothing * scale_height;
    return smooth;
}


///
/// \brief compute_smoothing_r
/// compute smoothing length at the bottom interface of the cell where v_radial
/// is defined
///
double compute_smoothing_r(t_data &data, const int n_radial,
			   const int n_azimuthal, const unsigned int nb)
{
	double rv;
	if (parameters::compatibility_no_star_smoothing && nb == 0) {
		return 0;
	}

	if (parameters::compatibility_smoothing_planetloc) {
		rv = compute_smoothing_iso_planet(data, nb);
	} else {
		rv = compute_smoothing_r_scaleheight(data, n_radial, n_azimuthal);
	}
	return rv;
}


inline static double compute_smoothing_az_scaleheight(t_data &data, const int n_radial,
			    const int n_azimuthal)
{
    int prev_n_az = n_azimuthal == 0
			? data[t_data::SCALE_HEIGHT].get_max_azimuthal()
			: n_azimuthal - 1;
    const double scale_height =
	0.5 * (data[t_data::SCALE_HEIGHT](n_radial, n_azimuthal) +
	       data[t_data::SCALE_HEIGHT](n_radial, prev_n_az));
    const double smooth = parameters::thickness_smoothing * scale_height;
    return smooth;
}

///
/// \brief compute_smoothing_az
/// compute smoothing length at the left interface of the cell where v_azimuthal
/// is defined
///
double compute_smoothing_az(t_data &data, const int n_radial,
			    const int n_azimuthal, const unsigned int nb)
{
	double rv;
	if (parameters::compatibility_no_star_smoothing && nb == 0) {
		return 0;
	}

	if (parameters::compatibility_smoothing_planetloc) {
		rv = compute_smoothing_iso_planet(data, nb);
	} else {
		rv = compute_smoothing_az_scaleheight(data, n_radial, n_azimuthal);
	}
	return rv;
}
