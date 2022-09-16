#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdlib.h>
#include <vector>

#include "Force.h"
#include "LowTasks.h"
#include "Pframeforce.h"
#include "RungeKutta.h"
#include "SideEuler.h"
#include "Theo.h"
#include "axilib.h"
#include "constants.h"
#include "find_cell_id.h"
#include "global.h"
#include "logging.h"
#include "parameters.h"
#include "selfgravity.h"
#include "units.h"
#include "util.h"
#include "viscosity/viscosity.h"
#include "quantities.h"
#include "frame_of_reference.h"

/* Below : work in non-rotating frame */
/**
 * @brief CalculatePotential: Nbody Potential caused by stars and planets
 * @param data
 */
void CalculateNbodyPotential(t_data &data, const double current_time)
{
    static const unsigned int N_planets =
	data.get_planetary_system().get_number_of_planets();

    // setup planet data
    for (unsigned int k = 0; k < N_planets; k++) {
	t_planet &planet = data.get_planetary_system().get_planet(k);
	g_mpl[k] = planet.get_rampup_mass(current_time);
	g_xpl[k] = planet.get_x();
	g_ypl[k] = planet.get_y();

	if (parameters::ASPECTRATIO_MODE == 1) {
		g_l1pl[k] = planet.get_dimensionless_roche_radius() *
		      planet.get_distance_to_primary();
	}
    }

    auto &pot = data[t_data::POTENTIAL];
    pot.clear();

    const unsigned int N_rad_max = pot.get_max_radial();
    const unsigned int N_az_max = pot.get_max_azimuthal();

	#pragma omp parallel for collapse(2)
    for (unsigned int n_rad = 0; n_rad <= N_rad_max; ++n_rad) {
	for (unsigned int n_az = 0; n_az <= N_az_max; ++n_az) {
	    const int cell = get_cell_id(n_rad, n_az);
	    const double x = CellCenterX->Field[cell];
	    const double y = CellCenterY->Field[cell];

	    for (unsigned int k = 0; k < N_planets; k++) {

		const double smooth = compute_smoothing(data, n_rad, n_az);
		const double dx = x - g_xpl[k];
		const double dy = y - g_ypl[k];
		const double dist_2 = std::pow(dx, 2) + std::pow(dy, 2);
		const double d_smoothed =
		    std::sqrt(dist_2 + std::pow(smooth, 2));

		double smooth_factor_klahr = 1.0;

		if (parameters::ASPECTRATIO_MODE == 1) {
		    /// scale height is reduced by the planets and can cause the
		    /// epsilon smoothing be not sufficient for numerical
		    /// stability. Thus we add the gravitational potential
		    /// smoothing proposed by Klahr & Kley 2005.
			if (std::sqrt(g_xpl[k]*g_xpl[k] + g_ypl[k]*g_ypl[k]) > 1.0e-10) { // only for non central objects
			// position of the l1 point between planet and central star.
			const double l1 = g_l1pl[k];
			const double r_sm =
			    l1 * parameters::klahr_smoothing_radius;

			const double dist = std::sqrt(dist_2);

			if (dist < r_sm) {
			    smooth_factor_klahr =
				(std::pow(dist / r_sm, 4.0) -
				 2.0 * std::pow(dist / r_sm, 3.0) +
				 2.0 * dist / r_sm);
			}
		    }
		}

		// direct term from planet
		pot(n_rad, n_az) +=
			-constants::G * g_mpl[k] / d_smoothed * smooth_factor_klahr;
	    }
	    // apply indirect term
	    // correct frame with contributions from disk and planets
	    pot(n_rad, n_az) += -refframe::IndirectTerm.x * x - refframe::IndirectTerm.y * y;
	}
    }
}

void CalculateAccelOnGas(t_data &data, const double current_time)
{

    static const unsigned int N_planets =
	data.get_planetary_system().get_number_of_planets();

    // setup planet data
    for (unsigned int k = 0; k < N_planets; k++) {
	t_planet &planet = data.get_planetary_system().get_planet(k);
	g_mpl[k] = planet.get_rampup_mass(current_time);
	g_xpl[k] = planet.get_x();
	g_ypl[k] = planet.get_y();

	if (parameters::ASPECTRATIO_MODE == 1) {
		g_l1pl[k] = planet.get_dimensionless_roche_radius() *
		      planet.get_distance_to_primary();
	}
    }

    double *acc_r = data[t_data::ACCEL_RADIAL].Field;
    const unsigned int N_az_max =
	data[t_data::ACCEL_RADIAL].get_size_azimuthal();
	const unsigned int Nr = data[t_data::ACCEL_RADIAL].get_size_radial() - 1;

	#pragma omp parallel for collapse(2)
    for (unsigned int n_rad = 1;
	 n_rad < Nr; ++n_rad) { // No need to compute Vr at the top of the outer ghost cells
	for (unsigned int n_az = 0; n_az < N_az_max; ++n_az) {

	    const double phi = (double)n_az * dphi;
	    const double r = Rinf[n_rad];

	    const double x = r * std::cos(phi);
	    const double y = r * std::sin(phi);
	    const double smooth = compute_smoothing_r(data, n_rad, n_az);

	    pair ar = refframe::IndirectTerm;
	    for (unsigned int k = 0; k < N_planets; k++) {

		const double dx = x - g_xpl[k];
		const double dy = y - g_ypl[k];
		const double dist_2 = std::pow(dx, 2) + std::pow(dy, 2);
		const double dist_2_sm = dist_2 + std::pow(smooth, 2);
		const double dist_sm = std::sqrt(dist_2_sm);
		const double dist_3_sm = dist_sm * dist_2_sm;
		const double inv_dist_3_sm = 1.0 / dist_3_sm;

		double smooth_factor_klahr = 1.0;

		if (parameters::ASPECTRATIO_MODE == 1) {
		    /// scale height is reduced by the planets and can cause the
		    /// epsilon smoothing be not sufficient for numerical
		    /// stability. Thus we add the gravitational potential
		    /// smoothing proposed by Klahr & Kley 2005; but the
		    /// derivative of it, since we apply it directly on the
		    /// force
			if (std::sqrt(g_xpl[k]*g_xpl[k] + g_ypl[k]*g_ypl[k]) > 1.0e-10) { // only for non central objects
			// position of the l1 point between planet and central
			// star.
			const double l1 = g_l1pl[k];
			const double r_sm =
			    l1 * parameters::klahr_smoothing_radius;

			const double dist = std::sqrt(dist_2);

			if (dist < r_sm) {
			    smooth_factor_klahr =
				-(3.0 * std::pow(dist / r_sm, 4.0) -
				  4.0 * std::pow(dist / r_sm, 3.0));
			}
		    }
		}

		// direct term from planet
		ar.x -= dx * constants::G * g_mpl[k] * inv_dist_3_sm *
			smooth_factor_klahr;
		ar.y -= dy * constants::G * g_mpl[k] * inv_dist_3_sm *
			smooth_factor_klahr;
	    }

	    const double accel = std::cos(phi) * ar.x + std::sin(phi) * ar.y;

	    /*
	    /// DEBUG, by comparing with accel from potential
	    const double gradphi = (data[t_data::POTENTIAL](n_rad, n_az) -
		       data[t_data::POTENTIAL](n_rad - 1, n_az)) *
		      InvDiffRmed[n_rad];

	    if(accel != 0.0){
	    const double rel_err = std::fabs((gradphi + accel)/accel);

	    const double dx = x - xpl[1];
	    const double dy = y - ypl[1];
	    const double dist = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));

	    double diff_phi = phi -
	    data.get_planetary_system().get_planet(1).get_phi(); const double
	    diff_r = r - data.get_planetary_system().get_planet(1).get_r();
	    diff_phi = std::min(diff_phi, std::fabs(diff_phi - 2 * M_PI));

	    if(rel_err > 1.e-1){
	    printf("(%d %d)	accel r	= (%.5e	%.5e)	%.5e	(%.5e
	    %.5e)\n", n_rad, n_az, gradphi, -accel, rel_err, diff_r, diff_phi);
	    }}
	    /// END DEBUG
	    */

	    const int cell_id = get_cell_id(n_rad, n_az);
	    acc_r[cell_id] = accel;
	}
    }

    double *accel_phi = data[t_data::ACCEL_AZIMUTHAL].Field;

	#pragma omp parallel for collapse(2)
    for (unsigned int n_rad = 0;
	 n_rad < Nr; ++n_rad) {
	for (unsigned int n_az = 0; n_az < N_az_max; ++n_az) {

	    const double phi = ((double)n_az - 0.5) * dphi;
	    const double r = Rmed[n_rad];

	    const double x = r * std::cos(phi);
	    const double y = r * std::sin(phi);
	    const double smooth = compute_smoothing_az(data, n_rad, n_az);

	    pair aphi = refframe::IndirectTerm;
	    for (unsigned int k = 0; k < N_planets; k++) {

		const double dx = x - g_xpl[k];
		const double dy = y - g_ypl[k];
		const double dist_2 = std::pow(dx, 2) + std::pow(dy, 2);
		const double dist_2_sm = dist_2 + std::pow(smooth, 2);
		const double dist_3_sm = std::sqrt(dist_2_sm) * dist_2_sm;
		const double inv_dist_3_sm = 1.0 / dist_3_sm;

		double smooth_factor_klahr = 1.0;

		if (parameters::ASPECTRATIO_MODE == 1) {
		    /// scale height is reduced by the planets and can cause the
		    /// epsilon smoothing be not sufficient for numerical
		    /// stability. Thus we add the gravitational potential
		    /// smoothing proposed by Klahr & Kley 2005; but the
		    /// derivative of it, since we apply it directly on the
		    /// force
			if (std::sqrt(g_xpl[k]*g_xpl[k] + g_ypl[k]*g_ypl[k]) > 1.0e-10) { // only for non central objects
			// position of the l1 point between planet and central
			// star.
			const double l1 = g_l1pl[k];
			const double r_sm =
			    l1 * parameters::klahr_smoothing_radius;

			const double dist = std::sqrt(dist_2);

			if (dist < r_sm) {
			    smooth_factor_klahr =
				-(3.0 * std::pow(dist / r_sm, 4.0) -
				  4.0 * std::pow(dist / r_sm, 3.0));
			}
		    }
		}

		// direct term from planet
		aphi.x -= dx * constants::G * g_mpl[k] * inv_dist_3_sm *
			  smooth_factor_klahr;
		aphi.y -= dy * constants::G * g_mpl[k] * inv_dist_3_sm *
			  smooth_factor_klahr;
	    }

	    const double accel =
		std::cos(phi) * aphi.y - std::sin(phi) * aphi.x;

	    /*
	    /// DEBUG, compare with accel from potential
	    const double invdxtheta = 1.0 / (dphi * Rmed[n_rad]);

	    // 1/r dPhi/dphi
	    const double gradphi =
	    (data[t_data::POTENTIAL](n_rad, n_az) -
	     data[t_data::POTENTIAL](
		     n_rad, n_az == 0
			       ? data[t_data::POTENTIAL].get_max_azimuthal()
			       : n_az - 1)) *
	    invdxtheta;

	    if(accel != 0.0){
	    const double rel_err = std::fabs((gradphi + accel)/accel);

	    const double dx = x - xpl[1];
	    const double dy = y - ypl[1];
	    const double dist = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));

	    double diff_phi = phi -
	    data.get_planetary_system().get_planet(1).get_phi(); diff_phi =
	    std::min(diff_phi, std::fabs(diff_phi - 2 * M_PI)); const double
	    diff_r = r - data.get_planetary_system().get_planet(1).get_r();
	    if(rel_err < 1e-5){
	    printf("(%d %d )	accel p	= (%.5e	%.5e)	%.5e	(%.5e %.5e)\n",
	    n_rad, n_az, gradphi, -accel, rel_err, diff_r, diff_phi);
	    }}
	    /// END DEBUG
	    */

	    const int cell_id = get_cell_id(n_rad, n_az);
	    accel_phi[cell_id] = accel;
	}
    }
}

/**
   Update disk forces onto planets if *diskfeedback* is turned on
*/
void ComputeDiskOnNbodyAccel(t_data &data)
{
    Pair accel;
    for (unsigned int k = 0;
	 k < data.get_planetary_system().get_number_of_planets(); k++) {
	t_planet &planet = data.get_planetary_system().get_planet(k);
	const double l1 = planet.get_dimensionless_roche_radius() *
		 planet.get_distance_to_primary();


	accel =
	    ComputeDiskOnPlanetAccel(data, planet.get_x(), planet.get_y(), l1);
	planet.set_disk_on_planet_acceleration(accel);

	const double torque =
	    (planet.get_x() * accel.y - planet.get_y() * accel.x) *
	    planet.get_mass();
	planet.set_torque(torque);
    }
}

/**
   Update mutual planet-planet accelerations
*/
void ComputeNbodyOnNbodyAccel(t_planetary_system &planetary_system)
{

    for (unsigned int npl = 0; npl < planetary_system.get_number_of_planets();
	 npl++) {
	t_planet &planet = planetary_system.get_planet(npl);
	const double x = planet.get_x();
	const double y = planet.get_y();
	double ax = 0.0;
	double ay = 0.0;
	for (unsigned int nother = 0;
	     nother < planetary_system.get_number_of_planets(); nother++) {
	    if (nother != npl) {
		t_planet &other_planet = planetary_system.get_planet(nother);
		const double xo = other_planet.get_x();
		const double yo = other_planet.get_y();
		const double mass = other_planet.get_mass();
		const double dist =
		    sqrt(std::pow(x - xo, 2) + std::pow(y - yo, 2));
		ax -= constants::G * mass / std::pow(dist, 3) * (x - xo);
		ay -= constants::G * mass / std::pow(dist, 3) * (y - yo);
	    }
	}
	planet.set_nbody_on_planet_acceleration_x(ax);
	planet.set_nbody_on_planet_acceleration_y(ay);
    }
}

static double q0[MAX1D], q1[MAX1D], PlanetMasses[MAX1D];
void ComputeNbodyOnNbodyAccelRK5(t_data &data, const double dt)
{
    unsigned int n = data.get_planetary_system().get_number_of_planets();

	if (parameters::integrate_planets) {
		for (unsigned int i = 0; i < data.get_planetary_system().get_number_of_planets(); i++) {
			q0[i+0*n] = data.get_planetary_system().get_planet(i).get_x();
			q0[i+1*n] = data.get_planetary_system().get_planet(i).get_y();
			q0[i+2*n] = data.get_planetary_system().get_planet(i).get_vx();
			q0[i+3*n] = data.get_planetary_system().get_planet(i).get_vy();
			PlanetMasses[i] = data.get_planetary_system().get_planet(i).get_mass();
		}

	RungeKutta(q0, dt, PlanetMasses, q1, n);

		for (unsigned int i = 0; i < data.get_planetary_system().get_number_of_planets(); i++) {
			data.get_planetary_system().get_planet(i).set_nbody_on_planet_acceleration_x(q1[i+2*n]);
			data.get_planetary_system().get_planet(i).set_nbody_on_planet_acceleration_y(q1[i+3*n]);
		}

	}
}

/**
	Updates planets velocities due to disk influence if "DiskFeedback" is
   set.
*/
void UpdatePlanetVelocitiesWithDiskForce(t_data &data, const double dt)
{
    for (unsigned int k = 0;
	 k < data.get_planetary_system().get_number_of_planets(); k++) {
	if (parameters::disk_feedback) {
	    t_planet &planet = data.get_planetary_system().get_planet(k);

	    const Pair gamma = planet.get_disk_on_planet_acceleration();
	    const double new_vx =
		planet.get_vx() + dt * gamma.x;
	    const double new_vy =
		planet.get_vy() + dt * gamma.y;
	    planet.set_vx(new_vx);
		planet.set_vy(new_vy);
	}
    }
}

double ConstructSequence(double *u, double *v, int n)
{
    int i;
    double lapl = 0.0;
    for (i = 1; i < n; i++) {
	u[i] = 2.0 * v[i] - u[i - 1];
    }
    for (i = 1; i < n - 1; i++) {
	lapl += fabs(u[i + 1] + u[i - 1] - 2.0 * u[i]);
    }
    return lapl;
}
