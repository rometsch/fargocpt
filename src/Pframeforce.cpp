#include <math.h>
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
#include "global.h"
#include "logging.h"
#include "parameters.h"
#include "selfgravity.h"
#include "units.h"
#include "util.h"
#include "viscosity.h"
#include "find_cell_id.h"

extern Pair IndirectTerm;
extern Pair IndirectTermDisk;
extern Pair IndirectTermPlanets;

/**
 * @brief ComputeIndirectTerm: IndirectTerm is the correction therm that needs
 * to be added to the accelerations.
 * @param force
 * @param data
 */
void ComputeIndirectTerm(t_data &data)
{
    IndirectTerm.x = 0.0;
    IndirectTerm.y = 0.0;
    IndirectTermDisk.x = 0.0;
    IndirectTermDisk.y = 0.0;
    IndirectTermPlanets.x = 0.0;
    IndirectTermPlanets.y = 0.0;

    // compute disk indirect term
    if (parameters::disk_feedback) {
	// add up contributions from disk on all bodies used to calculate the
	// center
	double mass_center = 0.0;
	for (unsigned int n = 0; n < parameters::n_bodies_for_hydroframe_center;
	     n++) {
	    t_planet &planet = data.get_planetary_system().get_planet(n);
	    const double mass = planet.get_mass();
	    const Pair accel = planet.get_disk_on_planet_acceleration();
	    IndirectTermDisk.x -= mass * accel.x;
	    IndirectTermDisk.y -= mass * accel.y;
	    mass_center += mass;
	}
	IndirectTermDisk.x /= mass_center;
	IndirectTermDisk.y /= mass_center;
    }

    // compute nbody indirect term
    // add up contributions from mutual interactions from all bodies used to
    // calculate the center
    double mass_center = 0.0;
    for (unsigned int n = 0; n < parameters::n_bodies_for_hydroframe_center;
	 n++) {
	t_planet &planet = data.get_planetary_system().get_planet(n);
	const double mass = planet.get_mass();
	const Pair accel = planet.get_nbody_on_planet_acceleration();
	IndirectTermPlanets.x -= mass * accel.x;
	IndirectTermPlanets.y -= mass * accel.y;
	mass_center += mass;
    }
    IndirectTermPlanets.x /= mass_center;
    IndirectTermPlanets.y /= mass_center;

    IndirectTerm.x = IndirectTermDisk.x + IndirectTermPlanets.x;
    IndirectTerm.y = IndirectTermDisk.y + IndirectTermPlanets.y;
}

/* Below : work in non-rotating frame */
/**
 * @brief CalculatePotential: Nbody Potential caused by stars and planets
 * @param data
 */
void CalculateNbodyPotential(t_data &data)
{
	const unsigned int N_planets = 
		data.get_planetary_system().get_number_of_planets();
    std::vector<double> xpl(N_planets);
    std::vector<double> ypl(N_planets);
    std::vector<double> mpl(N_planets);
    std::vector<double> smooth_pl(N_planets);

    // setup planet data
    for (unsigned int k = 0; k < N_planets; k++) {
	t_planet &planet = data.get_planetary_system().get_planet(k);
	mpl[k] = planet.get_rampup_mass();
	xpl[k] = planet.get_x();
	ypl[k] = planet.get_y();
    }

	auto& pot = data[t_data::POTENTIAL];
    pot.clear();
	
	const unsigned int N_rad_max = pot.get_max_radial();
	const unsigned int N_az_max = pot.get_max_azimuthal();

    for (unsigned int n_rad = 0; n_rad <= N_rad_max; ++n_rad) {
	for (unsigned int n_az = 0; n_az <= N_az_max; ++n_az) {
	    const double angle = (double)n_az /
		    (double)pot.get_size_azimuthal() * 2.0 * M_PI;
		const double x = Rmed[n_rad] * std::cos(angle);
		const double y = Rmed[n_rad] * std::sin(angle);

	    for (unsigned int k = 0; k < N_planets; k++) {

		const double smooth = compute_smoothing(Rmed[n_rad], data, n_rad, n_az);
		const double dx = x - xpl[k];
		const double dy = y - ypl[k];
		const double d_smoothed = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2) + std::pow(smooth, 2));

		// direct term from planet
		pot(n_rad, n_az) += -constants::G * mpl[k] / d_smoothed;
	    }
	    // apply indirect term
	    // correct frame with contributions from disk and planets
	    pot(n_rad, n_az) += -IndirectTerm.x * x - IndirectTerm.y * y;
	}
    }
}


void CalculateAccelOnGas(t_data &data)
{
	const int ns = data[t_data::DENSITY].Nsec;

	const unsigned int N_planets =
		data.get_planetary_system().get_number_of_planets();
	std::vector<double> xpl(N_planets);
	std::vector<double> ypl(N_planets);
	std::vector<double> mpl(N_planets);
	std::vector<double> smooth_pl(N_planets);

	// setup planet data
	for (unsigned int k = 0; k < N_planets; k++) {
	t_planet &planet = data.get_planetary_system().get_planet(k);
	mpl[k] = planet.get_rampup_mass();
	xpl[k] = planet.get_x();
	ypl[k] = planet.get_y();
	}

	double* acc_r = data[t_data::ACCEL_RADIAL].Field;
	const unsigned int N_az_max = data[t_data::ACCEL_RADIAL].get_size_azimuthal();

	for (unsigned int n_rad = 1; n_rad < data[t_data::ACCEL_RADIAL].get_size_radial()-1; ++n_rad) {
	for (unsigned int n_az = 0; n_az < N_az_max; ++n_az) {

		const double phi = (double)n_az * dphi;
		const double r = Rinf[n_rad];

		const double x = r * std::cos(phi);
		const double y = r * std::sin(phi);
		const double smooth = compute_smoothing_r(data, n_rad, n_az);

		pair ar = IndirectTerm;
		for (unsigned int k = 0; k < N_planets; k++) {

		const double dx = x - xpl[k];
		const double dy = y - ypl[k];
		const double dist_2 = std::pow(dx, 2) + std::pow(dy, 2);
		const double dist_2_sm = dist_2 + std::pow(smooth, 2);
		const double dist_3_sm = std::sqrt(dist_2) * dist_2_sm;
		const double inv_dist_3_sm = 1.0 / dist_3_sm;

		// direct term from planet
		ar.x -= dx * constants::G * mpl[k] * inv_dist_3_sm;
		ar.y -= dy * constants::G * mpl[k] * inv_dist_3_sm;
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

		double diff_phi = phi - data.get_planetary_system().get_planet(1).get_phi();
		const double diff_r = r - data.get_planetary_system().get_planet(1).get_r();
		diff_phi = std::min(diff_phi, std::fabs(diff_phi - 2 * M_PI));

		if(rel_err > 1.e-1){
		printf("(%d %d)	accel r	= (%.5e	%.5e)	%.5e	(%.5e	%.5e)\n", n_rad, n_az, gradphi, -accel, rel_err, diff_r, diff_phi);
		}}
		/// END DEBUG
		*/

		const int cell_id = n_az + n_rad * ns;
		acc_r[cell_id] = accel;
	}
	}

	double* accel_phi = data[t_data::ACCEL_AZIMUTHAL].Field;
	for (unsigned int n_rad = 0; n_rad < data[t_data::ACCEL_AZIMUTHAL].get_size_radial(); ++n_rad) {
	for (unsigned int n_az = 0; n_az < N_az_max; ++n_az) {

		const double phi = ((double)n_az - 0.5) * dphi;
		const double r = Rmed[n_rad];

		const double x = r * std::cos(phi);
		const double y = r * std::sin(phi);
		const double smooth = compute_smoothing_az(data, n_rad, n_az);

		pair aphi = IndirectTerm;
		for (unsigned int k = 0; k < N_planets; k++) {

		const double dx = x - xpl[k];
		const double dy = y - ypl[k];
		const double dist_2 = std::pow(dx, 2) + std::pow(dy, 2);
		const double dist_2_sm = dist_2 + std::pow(smooth, 2);
		const double dist_3_sm = std::sqrt(dist_2) * dist_2_sm;
		const double inv_dist_3_sm = 1.0 / dist_3_sm;

		// direct term from planet
		aphi.x -= dx * constants::G * mpl[k] * inv_dist_3_sm;
		aphi.y -= dy * constants::G * mpl[k] * inv_dist_3_sm;
		}

		const double accel = std::cos(phi) * aphi.y - std::sin(phi) * aphi.x;


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

		double diff_phi = phi - data.get_planetary_system().get_planet(1).get_phi();
		diff_phi = std::min(diff_phi, std::fabs(diff_phi - 2 * M_PI));
		const double diff_r = r - data.get_planetary_system().get_planet(1).get_r();
		if(rel_err < 1e-5){
		printf("(%d %d )	accel p	= (%.5e	%.5e)	%.5e	(%.5e %.5e)\n", n_rad, n_az, gradphi, -accel, rel_err, diff_r, diff_phi);
		}}
		/// END DEBUG
		*/

		const int cell_id = n_az + n_rad * ns;
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
		accel = ComputeDiskOnNbodyAccel(data, planet.get_x(), planet.get_y());
		planet.set_disk_on_planet_acceleration(accel);

		const double torque = (planet.get_x() * accel.y - planet.get_y() * accel.x)*planet.get_mass();
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
		const double dist = sqrt(std::pow(x - xo, 2) + std::pow(y - yo, 2));
		ax -= constants::G * mass / std::pow(dist, 3) * (x - xo);
		ay -= constants::G * mass / std::pow(dist, 3) * (y - yo);
	    }
	}
	planet.set_nbody_on_planet_acceleration_x(ax);
	planet.set_nbody_on_planet_acceleration_y(ay);
    }
}

/**
	Updates planets velocities due to disk influence if "DiskFeedback" is
   set.
*/
void UpdatePlanetVelocitiesWithDiskForce(t_data &data, double dt)
{
    for (unsigned int k = 0;
	 k < data.get_planetary_system().get_number_of_planets(); k++) {
	if (parameters::disk_feedback) {
	    t_planet &planet = data.get_planetary_system().get_planet(k);

		const Pair gamma = planet.get_disk_on_planet_acceleration();
	    const double new_vx =
		planet.get_vx() + dt * gamma.x + dt * IndirectTermDisk.x;
	    const double new_vy =
		planet.get_vy() + dt * gamma.y + dt * IndirectTermDisk.y;
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
