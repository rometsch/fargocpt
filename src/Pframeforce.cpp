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
extern boolean AllowAccretion, Corotating, Cooling;

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
	    double mass = planet.get_mass();
	    Pair accel = planet.get_disk_on_planet_acceleration();
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
	double mass = planet.get_mass();
	Pair accel = planet.get_nbody_on_planet_acceleration();
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
void CalculatePotential(t_data &data)
{
    double x, y, angle, distancesmooth;
    double smooth = 0.0;
    unsigned int number_of_planets =
	data.get_planetary_system().get_number_of_planets();
    std::vector<double> xpl(number_of_planets);
    std::vector<double> ypl(number_of_planets);
    std::vector<double> mpl(number_of_planets);
    std::vector<double> smooth_pl(number_of_planets);

    // setup planet data
    for (unsigned int k = 0; k < number_of_planets; k++) {
	t_planet &planet = data.get_planetary_system().get_planet(k);
	mpl[k] = data.get_planetary_system().get_planet(k).get_rampup_mass();
	xpl[k] = planet.get_x();
	ypl[k] = planet.get_y();
	if (RocheSmoothing) {
	    double r_hill = pow(mpl[k] / (3.0 * (M + mpl[k])), 1.0 / 3.0) *
			    planet.get_semi_major_axis();
	    smooth_pl[k] = pow2(r_hill * ROCHESMOOTHING);
	} else {
	    smooth_pl[k] = pow2(compute_smoothing_isothermal(planet.get_r()));
	}
    }

    data[t_data::POTENTIAL].clear();

    // gravitational potential from planets on gas
    for (unsigned int n_radial = 0;
	 n_radial <= data[t_data::POTENTIAL].get_max_radial(); ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::POTENTIAL].get_max_azimuthal();
	     ++n_azimuthal) {
	    angle = (double)n_azimuthal /
		    (double)data[t_data::POTENTIAL].get_size_azimuthal() * 2.0 *
		    PI;
	    x = Rmed[n_radial] * cos(angle);
	    y = Rmed[n_radial] * sin(angle);

	    for (unsigned int k = 0; k < number_of_planets; k++) {
		if (!ThicknessSmoothingAtCell) {
		    smooth = smooth_pl[k];
		}
		if (ThicknessSmoothingAtCell) {
		    smooth = pow2(compute_smoothing(Rmed[n_radial], data,
						    n_radial, n_azimuthal));
		}
		double distance2 = pow2(x - xpl[k]) + pow2(y - ypl[k]);
		if (k == 0)
		    smooth = 0.0;
		distancesmooth = sqrt(distance2 + smooth);
		// direct term from planet
		data[t_data::POTENTIAL](n_radial, n_azimuthal) +=
		    -constants::G * mpl[k] / distancesmooth;
	    }
	    // apply indirect term
	    // correct frame with contributions from disk and planets
	    data[t_data::POTENTIAL](n_radial, n_azimuthal) +=
		-IndirectTerm.x * x - IndirectTerm.y * y;
	}
    }
}

/**
   Update disk forces onto planets if *diskfeedback* is turned on
*/
void ComputeDiskOnNbodyAccel(Force *force, t_data &data)
{
    Pair accel;
    for (unsigned int k = 0;
	 k < data.get_planetary_system().get_number_of_planets(); k++) {
	if (parameters::disk_feedback) {
	    t_planet &planet = data.get_planetary_system().get_planet(k);
	    accel = ComputeAccel(force, data, planet.get_x(), planet.get_y(),
				 planet.get_mass());
	    planet.set_disk_on_planet_acceleration(accel);
	}
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
	double x = planet.get_x();
	double y = planet.get_y();
	double ax = 0.0;
	double ay = 0.0;
	for (unsigned int nother = 0;
	     nother < planetary_system.get_number_of_planets(); nother++) {
	    if (nother != npl) {
		t_planet &other_planet = planetary_system.get_planet(nother);
		double xo = other_planet.get_x();
		double yo = other_planet.get_y();
		double mass = other_planet.get_mass();
		double dist = sqrt(pow2(x - xo) + pow2(y - yo));
		ax -= constants::G * mass / pow3(dist) * (x - xo);
		ay -= constants::G * mass / pow3(dist) * (y - yo);
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
    Pair gamma;

    for (unsigned int k = 0;
	 k < data.get_planetary_system().get_number_of_planets(); k++) {
	if (parameters::disk_feedback) {
	    t_planet &planet = data.get_planetary_system().get_planet(k);

	    gamma = planet.get_disk_on_planet_acceleration();
	    double new_vx =
		planet.get_vx() + dt * gamma.x + dt * IndirectTermDisk.x;
	    double new_vy =
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

void AccreteOntoPlanets(t_data &data, double dt)
{
    bool masses_changed = false;
    double RRoche, Rplanet, distance, dx, dy, deltaM, angle, temp;
    int j_min, j_max, j, l, jf, ns, lip, ljp;
    unsigned int i, i_min, i_max, nr;
    double Xplanet, Yplanet, Mplanet, VXplanet, VYplanet;
    double facc, facc1, facc2, frac1,
	frac2; /* We adopt the same notations as W. Kley */
    double *dens, *vrad, *vtheta;
    double PxPlanet, PyPlanet, vrcell, vtcell, vxcell, vycell, xc, yc;
    double dPxPlanet, dPyPlanet, dMplanet;
    nr = data[t_data::DENSITY].Nrad;
    ns = data[t_data::DENSITY].Nsec;
    dens = data[t_data::DENSITY].Field;
    const double* cell_center_x = CellCenterX->Field;
    const double* cell_center_y = CellCenterY->Field;
    vrad = data[t_data::V_RADIAL].Field;
    vtheta = data[t_data::V_AZIMUTHAL].Field;

    auto &planetary_system = data.get_planetary_system();

    for (unsigned int k = 0; k < planetary_system.get_number_of_planets();
	 k++) {
	auto &planet = planetary_system.get_planet(k);
	if (planet.get_acc() > 1e-10) {
	    dMplanet = dPxPlanet = dPyPlanet = 0.0;
	    // Hereafter : initialization of W. Kley's parameters
	    // remove a ratio of facc = planet.get_acc() of the mass inside the
	    // Hill sphere every free fall time at the Hill radius
	    facc = dt * planet.get_acc() * planet.get_omega() * sqrt(12.0) /
		   2.0 / PI;
	    facc1 = 1.0 / 3.0 * facc;
	    facc2 = 2.0 / 3.0 * facc;
	    frac1 = 0.75;
	    frac2 = 0.45;

	    // W. Kley's parameters initialization finished
	    Xplanet = planet.get_x();
	    Yplanet = planet.get_y();
	    VXplanet = planet.get_vx();
	    VYplanet = planet.get_vy();
	    Mplanet = planet.get_mass();
	    Rplanet = sqrt(Xplanet * Xplanet + Yplanet * Yplanet);
	    RRoche = pow((1.0 / 3.0 * Mplanet), (1.0 / 3.0)) * Rplanet;

	    // Central mass is 1.0
	    i_min = 0;
	    i_max = nr - 1;
	    while ((Rsup[i_min] < Rplanet - RRoche) && (i_min < nr))
		i_min++;
	    while ((Rinf[i_max] > Rplanet + RRoche) && (i_max > 0))
		i_max--;
	    angle = atan2(Yplanet, Xplanet);
		j_min =
		(int)((double)ns / 2.0 / PI * (angle - 2.0 * RRoche / Rplanet));
	    j_max =
		(int)((double)ns / 2.0 / PI * (angle + 2.0 * RRoche / Rplanet));

		unsigned int i_min_test = clamp_id_to_grid(get_rinf_id(parameters::radial_grid_type, Rplanet - RRoche));
		unsigned int i_max_test = clamp_id_to_grid(get_rinf_id(parameters::radial_grid_type, Rplanet + RRoche));

		unsigned int j_min_test = get_inf_azimuthal_id(angle - 2.0 * RRoche / Rplanet);
		unsigned int j_max_test = get_inf_azimuthal_id(angle + 2.0 * RRoche / Rplanet);

		if(i_min != i_min_test)
			printf("i_min is wrong: %d	%d\n", i_min, i_min_test);

		if(i_max != i_max_test)
			printf("i_max is wrong: %d	%d\n", i_max, i_max_test);


		double dj_min = ((double)ns / 2.0 / PI * (angle - 2.0 * RRoche / Rplanet));
		double djf_min = std::floor(dj_min);

		printf("Floor : %g	%g\n", dj_min, djf_min);

		if(j_min != j_min_test)
			printf("j_min is wrong: %d	%d\n", j_min, j_min_test);

		if(j_max != j_max_test)
			printf("j_max is wrong: %d	%d\n", j_max, j_max_test);



	    PxPlanet = Mplanet * VXplanet;
	    PyPlanet = Mplanet * VYplanet;
	    for (i = i_min; i <= i_max; i++) {
		for (j = j_min; j <= j_max; j++) {
		    jf = j;
		    while (jf < 0)
			jf += ns;
		    while (jf >= ns)
			jf -= ns;
		    l = jf + i * ns;
		    lip = l + ns;
		    ljp = l + 1;
		    if (jf == ns - 1)
			ljp = i * ns;
		    xc = cell_center_x[l];
		    yc = cell_center_y[l];
		    dx = Xplanet - xc;
		    dy = Yplanet - yc;
		    distance = sqrt(dx * dx + dy * dy);
		    vtcell =
			0.5 * (vtheta[l] + vtheta[ljp]) + Rmed[i] * OmegaFrame;
		    vrcell = 0.5 * (vrad[l] + vrad[lip]);
		    vxcell = (vrcell * xc - vtcell * yc) / Rmed[i];
		    vycell = (vrcell * yc + vtcell * xc) / Rmed[i];
		    if (distance < frac1 * RRoche) {
			deltaM = facc1 * dens[l] * Surf[i];
			if (i < Zero_or_active)
			    deltaM = 0.0;
			if (i >= Max_or_active)
			    deltaM = 0.0;
			dens[l] *= (1.0 - facc1);
			dPxPlanet += deltaM * vxcell;
			dPyPlanet += deltaM * vycell;
			dMplanet += deltaM;
		    }
		    if (distance < frac2 * RRoche) {
			deltaM = facc2 * dens[l] * Surf[i];
			if (i < Zero_or_active)
			    deltaM = 0.0;
			if (i >= Max_or_active)
			    deltaM = 0.0;
			dens[l] *= (1.0 - facc2);
			dPxPlanet += deltaM * vxcell;
			dPyPlanet += deltaM * vycell;
			dMplanet += deltaM;
		    }
		}
	    }
	    MPI_Allreduce(&dMplanet, &temp, 1, MPI_DOUBLE, MPI_SUM,
			  MPI_COMM_WORLD);
	    dMplanet = temp;
	    MPI_Allreduce(&dPxPlanet, &temp, 1, MPI_DOUBLE, MPI_SUM,
			  MPI_COMM_WORLD);
	    dPxPlanet = temp;
	    MPI_Allreduce(&dPyPlanet, &temp, 1, MPI_DOUBLE, MPI_SUM,
			  MPI_COMM_WORLD);
	    dPyPlanet = temp;
	    PxPlanet += dPxPlanet;
	    PyPlanet += dPyPlanet;
	    Mplanet += dMplanet;
	    if (parameters::disk_feedback) {
		planet.set_vx(PxPlanet / Mplanet);
		planet.set_vy(PyPlanet / Mplanet);
	    }
	    planet.set_mass(Mplanet);
	    if (!masses_changed)
		masses_changed = true;
	    // update hydro center mass
	}
    }
    if (masses_changed) {
	planetary_system.update_global_hydro_frame_center_mass();
    }
}
