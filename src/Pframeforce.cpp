/**
	\file Pframeforce.cpp

	Functions that evaluate the %force between the planets and the disk. The FillForcesArrays() function is ill-named: it rather fill an array of the potential, that is later use to derive the force acting on the disk at every zone. The name of this file is due to the fact that we work in the frame centered on the primary (which is therefore not inertial). Old versions of fargo also featured the file Gframeforce.c, in which we worked in the frame centered on the center of gravity of the system.  The present file also contains the functions necessary to update the planets' position and velocities (taking into account, or not, the indirect term, ie the term arising from the fact that the frame is not inertial), as well as the function that initializes the hydrodynamics fields with analytic prescription.
*/

#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "parameters.h"
#include "Pframeforce.h"
#include "global.h"
#include "SideEuler.h"
#include "logging.h"
#include "viscosity.h"
#include "Theo.h"
#include "RungeKutta.h"
#include "Force.h"
#include "selfgravity.h"
#include "axilib.h"
#include "constants.h"
#include "units.h"
#include "LowTasks.h"
#include "util.h"

extern Pair IndirectTerm;
extern Pair IndirectTermDisk;
extern Pair IndirectTermPlanets;
extern boolean AllowAccretion, Corotating, Cooling;
static double q0[MAX1D], q1[MAX1D], PlanetMasses[MAX1D];
static int FeelOthers[MAX1D];

void ComputeIndirectTerm(Force* force,t_data &data) {
	IndirectTerm.x = 0.0;
	IndirectTerm.y = 0.0;
	IndirectTermDisk.x = 0.0;
	IndirectTermDisk.y = 0.0;
	IndirectTermPlanets.x = 0.0;
	IndirectTermPlanets.y = 0.0;

	// compute disk indirect term
	if (parameters::disk_feedback) {
		if (parameters::no_default_star) {
			// add up contributions from disk on all bodies used to calculate the center
			double mass_center = 0.0;
			for (unsigned int n=0; n<parameters::n_bodies_for_barycenter; n++) {
				t_planet &planet = data.get_planetary_system().get_planet(n);
				double mass =  planet.get_mass();
				Pair &accel = planet.get_disk_on_planet_acceleration();
				IndirectTermDisk.x -= mass*accel.x;
				IndirectTermDisk.y -= mass*accel.y;
				mass_center += mass;
			}
			IndirectTermDisk.x /= mass_center;
			IndirectTermDisk.y /= mass_center;
		} else {
			// default mode with primary star in coordinate center
			IndirectTermDisk = ComputeAccel(force, data, 0.0, 0.0, 0.0);
			IndirectTermDisk.x = -IndirectTermDisk.x;
			IndirectTermDisk.y = -IndirectTermDisk.y;
		}
	}

	// compute nbody indirect term
	if (parameters::no_default_star) {
		// add up contributions from mutual interactions from all bodies used to calculate the center
		double mass_center = 0.0;
		for (unsigned int n=0; n<parameters::n_bodies_for_barycenter; n++) {
			t_planet &planet = data.get_planetary_system().get_planet(n);
			double mass =  planet.get_mass();
			Pair &accel = planet.get_nbody_on_planet_acceleration();
			IndirectTermPlanets.x -= mass*accel.x;
			IndirectTermPlanets.y -= mass*accel.y;
			mass_center += mass;
		}
		IndirectTermPlanets.x /= mass_center;
		IndirectTermPlanets.y /= mass_center;
	} else {
		// default mode with primary star in coordinate center
		for (unsigned int k = 0; k < data.get_planetary_system().get_number_of_planets(); k++) {
			t_planet &planet = data.get_planetary_system().get_planet(k);
			double InvPlanetDistance3 =  1.0/pow3(planet.get_distance());
			double mplanet = data.get_planetary_system().get_planet(k).get_mass();
			IndirectTermPlanets.x = -constants::G*mplanet*InvPlanetDistance3*planet.get_x();
			IndirectTermPlanets.y = -constants::G*mplanet*InvPlanetDistance3*planet.get_y();
		}
	}

	IndirectTerm.x = IndirectTermDisk.x + IndirectTermPlanets.x;
	IndirectTerm.y = IndirectTermDisk.y + IndirectTermPlanets.y;
}

/* Below : work in non-rotating frame */
void FillForcesArrays(t_data &data)
{
	double x, y, angle, distancesmooth;;
	double smooth, mplanet;
	double pot;
	double InvDistance;
	//double xbin, ybin, mbin, distbin, Invdistbin3;

	data[t_data::POTENTIAL].clear();

	// gravitational potential from planets on gas
	for (unsigned int k = 0;k  < data.get_planetary_system().get_number_of_planets(); k++) {
		t_planet &planet = data.get_planetary_system().get_planet(k);
		mplanet = data.get_planetary_system().get_planet(k).get_mass();

		if (data.get_planetary_system().get_planet(k).get_rampuptime() > 0) {
			double ramping = 1.0;
			if (PhysicalTime < data.get_planetary_system().get_planet(k).get_rampuptime()*DT) {
				ramping = 1.0-pow2(cos(PhysicalTime*PI/2.0/(data.get_planetary_system().get_planet(k).get_rampuptime()*DT)));
			}
			mplanet *= ramping;
		}

		// hill radius
		double r_hill = pow(planet.get_mass()/(3.0*(M+planet.get_mass())),1.0/3.0)*planet.get_semi_major_axis();

		// calculate smoothing length only once if not dependend on radius
		if (RocheSmoothing) {
			smooth = pow2(r_hill*ROCHESMOOTHING);
		} else {
			smooth = pow2(compute_smoothing(planet.get_distance()));
		}

		for (unsigned int n_radial = 0; n_radial <= data[t_data::POTENTIAL].get_max_radial(); ++n_radial) {
			InvDistance = 1.0/Rmed[n_radial];
			// calculate smoothing length if dependend on radius
			// i.e. for thickness smoothing with scale height at cell location
			if (ThicknessSmoothingAtCell) {
				smooth = pow2(compute_smoothing(Rmed[n_radial]));
			}
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::POTENTIAL].get_max_azimuthal(); ++n_azimuthal) {
				angle = (double)n_azimuthal/(double)data[t_data::POTENTIAL].get_size_azimuthal()*2.0*PI;
				x = Rmed[n_radial]*cos(angle);
				y = Rmed[n_radial]*sin(angle);
				double distance2 = pow2(x-planet.get_x())+pow2(y-planet.get_y());
				distancesmooth = sqrt(distance2+smooth);

				pot = -constants::G*mplanet/distancesmooth;
				data[t_data::POTENTIAL](n_radial,n_azimuthal) += pot;
			}
		}
	}

	// gravitational potential from star on gas
	for (unsigned int n_radial = 0; n_radial <= data[t_data::POTENTIAL].get_max_radial(); ++n_radial) {
		InvDistance = 1.0/Rmed[n_radial];
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::POTENTIAL].get_max_azimuthal(); ++n_azimuthal) {
			angle = (double)n_azimuthal/(double)data[t_data::POTENTIAL].get_size_azimuthal()*2.0*PI;
			x = Rmed[n_radial]*cos(angle);
			y = Rmed[n_radial]*sin(angle);

			// direct term from star
			pot = -constants::G*1.0*InvDistance;
			// correct frame with contributions from disk and planets
			pot -= IndirectTerm.x*x + IndirectTerm.y*y;

			data[t_data::POTENTIAL](n_radial,n_azimuthal) += pot;
		}
	}
}

/**
   Update disk forces onto planets if *diskfeedback* is turned on
*/
void ComputeDiskOnNbodyAccel(Force* force,t_data &data)
{
	Pair accel;
	for (unsigned int k = 0; k < data.get_planetary_system().get_number_of_planets(); k++) {
		if (data.get_planetary_system().get_planet(k).get_feeldisk()) {
			t_planet &planet = data.get_planetary_system().get_planet(k);
			accel = ComputeAccel(force, data, planet.get_x(), planet.get_y(), planet.get_mass());
			planet.set_disk_on_planet_acceleration(accel);
		}
	}
}

/**
   Update mutual planet-planet accelerations
*/
void ComputeNbodyOnNbodyAccel(t_planetary_system &planetary_system)
{

	for (unsigned int npl = 0; npl < data.get_planetary_system().get_number_of_planets(); npl++) {
		t_planet &planet = data.get_planetary_system().get_planet(npl);
		double x = planet.get_x();
		double y = planet.get_y();
		double ax = 0.0;
		double ay = 0.0;
		for (unsigned int nother = 0; nother < data.get_planetary_system().get_number_of_planets(); nother++) {
			if (nother != npl) {
				t_planet &other_planet = data.get_planetary_system().get_planet(nother);
				double xo = other_planet.get_x();
				double yo = other_planet.get_y();
				double mass = other_planet.get_mass();
				double dist = sqrt( pow2(x - xo) + pow2(y - yo) );
				ax -= constants::G*mass/pow3(dist)*(x-xo);
				ay -= constants::G*mass/pow3(dist)*(y-yo);
			}
		}
		planet.set_nbody_on_planet_acceleration_x(ax);
		planet.set_nbody_on_planet_acceleration_y(ay);
	}
}

/**
	Updates planets velocities due to disk influence if "feeldisk" is set for the planet
*/
void AdvanceSystemFromDisk(Force* force, t_data &data, double dt)
{
	Pair gamma;

	for (unsigned int k = 0; k < data.get_planetary_system().get_number_of_planets(); k++) {
		if (data.get_planetary_system().get_planet(k).get_feeldisk()) {
			t_planet &planet = data.get_planetary_system().get_planet(k);

			gamma = planet.get_disk_on_planet_acceleration();
			double new_vx = planet.get_vx() + dt * gamma.x + dt * IndirectTermDisk.x;
			double new_vy = planet.get_vy() + dt * gamma.y + dt * IndirectTermDisk.y;
			planet.set_vx(new_vx);
			planet.set_vy(new_vy);
		}
	}
}

void AdvanceSystemRK5(t_data &data, double dt)
{
	unsigned int n = data.get_planetary_system().get_number_of_planets();

	if (parameters::integrate_planets) {
		for (unsigned int i = 0; i < data.get_planetary_system().get_number_of_planets(); i++) {
			q0[i+0*n] = data.get_planetary_system().get_planet(i).get_x();
			q0[i+1*n] = data.get_planetary_system().get_planet(i).get_y();
			q0[i+2*n] = data.get_planetary_system().get_planet(i).get_vx();
			q0[i+3*n] = data.get_planetary_system().get_planet(i).get_vy();
			PlanetMasses[i] = data.get_planetary_system().get_planet(i).get_mass();
			FeelOthers[i] = data.get_planetary_system().get_planet(i).get_feelother();
		}

		RungeKutta(q0, dt, PlanetMasses, q1, n, FeelOthers);

		for (unsigned int i = 0; i < data.get_planetary_system().get_number_of_planets(); i++) {
			data.get_planetary_system().get_planet(i).set_x(q1[i+0*n]);
			data.get_planetary_system().get_planet(i).set_y(q1[i+1*n]);
			data.get_planetary_system().get_planet(i).set_vx(q1[i+2*n]);
			data.get_planetary_system().get_planet(i).set_vy(q1[i+3*n]);
		}

	}
}

void SolveOrbits(t_data &data)
{
	double x, y, vx, vy;
	for (unsigned int i = 0; i < data.get_planetary_system().get_number_of_planets(); i++) {
		x = data.get_planetary_system().get_planet(i).get_x();
		y = data.get_planetary_system().get_planet(i).get_y();
		vx = data.get_planetary_system().get_planet(i).get_vx();
		vy = data.get_planetary_system().get_planet(i).get_vy();
		FindOrbitalElements(x, y, vx, vy, 1.0+data.get_planetary_system().get_planet(i).get_mass(), i);
	}
}

double ConstructSequence(double* u, double* v, int n)
{
	int i;
	double lapl=0.0;
	for (i = 1; i < n; i++) {
		u[i] = 2.0*v[i]-u[i-1];
	}
	for (i = 1; i < n-1; i++) {
		lapl += fabs(u[i+1]+u[i-1]-2.0*u[i]);
	}
	return lapl;
}

void AccreteOntoPlanets(t_data &data, double dt)
{
	double RRoche, Rplanet, distance, dx, dy, deltaM, angle, temp;
	int j_min, j_max, j, l, jf, ns, lip, ljp;
	unsigned int i, i_min,i_max, nr;
	double Xplanet, Yplanet, Mplanet, VXplanet, VYplanet;
	double facc, facc1, facc2, frac1, frac2; /* We adopt the same notations as W. Kley */
	double *dens, *abs, *ord, *vrad, *vtheta;
	double PxPlanet, PyPlanet, vrcell, vtcell, vxcell, vycell, xc, yc;
	double dPxPlanet, dPyPlanet, dMplanet;
	nr     = data[t_data::DENSITY].Nrad;
	ns     = data[t_data::DENSITY].Nsec;
	dens   = data[t_data::DENSITY].Field;
	abs    = CellAbscissa->Field;
	ord    = CellOrdinate->Field;
	vrad   = data[t_data::V_RADIAL].Field;
	vtheta = data[t_data::V_AZIMUTHAL].Field;

	for (unsigned int k=0; k < data.get_planetary_system().get_number_of_planets(); k++) {
		if (data.get_planetary_system().get_planet(k).get_acc() > 1e-10) {
			dMplanet = dPxPlanet = dPyPlanet = 0.0;
			// Hereafter : initialization of W. Kley's parameters
			facc = dt*data.get_planetary_system().get_planet(k).get_acc();
			facc1 = 1.0/3.0*facc;
			facc2 = 2.0/3.0*facc;
			frac1 = 0.75;
			frac2 = 0.45;

			// W. Kley's parameters initialization finished
			Xplanet = data.get_planetary_system().get_planet(k).get_x();
			Yplanet = data.get_planetary_system().get_planet(k).get_y();
			VXplanet = data.get_planetary_system().get_planet(k).get_vx();
			VYplanet = data.get_planetary_system().get_planet(k).get_vy();
			Mplanet = data.get_planetary_system().get_planet(k).get_mass();
			Rplanet = sqrt(Xplanet*Xplanet+Yplanet*Yplanet);
			RRoche = pow((1.0/3.0*Mplanet),(1.0/3.0))*Rplanet;

			// Central mass is 1.0
			i_min=0;
			i_max=nr-1;
			while ((Rsup[i_min] < Rplanet-RRoche) && (i_min < nr)) i_min++;
			while ((Rinf[i_max] > Rplanet+RRoche) && (i_max > 0)) i_max--;
			angle = atan2 (Yplanet, Xplanet);
			j_min =(int)((double)ns/2.0/PI*(angle - 2.0*RRoche/Rplanet));
			j_max =(int)((double)ns/2.0/PI*(angle + 2.0*RRoche/Rplanet));
			PxPlanet = Mplanet*VXplanet;
			PyPlanet = Mplanet*VYplanet;
			for (i = i_min; i <= i_max; i++) {
				for (j = j_min; j <= j_max; j++) {
					jf = j;
					while (jf <  0)  jf += ns;
					while (jf >= ns) jf -= ns;
					l   = jf+i*ns;
					lip = l+ns;
					ljp = l+1;
					if (jf == ns-1) ljp = i*ns;
					xc = abs[l];
					yc = ord[l];
					dx = Xplanet-xc;
					dy = Yplanet-yc;
					distance = sqrt(dx*dx+dy*dy);
					vtcell=0.5*(vtheta[l]+vtheta[ljp])+Rmed[i]*OmegaFrame;
					vrcell=0.5*(vrad[l]+vrad[lip]);
					vxcell=(vrcell*xc-vtcell*yc)/Rmed[i];
					vycell=(vrcell*yc+vtcell*xc)/Rmed[i];
					if (distance < frac1*RRoche) {
						deltaM = facc1*dens[l]*Surf[i];
						if (i < Zero_or_active) deltaM = 0.0;
						if (i >= Max_or_active) deltaM = 0.0;
						dens[l] *= (1.0 - facc1);
						dPxPlanet    += deltaM*vxcell;
						dPyPlanet    += deltaM*vycell;
						dMplanet     += deltaM;
					}
					if (distance < frac2*RRoche) {
						deltaM = facc2*dens[l]*Surf[i];
						if (i < Zero_or_active) deltaM = 0.0;
						if (i >= Max_or_active) deltaM = 0.0;
						dens[l] *= (1.0 - facc2);
						dPxPlanet    += deltaM*vxcell;
						dPyPlanet    += deltaM*vycell;
						dMplanet     += deltaM;
					}
				}
			}
			MPI_Allreduce (&dMplanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			dMplanet = temp;
			MPI_Allreduce (&dPxPlanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			dPxPlanet = temp;
			MPI_Allreduce (&dPyPlanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			dPyPlanet = temp;
			PxPlanet += dPxPlanet;
			PyPlanet += dPyPlanet;
			Mplanet  += dMplanet;
			if (data.get_planetary_system().get_planet(k).get_feeldisk()) {
				data.get_planetary_system().get_planet(k).set_vx(PxPlanet/Mplanet);
				data.get_planetary_system().get_planet(k).set_vy(PyPlanet/Mplanet);
			}
			data.get_planetary_system().get_planet(k).set_mass(Mplanet);
		}
	}
}

void FindOrbitalElements(double x, double y, double vx, double vy, double m, int n)
{
	// m is mstar + mplanet
	double Ax, Ay, e, d, h, a, E, M, V;
	double PerihelionPA;
	FILE *output;
	double temp;
	char name[256];

	if (CPU_Rank != CPU_Highest) {
		return;
	}

	sprintf(name, "%sorbit%d.dat", OUTPUTDIR, n);
	output = fopen (name, "a");

	if (output == NULL) {
		logging::print(LOG_ERROR "Can't open 'orbit%d.dat'. Exited.\n",n);
		PersonalExit(1);
	}

	h = x*vy-y*vx;
	d = sqrt(x*x+y*y);
	Ax = x*vy*vy-y*vx*vy-constants::G*m*x/d;
	Ay = y*vx*vx-x*vx*vy-constants::G*m*y/d;
	e = sqrt(Ax*Ax+Ay*Ay)/constants::G/m;
	a = h*h/constants::G/m/(1.0-e*e);

	if (e != 0.0) {
		temp = (1.0-d/a)/e;
		if (temp > 1.0) {
			// E = acos(1)
			E = 0.0;
		} else if (temp < -1.0) {
			// E = acos(-1)
			E = PI;
		} else {
			E = acos(temp);
		}
	} else {
		E = 0.0;
	}

	if ((x*y*(vy*vy-vx*vx)+vx*vy*(x*x-y*y)) < 0) {
		E= -E;
	}

	M = E-e*sin(E);

	if (e != 0.0) {
		temp = (a*(1.0-e*e)/d-1.0)/e;
		if (temp > 1.0) {
			// V = acos(1)
			V = 0.0;
		} else if (temp < -1.0) {
			// V = acos(-1)
			V = PI;
		} else {
			V = acos(temp);
		}
	} else {
		V = 0.0;
	}

	if (E < 0.0) {
		V = -V;
	}

	if (e != 0.0) {
		PerihelionPA=atan2(Ay,Ax);
	} else {
		PerihelionPA=atan2(y,x);
	}

	fprintf (output, "%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n", PhysicalTime, e, a, M, V, PerihelionPA);
	fclose (output);
}
