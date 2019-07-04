/**
	\file Force.cpp

	Contains the function used to evaluate the %force due to disk, and the function that writes the 'tqwk' log files. Although the planet mass is given as an argument to ComputeForce(), this mass is used only to specify the distance cutoff in the case of the Hill sphere avoidance. The force returned is a specific force. It has therefore the dimension of an acceleration (LT^-2).
**/

#include <stdio.h>
#include <math.h>
#include <float.h>

#include "Force.h"
#include "constants.h"
#include "viscosity.h"
#include "LowTasks.h"
#include "logging.h"
#include "global.h"
#include "parameters.h"

Force *AllocateForce(int dimfxy)
{
	int i;
	Force *force;
	double *globalforce;

	force = (Force*)malloc(sizeof(Force));
	globalforce = (double*)malloc(sizeof(double) * 4 * dimfxy);

	if ( (force == NULL) || (globalforce == NULL) ) {
		logging::print(LOG_ERROR "Not Enough Memory.\n");
		PersonalExit(1);
	}

	for (i = 0; i < 4*dimfxy; i++)
		globalforce[i] = 0.;
	force->GlobalForce = globalforce;

	return force;
}

void FreeForce(Force* force)
{
	free(force->GlobalForce);
	free(force);
}

/**
	Computes the force due to the disk on an object at position (x,y)
*/
void ComputeForce(t_data &data, Force* force, double x, double y, double mass)
{
	int l, ns;
	double localforce[8]={0.,0.,0.,0.,0.,0.,0.,0.}, globalforce[8]={0.,0.,0.,0.,0.,0.,0.,0.};
	double xc, yc, cellmass, dx, dy, distance, dist2, rh, a;
	double InvDist3, fxi, fyi, fxhi, fyhi, fxo, fyo, fxho, fyho, hill_cut;
	double *abs, *ord;
	double rsmoothing = 0.0;

	ns = data[t_data::DENSITY].Nsec;
	abs = CellAbscissa->Field;
	ord = CellOrdinate->Field;
	fxi = fyi = fxhi = fyhi = fxo = fyo = fxho = fyho = 0.0;
	a = sqrt(x*x+y*y);
	rh = pow(mass/3.0, 1./3.)*a+DBL_EPSILON;

	// calculate smoothing length only once if not dependend on radius
	if (RocheSmoothing) {
		rsmoothing = rh*ROCHESMOOTHING;
	} else {
		// Thickness smoothing = smoothing with scale height
		rsmoothing = compute_smoothing(a);
	}
	// consider calculation of force on primary, don't need to recalculate rsmoothing then
	bool SmoothingEnabled = (a != 0.0);

	for (unsigned int n_radial = Zero_or_active; n_radial < Max_or_active; ++n_radial) {
		// calculate smoothing length if dependend on radius
		// i.e. for thickness smoothing with scale height at cell location
		if (SmoothingEnabled && ThicknessSmoothingAtCell) {
			rsmoothing = compute_smoothing(Rmed[n_radial]);
		}
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
			l = n_azimuthal+n_radial*ns;
			xc = abs[l];
			yc = ord[l];
			cellmass = Surf[n_radial]*data[t_data::DENSITY](n_radial,n_azimuthal);
			dx = xc-x;
			dy = yc-y;
			dist2 = dx*dx+dy*dy;
			hill_cut = 1.-exp(-dist2/(rh*rh));
			dist2 += rsmoothing*rsmoothing;
			distance = sqrt(dist2);
			InvDist3 = 1.0/dist2/distance;
			if (Rmed[n_radial] < a) {
				fxi += constants::G*cellmass*dx*InvDist3;
				fyi += constants::G*cellmass*dy*InvDist3;
				fxhi+= constants::G*cellmass*dx*InvDist3*hill_cut;
				fyhi+= constants::G*cellmass*dy*InvDist3*hill_cut;
			} else {
				fxo += constants::G*cellmass*dx*InvDist3;
				fyo += constants::G*cellmass*dy*InvDist3;
				fxho+= constants::G*cellmass*dx*InvDist3*hill_cut;
				fyho+= constants::G*cellmass*dy*InvDist3*hill_cut;
			}
		}
	}

	localforce[0] = fxi;
	localforce[1] = fyi;
	localforce[2] = fxhi;
	localforce[3] = fyhi;
	localforce[4] = fxo;
	localforce[5] = fyo;
	localforce[6] = fxho;
	localforce[7] = fyho;
	MPI_Allreduce(&localforce, &globalforce, 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	force->fx_inner = globalforce[0];
	force->fy_inner    = globalforce[1];
	force->fx_ex_inner = globalforce[2];
	force->fy_ex_inner = globalforce[3];
	force->fx_outer    = globalforce[4];
	force->fy_outer    = globalforce[5];
	force->fx_ex_outer = globalforce[6];
	force->fy_ex_outer = globalforce[7];
}

double compute_smoothing(double r)
{
	double smooth;
	smooth = parameters::thickness_smoothing * ASPECTRATIO * pow(r, 1.0+FLARINGINDEX);
	return smooth;
}

void UpdateLog(t_data &data, Force* fc, int outputnb, double time)
{
	double x, y, m, vx, vy;
	FILE *out;
	char filename[255];
	for (unsigned int i = 0; i < data.get_planetary_system().get_number_of_planets(); i++) {
		x = data.get_planetary_system().get_planet(i).get_x();
		y = data.get_planetary_system().get_planet(i).get_y();
		vx = data.get_planetary_system().get_planet(i).get_vx();
		vy = data.get_planetary_system().get_planet(i).get_vy();
		//r = sqrt(x*x+y*y);
		m = data.get_planetary_system().get_planet(i).get_mass();
		// if (RocheSmoothing)
		// 	smoothing = r*pow(m/3.,1./3.)*ROCHESMOOTHING;
		// else
		// 	smoothing = compute_smoothing(r);
        ComputeForce(data, fc, x, y, m);
		if (CPU_Rank == CPU_Number-1) {
			sprintf (filename, "%stqwk%d.dat", OUTPUTDIR, i);
			out = fopen (filename, "a");
			if (out == NULL) {
				fprintf (stderr, "Can't open %s\n", filename);
				fprintf (stderr, "Aborted.\n");
			}
			fprintf (out, "%d\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\n", outputnb,
				x*fc->fy_inner-y*fc->fx_inner,
				x*fc->fy_outer-y*fc->fx_outer,
				x*fc->fy_ex_inner-y*fc->fx_ex_inner,
				x*fc->fy_ex_outer-y*fc->fx_ex_outer,
				vx*fc->fx_inner+vy*fc->fy_inner,
				vx*fc->fx_outer+vy*fc->fy_outer,
				vx*fc->fx_ex_inner+vy*fc->fy_ex_inner,
				vx*fc->fx_ex_outer+vy*fc->fy_ex_outer, time);
			fclose (out);
		}
	}
}
