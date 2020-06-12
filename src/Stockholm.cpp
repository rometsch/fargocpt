#include <math.h>
#include <stdio.h>

#include "Force.h"
#include "Stockholm.h"
#include "constants.h"
#include "global.h"
#include "parameters.h"

Force ComputeForceStockholm(t_data &data, double x, double y, double rsmoothing,
			    double mass)
{
    int j, l, ns;
    double localforce[8] = {0., 0., 0., 0., 0., 0., 0., 0.},
	   globalforce[8] = {0., 0., 0., 0., 0., 0., 0., 0.};
    double xc, yc, cellmass, dx, dy, distance, dist2, rh, a;
    double InvDist3, fxi, fyi, fxhi, fyhi, fxo, fyo, fxho, fyho, outside_hill,
	inside_hill;
    unsigned int i;

    Force Force;
    ns = data[t_data::DENSITY].Nsec;
    const double* dens = data[t_data::DENSITY].Field;
    const double* cell_center_x = CellCenterX->Field;
    const double* cell_center_y = CellCenterY->Field;
    fxi = fyi = fxhi = fyhi = fxo = fyo = fxho = fyho = 0.0;
    a = sqrt(x * x + y * y);
    rh = pow(mass / 3., 1. / 3.) * a + 1e-15;

    for (i = Zero_or_active; i < Max_or_active; i++) {
	for (j = 0; j < ns; j++) {
	    l = j + i * ns;
	    xc = cell_center_x[l];
	    yc = cell_center_y[l];
	    cellmass = Surf[i] * dens[l];
	    dx = xc - x;
	    dy = yc - y;
	    dist2 = dx * dx + dy * dy;
	    outside_hill = (dist2 >= rh * rh ? 1.0 : 0.0);
	    inside_hill =
		(((dist2 >= 0.25 * rh * rh) && (dist2 < rh * rh)) ? 1.0 : 0.0);
	    dist2 += rsmoothing * rsmoothing;
	    distance = sqrt(dist2);
	    InvDist3 = 1.0 / dist2 / distance;
	    if (Rmed[i] < a) {
		fxhi += constants::G * cellmass * dx * InvDist3 * inside_hill;
		fyhi += constants::G * cellmass * dy * InvDist3 * inside_hill;
		fxi += constants::G * cellmass * dx * InvDist3 * outside_hill;
		fyi += constants::G * cellmass * dy * InvDist3 * outside_hill;
	    } else {
		fxho += constants::G * cellmass * dx * InvDist3 * inside_hill;
		fyho += constants::G * cellmass * dy * InvDist3 * inside_hill;
		fxo += constants::G * cellmass * dx * InvDist3 * outside_hill;
		fyo += constants::G * cellmass * dy * InvDist3 * outside_hill;
	    }
	}
    }

    globalforce[0] = fxi;
    globalforce[1] = fyi;
    globalforce[2] = fxhi;
    globalforce[3] = fyhi;
    globalforce[4] = fxo;
    globalforce[5] = fyo;
    globalforce[6] = fxho;
    globalforce[7] = fyho;

    MPI_Allreduce(&localforce, &globalforce, 8, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);

    Force.fx_inner = globalforce[0];
    Force.fy_inner = globalforce[1];
    Force.fx_ex_inner = globalforce[2];
    Force.fy_ex_inner = globalforce[3];
    Force.fx_outer = globalforce[4];
    Force.fy_outer = globalforce[5];
    Force.fx_ex_outer = globalforce[6];
    Force.fy_ex_outer = globalforce[7];

    Force.GlobalForce = NULL;

    return Force;
}

Pair MassInOut(t_data &data, double a)
{
    int j, l, ns;
    unsigned int i;
    double localmass[2] = {0., 0.}, globalmass[2] = {0., 0.};
    double massin = 0.0, massout = 0.0, outside = 0.0, inside = 0.0;
    double *dens, cellmass;
    Pair massinout;

    ns = data[t_data::DENSITY].Nsec;
    dens = data[t_data::DENSITY].Field;

    for (i = Zero_or_active; i < Max_or_active; i++) {
	if (Rinf[i] > a) {
	    outside = 1.0;
	    inside = 0.0;
	}
	if (Rsup[i] < a) {
	    outside = 0.0;
	    inside = 1.0;
	}
	if ((Rinf[i] < a) && (Rsup[i] > a)) {
	    outside = (Rsup[i] * Rsup[i] - a * a) /
		      (Rsup[i] * Rsup[i] - Rinf[i] * Rinf[i]);
	    inside = (-Rinf[i] * Rinf[i] + a * a) /
		     (Rsup[i] * Rsup[i] - Rinf[i] * Rinf[i]);
	}
	for (j = 0; j < ns; j++) {
	    l = j + i * ns;
	    cellmass = Surf[i] * dens[l];
	    massin += cellmass * inside;
	    massout += cellmass * outside;
	}
    }

    localmass[0] = massin;
    localmass[1] = massout;
    MPI_Allreduce(&localmass, &globalmass, 2, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);

    massinout.x = globalmass[0];
    massinout.y = globalmass[1];
    return massinout;
}

void UpdateLogStockholm(t_data &data, double time)
{
    double x, y, r, m, smoothing /*, iplanet, cs, frac*/;
    Force fc;
    Pair massinout;
    FILE *out;
    char filename[255];
    for (unsigned int i = 0;
	 i < data.get_planetary_system().get_number_of_planets(); i++) {
	x = data.get_planetary_system().get_planet(i).get_x();
	y = data.get_planetary_system().get_planet(i).get_y();
	r = sqrt(x * x + y * y);
	m = data.get_planetary_system().get_planet(i).get_mass();
	if (parameters::roche_smoothing_enabled)
	    smoothing = r * pow(m / 3., 1. / 3.) * parameters::roche_smoothing_factor;
	else
	    smoothing = compute_smoothing_isothermal(r);
	fc = ComputeForceStockholm(data, x, y, smoothing, m);
	massinout = MassInOut(data, r);
	if (CPU_Rank == CPU_Number - 1) {
	    sprintf(filename, "%storque%d.dat", OUTPUTDIR.c_str(), i);
	    out = fopen(filename, "a");
	    if (out == NULL) {
		fprintf(stderr, "Can't open %s\n", filename);
		fprintf(stderr, "Aborted.\n");
	    }
	    fprintf(out, "%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\n",
		    time, massinout.x, massinout.y,
		    x * fc.fy_inner - y * fc.fx_inner,
		    x * fc.fy_outer - y * fc.fx_outer,
		    x * fc.fy_ex_inner - y * fc.fx_ex_inner,
		    x * fc.fy_ex_outer - y * fc.fx_ex_outer);
	    fclose(out);
	}
    }
}
