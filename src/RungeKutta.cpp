#include <math.h>
#include "util.h"

#include "RungeKutta.h"
#include "constants.h"
#include "types.h"

#define RKACCURACY 1e-1

static double k1[MAX1D], k2[MAX1D], k3[MAX1D], k4[MAX1D], k5[MAX1D], k6[MAX1D];
static double Dist[MAX1D];

void DerivMotionRK5(double* q_init, double* masses, double* deriv, int n, double dt, int* feelothers)
{
	double *x,*y,*vx,*vy, dist;
	double *derivx, *derivy, *derivvx, *derivvy;
	int i, j;

	x = q_init;
	y = x+n;
	vx = y+n;
	vy = vx+n;
	derivx = deriv;
	derivy = deriv+n;
	derivvx = derivy+n;
	derivvy = derivvx+n;
	for (i = 0; i < n; i++)
		Dist[i] = sqrt(pow2(x[i])+pow2(y[i]));
	for (i = 0; i < n; i++) {
		derivx[i] = vx[i];
		derivy[i] = vy[i];
		derivvx[i] = -constants::G*1.0/Dist[i]/Dist[i]/Dist[i]*x[i];
		derivvy[i] = -constants::G*1.0/Dist[i]/Dist[i]/Dist[i]*y[i];
		for (j = 0; j < n; j++) {
			// apply correction term caused by non barycenter systems
			derivvx[i] -= constants::G*masses[j]/pow3(Dist[j])*x[j];
			derivvy[i] -= constants::G*masses[j]/pow3(Dist[j])*y[j];
			if ((j != i) && (feelothers[i])) {
				dist = pow2(x[i]-x[j])+pow2(y[i]-y[j]);
				dist = sqrt(dist);
				derivvx[i] += constants::G*masses[j]/pow3(dist)*(x[j]-x[i]);
				derivvy[i] += constants::G*masses[j]/pow3(dist)*(y[j]-y[i]);
			}
		}
	}
	for (i = 0; i < 4*n; i++)
		deriv[i] *= dt;
}

void TranslatePlanetRK5(double* qold, double c1, double c2, double c3, double c4, double c5, double* qnew, int n)
{
	int i;

	for (i = 0; i < 4*n; i++)
		qnew[i] = qold[i]+c1*k1[i]+c2*k2[i]+c3*k3[i]+c4*k4[i]+c5*k5[i];
}

void RungeKutta(double* q0, double dt, double* masses, double* q1, int n, int* feelothers)
{
	int i;
	double timestep;

	timestep = dt;
	for (i = 0; i < n*4; i++) {
		k1[i] = k2[i] = k3[i] = k4[i] = k5[i] = k6[i];
	}
	DerivMotionRK5(q0, masses, k1, n, timestep, feelothers);
	TranslatePlanetRK5(q0, 0.2, 0.0, 0.0, 0.0, 0.0, q1, n);
	DerivMotionRK5(q1, masses, k2, n, timestep, feelothers);
	TranslatePlanetRK5(q0, 0.075, 0.225, 0.0, 0.0, 0.0, q1, n);
	DerivMotionRK5(q1, masses, k3, n, timestep, feelothers);
	TranslatePlanetRK5(q0, 0.3, -0.9, 1.2, 0.0, 0.0, q1, n);
	DerivMotionRK5(q1, masses, k4, n, timestep, feelothers);
	TranslatePlanetRK5(q0, -11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0, 0.0, q1, n);
	DerivMotionRK5(q1, masses, k5, n, timestep, feelothers);
	TranslatePlanetRK5(q0, 1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0, q1, n);
	DerivMotionRK5(q1, masses, k6, n, timestep, feelothers);
	for (i = 0; i < 4*n; i++) {
		q1[i]=q0[i]+37.0/378.0*k1[i]+250.0/621.0*k3[i]+125.0/594.0*k4[i]+512.0/1771.0*k6[i];
	}
}
