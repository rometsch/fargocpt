#include <math.h>
#include "util.h"

#include "RungeKutta.h"
#include "constants.h"
#include "parameters.h"
#include "types.h"
#include <vector>

#define RKACCURACY 1e-1

static double k1[MAX1D], k2[MAX1D], k3[MAX1D], k4[MAX1D], k5[MAX1D], k6[MAX1D];
static double Dist[MAX1D];

void DerivMotionRK5(double* q_init, double* masses, double* deriv,unsigned int n, double dt, int* feelothers)
{
	double *x,*y,*vx,*vy, dist;
	double *deriv_x, *deriv_y, *deriv_vx, *deriv_vy;
	double indirect_term_x = 0.0;
	double indirect_term_y = 0.0;
	std::vector<double> accel_x(n);
	std::vector<double> accel_y(n);
	unsigned int number_of_planets = n;

	// select range of q_init array for local variables
	x = q_init;
	y = x+n;
	vx = y+n;
	vy = vx+n;
	deriv_x = deriv;
	deriv_y = deriv+n;
	deriv_vx = deriv_y+n;
	deriv_vy = deriv_vx+n;

	// calculate mutual forces
	for (unsigned int npl = 0; npl < number_of_planets; npl++) {
		double ax = 0.0;
		double ay = 0.0;
		for (unsigned int nother = 0; nother < number_of_planets; nother++) {
			if (nother != npl) {
				dist = sqrt( pow2(x[npl] - x[nother]) + pow2(y[npl] - y[nother]) );
				ax += -constants::G*masses[nother]/pow3(dist)*(x[npl]-x[nother]);
				ay += -constants::G*masses[nother]/pow3(dist)*(y[npl]-y[nother]);
			}
		}
		accel_x[npl] = ax;
		accel_y[npl] = ay;
	}
	// calculate indirect term
	double mass_center = 0.0;
	for (unsigned int i=0; i<parameters::n_bodies_for_hydroframe_center; i++) {
		indirect_term_x -= masses[i]*accel_x[i];
		indirect_term_y -= masses[i]*accel_y[i];
		mass_center += masses[i];
	}
	indirect_term_x /= mass_center;
	indirect_term_y /= mass_center;

	// apply accelerations
	for (unsigned int i=0; i<number_of_planets; i++) {
		deriv_x[i] = vx[i];
		deriv_y[i] = vy[i];
		deriv_vx[i] = accel_x[i] + indirect_term_x;
		deriv_vy[i] = accel_y[i] + indirect_term_y;
	}

	// multiply all derivatives by dt
	for (unsigned int i = 0; i < 4*n; i++)
		deriv[i] *= dt;
}

void TranslatePlanetRK5(double* qold, double c1, double c2, double c3, double c4, double c5, double* qnew,unsigned int n)
{
	for (unsigned int i = 0; i < 4*n; i++)
		qnew[i] = qold[i]+c1*k1[i]+c2*k2[i]+c3*k3[i]+c4*k4[i]+c5*k5[i];
}

void RungeKutta(double* q0, double dt, double* masses, double* q1,unsigned int n, int* feelothers)
{
	double timestep;

	timestep = dt;
	for (unsigned int i = 0; i < n*4; i++) {
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
	for (unsigned int i = 0; i < 4*n; i++) {
		q1[i]=q0[i]+37.0/378.0*k1[i]+250.0/621.0*k3[i]+125.0/594.0*k4[i]+512.0/1771.0*k6[i];
	}
}
