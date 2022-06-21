#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

void RungeKutta(double *q0, double dt, double *masses, double *q1,
		unsigned int n);
void TranslatePlanetRK5(double *qold, const double c1, const double c2,
			const double c3, const double c4, const double c5,
			double *qnew, const unsigned int n, const double dt);
void DerivMotionRK5(double *q_init, double *masses, double *deriv,
		    unsigned int n);

#endif // RUNGEKUTTA_H
