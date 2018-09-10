#ifndef _RUNGEKUTTA_H_
#define _RUNGEKUTTA_H_

void RungeKutta(double* q0, double dt, double* masses, double* q1, int n, int* feelothers);
void TranslatePlanetRK5(double* qold, double c1, double c2, double c3, double c4, double c5, double* qnew, int n);
void DerivMotionRK5(double* q_init, double* masses, double* deriv, int n, double dt, int* feelothers);

#endif
