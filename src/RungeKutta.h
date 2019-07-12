#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

void RungeKutta(double* q0, double dt, double* masses, double* q1, unsigned int n, int* feelothers);
void TranslatePlanetRK5(double* qold, double c1, double c2, double c3, double c4, double c5, double* qnew,unsigned int n);
void DerivMotionRK5(double* q_init, double* masses, double* deriv,unsigned int n, double dt, int* feelothers);

#endif // RUNGEKUTTA_H
