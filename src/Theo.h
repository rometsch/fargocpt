#ifndef THEO_H
#define THEO_H

#include "constants.h"
#include "global.h"
#include "types.h"
#include <cmath>

double calculate_omega_kepler(double r);
void RefillSigma(t_polargrid *Density);
void RefillEnergy(t_polargrid *Energy);
double eggleton_1983(const double q, const double r);
double init_l1(const double central_star_mass, const double other_star_mass);
void update_l1(const double central_star_mass, const double m2, double &l1);

#endif // THEO_H
