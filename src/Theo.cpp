#include "Theo.h"
#include <cassert>

/**
	\param Density
*/
void RefillSigma(t_polargrid *Density)
{
    unsigned int nRadial, nAzimuthal, cell;
    double mean;

    for (nRadial = 0; nRadial < Density->Nrad; ++nRadial) {
	mean = 0.0;
	for (nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
	    cell = nAzimuthal + nRadial * Density->Nsec;
	    mean += Density->Field[cell];
	}
	mean /= (double)(Density->Nsec);
	SigmaMed[nRadial] = mean;
    }

    SigmaInf[0] = SigmaMed[0];

    for (nRadial = 1; nRadial < Density->Nrad; ++nRadial) {
	SigmaInf[nRadial] =
	    (SigmaMed[nRadial - 1] * (Rmed[nRadial] - Rinf[nRadial]) +
	     SigmaMed[nRadial] * (Rinf[nRadial] - Rmed[nRadial - 1])) /
	    (Rmed[nRadial] - Rmed[nRadial - 1]);
    }
}

/**

*/
void RefillEnergy(t_polargrid *Energy)
{
    unsigned int nRadial, nAzimuthal, cell;
    double mean;

    for (nRadial = 0; nRadial < Energy->Nrad; ++nRadial) {
	mean = 0.0;
	for (nAzimuthal = 0; nAzimuthal < Energy->Nsec; ++nAzimuthal) {
	    cell = nAzimuthal + nRadial * Energy->Nsec;
	    mean += Energy->Field[cell];
	}
	mean /= (double)Energy->Nsec;
	EnergyMed[nRadial] = mean;
    }
}

double calculate_omega_kepler(double r)
{
    return sqrt(constants::G * hydro_center_mass / (r * r * r));
}

double init_l1(const double central_star_mass, const double other_star_mass){
	const double q = central_star_mass / (central_star_mass+other_star_mass);

	double x = std::pow(other_star_mass/(3.0*central_star_mass), 1.0/3.0);

	// Newton Raphson
	double f;
	double df;

	int counter = 0;
	do{
		counter ++;
		if(counter > 10){
			break;
		}

		f = q/std::pow(1.0 - x, 2) - (1.0 - q) / std::pow(x, 2) - q + x;
		df = 2.0*q/std::pow(1.0 - x, 3) + 2.0*(1.0 - q) / std::pow(x, 3) + 1.0;

		x = x - f/df;

	} while (std::fabs(f) > 1e-14);

	return x;
}


void update_l1(const double central_star_mass, const double other_star_mass, double &l1){
	const double q = central_star_mass / (central_star_mass+other_star_mass);

	double x = l1;

	// Newton Raphson, one iteration
	double f = q/std::pow(1.0 - x, 2) - (1.0 - q) / std::pow(x, 2) - q + x;
	double df = 2.0*q/std::pow(1.0 - x, 3) + 2.0*(1.0 - q) / std::pow(x, 3) + 1.0;

	assert(std::fabs(f) < 1e-8);

	x = x - f/df;

	l1 = x;
}


double eggleton_1983(const double q, const double r)
{
	const double rL = 0.49 * std::pow(q, 2.0/3.0) / (0.6 * std::pow(q, 2.0/3.0) + std::log(1.0 + std::pow(q, 1.0/3.0)));
	return rL * r;
}
