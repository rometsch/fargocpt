#ifndef PVTE_LAW_H
#define PVTE_LAW_H

#include <vector>

#include "data.h"
#include "global.h"

// double Ni;

// double Nj;

// double rhomin;
// double rhomax;

// double emin;
// double emax;

// double deltaLogRho;
// double deltaLogE;

namespace pvte
{

// struct for storing mean molecular weight, gamma_eff and gamma1
typedef struct {
    double mow, geff, g1;
} t_eosQuantities;

void makeZetaTables();

double get_funcDum(double temperatureCGS);

void initializeLookupTables();

// interpolate values from lookup table
double interpolate(const std::vector<double> &table, const int i, const int j,
		   const double x, const double y);

t_eosQuantities lookup(const double densityCGS, const double energyCGS);

// hydrogen ionization fraction
double H_ionization_fraction(const double densityCGS, const double temperatureCGS);

// hydrogen dissociation fraction
double H_dissociation_fraction(const double densityCGS, const double temperatureCGS);

// mean molecular weight mu
double mean_molecular_weight(const double temperatureCGS,
			     const double densityCGS);

// energy contributions to the internal energy of the gas
double gasEnergyContributions(const double x, const double y,
			      const double temperatureCGS);

// effective adiabatic index to relate pressure and internal energy
double gamma_eff(const double temperatureCGS, const double densityCGS);

// first adiabatic index to calculate the speed of sound
double gamma1(const double temperatureCGS, const double densityCGS);

// root finding problem for the calculation of the temperature
double gamma_mu_root(const double temperatureCGS, const double densityCGS,
		     const double energyCGS);

// solving the root finding problem
double energy_to_temperature(double energyCGS, double densityCGS);

double temperature_to_energy(const double temperatureCGS,
			     const double densityCGS);

void compute_gamma_mu(t_data &data);

double get_gamma_eff(t_data &data, const int n_radial, const int n_azimuthal);

double get_mu(t_data &data, const int n_radial, const int n_azimuthal);

double get_gamma1(t_data &data, const int n_radial, const int n_azimuthal);
} // namespace pvte

#endif // PVTE_LAW_H
