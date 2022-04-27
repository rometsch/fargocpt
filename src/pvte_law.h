#ifndef PVTE_LAW_H
#define PVTE_LAW_H

#include <vector>

#include "global.h"
#include "data.h"

//double Ni;

//double Nj;

//double rhomin;
//double rhomax;

//double emin; 
//double emax;

//double deltaLogRho; 
//double deltaLogE;



namespace pvte
{

//struct for storing mean molecular weight, gammaeff and gamma1
typedef struct t_eosQuantities {
    double mow, geff, g1;
};

void  makeZetaTables();

double get_funcDum(double temperatureCGS);

void initializeLookupTables(std::vector<double> &mu_table, std::vector<double> &gammaeff_table,
std::vector<double> &gamma1_table);

//interpolate values from lookup table
double interpolate(std::vector<double> &table,int i ,int j,double x,double y);

t_eosQuantities lookup(const std::vector<double>& mu_tab, const std::vector<double>& gammeff_tab,
const std::vector<double>& gamm1_tab, double densityCGS, double energyCGS);

//hydrogen ionization fraction
double Hfraction (const double densityCGS, const double temperatureCGS);

//hydrogen dissociation fraction
double H2fraction (double densityCGS, double temperatureCGS);

//mean molecular weight mu
double mean_molecular_weight (double temperatureCGS, double densityCGS);

//energy contributions to the internal energy of the gas
double gasEnergyContributions(double xMF, double x, 
double y, double temperatureCGS);

//effective adiabatic index to relate pressure and internal energy
double gammaeff(double temperatureCGS, double densityCGS);

//first adiabatic index to calculate the speed of sound
double gamma1(double temperatureCGS, double densityCGS);

//root finding problem for the calculation of the temperature
double gamma_mu_root(double temperatureCGS, double densityCGS, double energyCGS);

//solving the root finding problem
double energy_to_temperature(double energyCGS, double densityCGS);

double temperature_to_energy(double temperatureCGS, double densityCGS);

void compute_gamma_mu(t_data &data);

double get_gammaeff(t_data &data, int n_radial, int n_azimuthal);

double get_mu(t_data &data, int n_radial, int n_azimuthal);

double get_gamma1(t_data &data, int n_radial, int n_azimuthal);
}



#endif // PVTE_LAW_H
