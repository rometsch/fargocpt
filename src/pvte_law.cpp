#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <vector>

#include "constants.h"
#include "data.h"
#include "logging.h"
#include "parameters.h"
#include "pvte_law.h"

namespace pvte
{
// Ni = number of density gridpoints
const int Ni = 1000;

// Nj = number of energy gridpoints
const int Nj = 1000;

// smallest and largest density for the lookup tables
const double rhomin = 1.0e-23;
const double rhomax = 1.0;

// smalles and largest energy for the lookup tables
const double emin = 1.0e8;
const double emax = 1.0e15;

// logarithmic grid spacing
const double deltaLogRho = std::log10(rhomax / rhomin) / (double)Ni;
const double deltaLogE = std::log10(emax / emin) / (double)Nj;

// parameters for zeta tables
const double ORTHO_PARA_MODE = 1;

const double THETA_V = 6140.0;
const double THETA_R = 85.5;

const int Nzeta = 5000;
const double Temp0 = 1.0;
const double Tmax = 1.0e12;

// lookup tables for the pve equation of state
std::vector<double> rho_table(Ni);
std::vector<double> e_table(Nj);
std::vector<double> mu_table(Ni *Nj);
std::vector<double> gamma_eff_table(Ni *Nj);
std::vector<double> gamma1_table(Ni *Nj);
std::vector<double> lnT(Nzeta);
std::vector<double> funcdum(Nzeta);

void makeZetaTables()
{
    int i, j;
    double dy = std::log(Tmax / Temp0) * (1. / (double)Nzeta);
    double a, b, b1, scrh, inv_T2;
    double zetaP, dzetaP, zetaR, dzetaR;
    double dum1, dum2, dum3;
    double alpha, beta, gamma;
    double db, sum1, sum2;

    logging::print_master(LOG_INFO " generating Zeta tables...\n");
    if (ORTHO_PARA_MODE == 0) {
	alpha = 1.0;
	beta = 0.0;
	gamma = 0.0;
    } else if (ORTHO_PARA_MODE == 2) {
	alpha = 0.25;
	beta = 0.75;
	gamma = 0.0;
    } else {
	alpha = 1.0;
	beta = 0.0;
	gamma = 1.0;
    }

    b1 = 2.0 * THETA_R;
    for (j = 0; j < Nzeta; j++) {
	const double T = Temp0 * std::exp(j * dy);
	inv_T2 = 1.0 / (T * T);
	zetaP = dzetaP = 0.0;
	sum1 = sum2 = 0.0;
	for (i = 0; i <= 10000; i++) {
	    a = 2 * i + 1;
	    b = i * (i + 1) * THETA_R;
	    if ((i % 2) == 0) {
		scrh = a * std::exp(-b / T);
		zetaP += scrh;
		dzetaP += scrh * b;
	    } else {
		db = b - b1;
		scrh = a * std::exp(-db / T);
		sum1 += scrh;
		sum2 += scrh * db;
	    }
	}

	dzetaP *= inv_T2;
	const double zetaO = std::exp(-b1 / T) * sum1;
	const double dzetaO = std::exp(-b1 / T) * (b1 * sum1 + sum2) * inv_T2;
	const double dzO_zO_m = sum2 / sum1 * inv_T2;
	lnT[j] = std::log(T);

	scrh = zetaO * std::exp(2.0 * THETA_R / T);

	zetaR =
	    std::pow(zetaP, alpha) * std::pow(scrh, beta) + 3.0 * gamma * zetaO;
	dzetaR = (zetaR - 3.0 * gamma * zetaO) *
		     (alpha * (dzetaP / zetaP) + beta * dzO_zO_m) +
		 3.0 * gamma * dzetaO;
	dum1 = THETA_V / T;
	dum2 = dum1 * std::exp(-dum1) / (1.0 - std::exp(-dum1));
	dum3 = (T / zetaR) * dzetaR;
	funcdum[j] = 1.5 + dum2 + dum3;
    }
}

double get_funcDum(double temperatureCGS)
{
    const double y = std::log(temperatureCGS);

    if (y > lnT[Nzeta - 2]) {
	const double funcdum_val = funcdum[Nzeta - 2];
	return funcdum_val;
    } else if (y < lnT[0]) {
	const double funcdum_val = funcdum[0];
	return funcdum_val;
    } else {
	const double dy = lnT[1] - lnT[0];
	const int indx = (int)(std::floor((y - lnT[0]) / dy));

	if (indx >= Nzeta || indx < 0) {
	    logging::print_master(
		LOG_ERROR "! GetFuncDum: indx out of range, indx = %d\n", indx);
	}
	const double funcdum_val = (funcdum[indx] * (lnT[indx + 1] - y) +
				    funcdum[indx + 1] * (y - lnT[indx])) /
				   dy;
	return funcdum_val;
    }
}

void initializeLookupTables()
{
    makeZetaTables();

    for (int i = 0; i < Ni; ++i) {
	for (int j = 0; j < Nj; ++j) {
	    double rhoi = std::pow(10.0, (deltaLogRho * i)) * rhomin;
	    double ej = std::pow(10.0, (deltaLogE * j)) * emin;
	    double T = energy_to_temperature(ej, rhoi);
	    double mu = mean_molecular_weight(T, rhoi);
	    double geff = gamma_eff(T, rhoi);
	    double g1 = gamma1(T, rhoi);

	    int index = j + i * Nj;
	    rho_table[i] = rhoi;
	    e_table[j] = ej;
	    mu_table[index] = mu;
	    gamma_eff_table[index] = geff;
	    gamma1_table[index] = g1;
	}
    }
}

// interpolate values from lookup table
double interpolate(const std::vector<double> &table, const int i, const int j,
		   const double x, const double y)
{
    const int ind1 = j + (i + 1) * Nj;
    const int ind2 = j + i * Nj;
    const int ind3 = j + 1 + (i + 1) * Nj;
    const int ind4 = j + 1 + i * Nj;
    double S_ij = table[ind1] * x + table[ind2] * (1.0 - x);
    double S_ijp1 = table[ind3] * x + table[ind4] * (1.0 - x);
    return S_ij * (1.0 - y) + S_ijp1 * y;
}

// get the interpolated quantities
t_eosQuantities lookup(const double densityCGS, const double energyCGS)
{
    int i = int(std::floor(std::log10(densityCGS / rhomin) / deltaLogRho));
    int j = int(std::floor(std::log10(energyCGS / emin) / deltaLogE));
    if (i >= Ni - 1) {
	i = Ni - 2;
    }
    if (i < 0) {
	i = 0;
    }
    if (j >= Nj - 1) {
	j = Nj - 2;
    }
    if (j < 0) {
	j = 0;
    }

    const double rhoi = rho_table[i];
    const double rhoip1 = rho_table[i + 1];
    const double ej = e_table[j];
    const double ejp1 = e_table[j + 1];
    const double x = (densityCGS - rhoi) / (rhoip1 - rhoi);
    const double y = (energyCGS - ej) / (ejp1 - ej);
    const double geff = interpolate(gamma_eff_table, i, j, x, y);
    const double mu = interpolate(mu_table, i, j, x, y);
    const double gamma1 = interpolate(gamma1_table, i, j, x, y);
    t_eosQuantities result;
    result.geff = geff;
    result.mow = mu;
    result.g1 = gamma1;
    return result;
}

// hydrogen ionization fraction
double Hfraction(const double densityCGS, const double temperatureCGS)
{
    double x = 1.0;
    double Ax = 5.3890039e-09 * std::pow(temperatureCGS, 1.5) *
		std::exp(-157821.45 / temperatureCGS) / densityCGS;
    if (Ax < 1.0e8) {
	x = 0.5 * (-Ax + std::sqrt(Ax * Ax + 4.0 * Ax));
    }
    return x;
}

// hydrogen dissociation fraction
double H2fraction(const double densityCGS, const double temperatureCGS)
{
    double y = 1.0;
    double Ay = 7.5035099e-05 * std::pow(temperatureCGS, 1.5) *
		std::exp(-51988.24 / temperatureCGS) / densityCGS;
    if (Ay < 1.0e8) {
	y = 0.5 * (-Ay + std::sqrt(Ay * Ay + 4.0 * Ay));
    }
    return y;
}

// mean molecular weight mu
double mean_molecular_weight(const double temperatureCGS,
			     const double densityCGS)
{
    double xMF = 0.75;

    double x = Hfraction(densityCGS, temperatureCGS);

    double y = H2fraction(densityCGS, temperatureCGS);

    return 4.0 / (2.0 * xMF * (1.0 + y + 2.0 * y * x) + 1.0 - xMF);
}

// energy contributions to the internal energy of the gas
double gasEnergyContributions(const double xMF, const double x, const double y,
			      const double temperatureCGS)
{
    const double epsH2 = 0.5 * xMF * (1.0 - y) * get_funcDum(temperatureCGS);
    const double epsHII = 157821.45 * xMF * x * y / temperatureCGS;
    const double epsHH = 25994.12 * xMF * y / temperatureCGS;
    const double epsHe = 0.375 * (1.0 - xMF);
    const double epsHI = 1.5 * xMF * (1.0 + x) * y;
    return epsH2 + epsHII + epsHH + epsHe + epsHI;
}

// effective adiabatic index to relate pressure and internal energy
double gamma_eff(const double temperatureCGS, const double densityCGS)
{
    // hydrogen mass fraction xMF
    const double xMF = 0.75;

    // hydrogen ionization fraction x
    const double x = Hfraction(densityCGS, temperatureCGS);

    // hydrogen dissociation fraction y
    const double y = H2fraction(densityCGS, temperatureCGS);

    const double mu = mean_molecular_weight(temperatureCGS, densityCGS);

    const double gamma_eff =
	1.0 + 1.0 / (mu * gasEnergyContributions(xMF, x, y, temperatureCGS));

    return gamma_eff;
}

// first adiabatic index to calculate the speed of sound
double gamma1(const double temperatureCGS, const double densityCGS)
{
    const double xMF = 0.75;
    const double epsilon = 1.0e-4;
    const double temperatureLeft = temperatureCGS * (1.0 - epsilon);
    const double temperatureRight = temperatureCGS * (1.0 + epsilon);
    const double deltaTemperature = temperatureLeft - temperatureRight;
    // hydrogen ionization fraction x
    double xL = Hfraction(densityCGS, temperatureLeft);
    double xR = Hfraction(densityCGS, temperatureRight);

    const double xc = Hfraction(densityCGS, temperatureCGS);

    // hydrogen dissociation fraction y
    double yL = H2fraction(densityCGS, temperatureLeft);
    double yR = H2fraction(densityCGS, temperatureRight);

    const double yc = H2fraction(densityCGS, temperatureCGS);

    // contributions to the internal energy
    const double eps = gasEnergyContributions(xMF, xc, yc, temperatureCGS);

    const double eL = (gasEnergyContributions(xMF, xL, yL, temperatureLeft)) *
		      temperatureLeft;
    const double eR = (gasEnergyContributions(xMF, xR, yR, temperatureRight)) *
		      temperatureRight;
    const double e = eps * temperatureCGS;

    const double cv = (eL - eR) / deltaTemperature;

    double muL = 4.0 / (2.0 * xMF * (1.0 + yL + 2.0 * yL * xL) + 1.0 - xMF);
    double muR = 4.0 / (2.0 * xMF * (1.0 + yR + 2.0 * yR * xR) + 1.0 - xMF);
    double muc = 4.0 / (2.0 * xMF * (1.0 + yc + 2.0 * yc * xc) + 1.0 - xMF);

    const double gamma_eff = 1.0 + 1.0 / (muc * eps);

    const double p = (gamma_eff - 1.0) * e;

    const double chiT =
	1.0 - temperatureCGS / muc * (muL - muR) / deltaTemperature;

    // Derivative with respect to the density
    const double rhoL = densityCGS * (1.0 - epsilon);
    const double rhoR = densityCGS * (1.0 + epsilon);
    const double deltaRho = rhoL - rhoR;

    // hydrogen ionization fraction x
    xL = Hfraction(rhoL, temperatureCGS);
    xR = Hfraction(rhoR, temperatureCGS);

    // hydrogen dissociation fraction y
    yL = H2fraction(rhoL, temperatureCGS);
    yR = H2fraction(rhoR, temperatureCGS);

    muL = 4.0 / (2.0 * xMF * (1.0 + yL + 2.0 * yL * xL) + 1.0 - xMF);
    muR = 4.0 / (2.0 * xMF * (1.0 + yR + 2.0 * yR * xR) + 1.0 - xMF);
    muc = 4.0 / (2.0 * xMF * (1.0 + yc + 2.0 * yc * xc) + 1.0 - xMF);

    double chiRho = 1.0 - densityCGS / muc * (muL - muR) / deltaRho;

    return p * std::pow(chiT, 2) / (cv * temperatureCGS) + chiRho;
}

// root finding problem for the calculation of the temperature
double gamma_mu_root(const double temperatureCGS, const double densityCGS,
		     const double energyCGS)
{
    // hydrogen mass fraction xMF
    const double xMF = 0.75;

    // hydrogen ionization fraction x
    const double x = Hfraction(densityCGS, temperatureCGS);

    // hydrogen dissociation fraction y
    const double y = H2fraction(densityCGS, temperatureCGS);

    const double mu = 4.0 / (2.0 * xMF * (1.0 + y + 2.0 * y * x) + 1.0 - xMF);

    const double gamma =
	1.0 + 1.0 / (mu * gasEnergyContributions(xMF, x, y, temperatureCGS));

    const double temperature =
	mu * energyCGS * (gamma - 1.0) / constants::_R.get_cgs_value();

    return temperature - temperatureCGS;
}

// solving the root finding problem
double energy_to_temperature(double energyCGS, double densityCGS)
{
    // Brent's root finding method
    const double delta = 1.0e-3;
    double a = 1.0e0; // parameters::minimum_temperature;
    double b = 1.0e7; // parameters::maximum_temperature;
    double c;
    double d;
    double s;
    double fa = gamma_mu_root(a, densityCGS, energyCGS);
    double fb = gamma_mu_root(b, densityCGS, energyCGS);
    double fs;
    volatile double
	fc; // otherwise fc is optimized away and the root finding crashes.

    if (std::abs(fa) < std::abs(fb)) {
	std::swap(a, b);
	std::swap(fa, fb);
    }
    c = a;
    fc = fa;
    bool mflag = true;
    while (std::abs(b - a) > delta) {
	if ((fa != fc) && (fb != fc)) {
	    s = a * fb * fc / ((fa - fb) * (fa - fc)) +
		b * fa * fc / ((fb - fa) * (fb - fc)) +
		c * fa * fb / ((fc - fa) * (fc - fb));
	} else {
	    s = b - fb * (b - a) / (fb - fa);
	}

	if (((s < std::min((3.0 * a + b) / 4.0, b)) &&
	     (s > std::max((3.0 * a + b) / 4.0, b))) ||
	    (mflag && (std::abs(s - b) >= std::abs(b - c) / 2.0)) ||
	    (!mflag && (std::abs(s - b) >= std::abs(c - d) / 2.0)) ||
	    (mflag && (std::abs(b - c) < delta)) ||
	    (!mflag and (std::abs(c - d) < delta))) {
	    s = (a + b) / 2.0;
	    mflag = true;
	} else {
	    mflag = false;
	}
	fs = gamma_mu_root(s, densityCGS, energyCGS);
	d = c;
	c = b;
	if (fa * fs < 0.0) {
	    b = s;
	} else {
	    a = s;
	}

	if (std::abs(fa) < std::abs(fb)) {
	    std::swap(a, b);
	    std::swap(fa, fb);
	}
    }
    return b;
}

double temperature_to_energy(const double temperatureCGS,
			     const double densityCGS)
{
    const double gamma = gamma_eff(temperatureCGS, densityCGS);
    const double mu = mean_molecular_weight(temperatureCGS, densityCGS);
    return temperatureCGS * constants::_R.get_cgs_value() /
	   (mu * (gamma - 1.0));
}

void compute_gamma_mu(t_data &data)
{
    /*
	static bool lookupTablesInitialized = false;

    //generate lookup tables once
	if (!lookupTablesInitialized){
		logging::print_master(LOG_INFO "Generating lookup tables \n");
		initializeLookupTables(mu_table, gammeff_table, gamma1_table);
		lookupTablesInitialized = true;
		logging::print_master(LOG_INFO "Lookup tables generated \n");
	}
    */
    for (unsigned int n_radial = 0;
	 n_radial < data[t_data::SIGMA].get_size_radial(); ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::SIGMA].get_size_azimuthal();
	     ++n_azimuthal) {

	    const double sigma = data[t_data::SIGMA](n_radial, n_azimuthal);

	    const double H = data[t_data::SCALE_HEIGHT](n_radial, n_azimuthal);

	    double densityCGS, energyCGS;

        // TODO: check units!
	    if (parameters::ShockTube > 0) {
		densityCGS = sigma * units::density.get_cgs_factor();
		energyCGS = data[t_data::ENERGY](n_radial, n_azimuthal) *
			    units::energy_density.get_cgs_factor() / densityCGS;
	    } else {
		densityCGS =
		    sigma / (parameters::density_factor * H) * units::density;

		energyCGS = data[t_data::ENERGY](n_radial, n_azimuthal) *
			    units::energy_density /
			    (sigma * units::surface_density);
	    }

	    t_eosQuantities q = lookup(densityCGS, energyCGS);

	    data[t_data::GAMMAEFF](n_radial, n_azimuthal) = q.geff;
	    data[t_data::MU](n_radial, n_azimuthal) = q.mow;
	    data[t_data::GAMMA1](n_radial, n_azimuthal) = q.g1;
	}
    }
}

double get_gamma_eff(t_data &data, const int n_radial, const int n_azimuthal)
{
    if (parameters::variableGamma) {
	return data[t_data::GAMMAEFF](n_radial, n_azimuthal);
    } else {
	return parameters::ADIABATICINDEX;
    }
}

double get_mu(t_data &data, const int n_radial, const int n_azimuthal)
{
    if (parameters::variableGamma) {
	return data[t_data::MU](n_radial, n_azimuthal);
    } else {
	return parameters::MU;
    }
}

double get_gamma1(t_data &data, const int n_radial, const int n_azimuthal)
{
    if (parameters::variableGamma) {
	return data[t_data::GAMMA1](n_radial, n_azimuthal);
    } else {
	return parameters::ADIABATICINDEX;
    }
}
} // namespace pvte
