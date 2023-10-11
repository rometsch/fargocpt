#pragma once

#include <vector>

#include "data.h"
#include "global.h"

namespace pvte
{

void initializeLookupTables();

// hydrogen ionization fraction
double H_ionization_fraction(const double densityCGS, const double temperatureCGS);

// hydrogen dissociation fraction
double H_dissociation_fraction(const double densityCGS, const double temperatureCGS);

void compute_gamma_mu(t_data &data);

double get_gamma_eff(t_data &data, const int n_radial, const int n_azimuthal);

double get_mu(t_data &data, const int n_radial, const int n_azimuthal);

double get_gamma1(t_data &data, const int n_radial, const int n_azimuthal);
} // namespace pvte
