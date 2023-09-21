#include "Theo.h"
#include <cassert>
#include <cmath>
#include "parameters.h"

/**
	\param Density
*/
void RefillSigma(t_polargrid *Density)
{

	#pragma omp parallel for
    for (unsigned int nRadial = 0; nRadial < Density->Nrad; ++nRadial) {
    double mean = 0.0;
    for (unsigned int nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
        unsigned int cell = nAzimuthal + nRadial * Density->Nsec;
	    mean += Density->Field[cell];
	}
	mean /= (double)(Density->Nsec);
	SigmaMed[nRadial] = mean;
    }

    SigmaInf[0] = SigmaMed[0];

    for (unsigned int nRadial = 1; nRadial < Density->Nrad; ++nRadial) {
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

	#pragma omp parallel for
    for (unsigned int nRadial = 0; nRadial < Energy->Nrad; ++nRadial) {
    double mean = 0.0;
    for (unsigned int nAzimuthal = 0; nAzimuthal < Energy->Nsec; ++nAzimuthal) {
        unsigned int cell = nAzimuthal + nRadial * Energy->Nsec;
	    mean += Energy->Field[cell];
	}
	mean /= (double)Energy->Nsec;
	EnergyMed[nRadial] = mean;
    }
}
void init_binary_quadropole_moment(const t_planetary_system &psys){

    if(parameters::n_bodies_for_hydroframe_center == 2){ // binary
    const double a_b = psys.get_planet(1).get_semi_major_axis();
    const double m1 = psys.get_planet(0).get_mass();
    const double m2 = psys.get_planet(1).get_mass();
    double q_b;

    // assume that binary mass ratio < 1
    if(m2 < m1){
        q_b = m2/m1;
    } else { // m1 < m2
        q_b = m1/m2;
    }

    const double e_b = psys.get_planet(1).get_eccentricity();
    binary_quadropole_moment = std::pow(a_b,2) / 4.0 * q_b / std::pow((1.0+q_b), 2) * (1.0 + 3.0/2.0 * std::pow(e_b, 2));
    } else {
    binary_quadropole_moment = 0.0;
    }
}

///
/// \brief initial_energy, computes energy for the locally isothermal profile at radius R for central mass M
/// \param R
/// \param M
/// \return energy
///
double initial_energy(const double R, const double M){

	const double h0 = parameters::ASPECTRATIO_REF;
	const double F = parameters::FLARINGINDEX;
	const double S = parameters::SIGMASLOPE;

	// energy = 1/(gamma - 1) * Sigma * (C_s,iso)^2
	// C_s,iso = h * v_k = h0 * r^F * sqrt(G M / r)
	const double energy =
	1.0 / (parameters::ADIABATICINDEX - 1.0) * parameters::sigma0 *
	std::pow(h0, 2) * std::pow(R, - S - 1.0 + 2.0 * F) * constants::G * M;

	return energy;
}

///
/// \brief computes the pressure supported azimuthal velocity
/// around mass M at distance R for the locally isothermal model
/// \param R
/// \return locally_isothermal_v_az
///
double initial_locally_isothermal_v_az(const double R, const double M){

	// r_sm = sqrt(r**2 + (eps * H)**2)
	// d r_sm / dr = r/r_sm  *  (1 + (F+1) * (eps * H / r)**2)

	// Phi = GMm / r_sm
	// vkep^2 / r = 1/Simga dP/dr + dPhi / dr
	const double h0 = parameters::ASPECTRATIO_REF;
	const double F = parameters::FLARINGINDEX;
	const double S = parameters::SIGMASLOPE;
	const double h = h0 * std::pow(R, F);
	// const double eps = parameters::thickness_smoothing;
	const double vk_2 = constants::G * M / R;
	const double pressure_support_2 = (2.0 * F - 1.0 - S) * std::pow(h, 2);

	// for normal pressure support, the derivative should be 1
	const double smoothing_derivative_2 = 1.0;

	const double v_az = std::sqrt(vk_2 * (smoothing_derivative_2 + pressure_support_2));

	return v_az;
}


double support_azi_pressure(const double R){
	const double h0 = parameters::ASPECTRATIO_REF;
	const double F = parameters::FLARINGINDEX;
	const double S = parameters::SIGMASLOPE;
	const double h = h0 * std::pow(R, F);
	const double rv = (2.0 * F - 1.0 - S) * std::pow(h, 2);
	return rv;
}

double support_azi_smoothing_derivative(const double R) {
	const double h0 = parameters::ASPECTRATIO_REF;
	const double F = parameters::FLARINGINDEX;
	const double h = h0 * std::pow(R, F);
	const double eps = parameters::thickness_smoothing;
	const double rv = (1.0 + (F+1.0) * std::pow(h * eps, 2))
			/ std::pow(std::sqrt(1 + std::pow(h * eps, 2)), 3);
	return rv;
}

double support_azi_quadrupole(const double R) {
	const double Q = binary_quadropole_moment;
	double rv = 0.0;
	if (Q > 0.0) {
		rv = 3.0 * Q / std::pow(R, 2);
	}
	return rv;
}

///
/// \brief computes the pressure supported azimuthal velocity
/// around mass M at distance R for the locally isothermal model
/// also inclues the derivatives of the potential smoothing
/// \param R
/// \return locally_isothermal_smoothed_v_az
///
double initial_locally_isothermal_smoothed_v_az(const double R, const double M){

	// r_sm = sqrt(r**2 + (eps * H)**2)
	// d r_sm / dr = r/r_sm  *  (1 + (F+1) * (eps * H / r)**2)

	// Phi = GMm / r_sm
	// vkep^2 / r = 1/Simga dP/dr + dPhi / dr
	const double smoothing_derivative_2 = support_azi_smoothing_derivative(R);
	const double pressure_support_2 = support_azi_pressure(R);

    const double support = smoothing_derivative_2 + pressure_support_2;

	const double vk_2 = constants::G * M / R;
	return std::sqrt(vk_2 * support);
}


double initial_locally_isothermal_smoothed_v_az_with_quadropole_moment(const double R, const double M){

    // r_sm = sqrt(r**2 + (eps * H)**2)
    // d r_sm / dr = r/r_sm  *  (1 + (F+1) * (eps * H / r)**2)

    // Phi = GMm / r_sm
    // vkep^2 / r = 1/Simga dP/dr + dPhi / dr

    const double pressure_support_2 = support_azi_pressure(R);
    const double quadropole_support = support_azi_quadrupole(R);
	const double smoothing_derivative_2 = support_azi_smoothing_derivative(R);

    const double support = quadropole_support + smoothing_derivative_2 + pressure_support_2;

	const double vk_2 = constants::G * M / R;
	return std::sqrt(vk_2 * support);
}

///
/// \brief compute_v_kepler, compute kepler velocity around mass M at distance R
/// \param R
/// \param M
/// \return kepler velocity
///
double compute_v_kepler(const double R, const double M){
	const double vk = std::sqrt(constants::G * M / R);

	return vk;
}

///
/// \brief compute_viscous_radial_speed, v_radial from steady, viscous accretion disk
/// around mass M at radius R
/// \param R
/// \param M
/// \return vr, assuming kepler velocity profile for the disk
///
double initial_viscous_radial_speed(const double R, const double M){

	if (parameters::ALPHAVISCOSITY > 0){
	double sqrt_gamma;
	if(parameters::Adiabatic){
		sqrt_gamma = std::sqrt(parameters::ADIABATICINDEX);
	} else {
		sqrt_gamma = 1.0;
	}
	const double v_k = std::sqrt(constants::G * M / R);
	const double h =
	parameters::ASPECTRATIO_REF * std::pow(R, parameters::FLARINGINDEX);
	const double cs = sqrt_gamma * h * v_k;
	const double H = h * R;
	const double nu = parameters::ALPHAVISCOSITY * cs * H;
	const double vr = -3.0 * nu / R *
		  (-parameters::SIGMASLOPE + 2.0 * parameters::FLARINGINDEX + 1.0);

	return vr;
	} else {
		const double nu = parameters::VISCOSITY;
		const double vr = -3.0 * nu / R * (-parameters::SIGMASLOPE + .5);
		return vr;
	}
}

double calculate_omega_kepler(double r)
{
    return sqrt(constants::G * hydro_center_mass / (r * r * r));
}

double init_l1(const double central_star_mass, const double other_star_mass)
{
    const double q = central_star_mass / (central_star_mass + other_star_mass);

    double x = std::pow(other_star_mass / (3.0 * central_star_mass), 1.0 / 3.0);

    // Newton Raphson
    double f;
    double df;

    int counter = 0;
    do {
	counter++;
	if (counter > 10) {
	    break;
	}

	f = q / std::pow(1.0 - x, 2) - (1.0 - q) / std::pow(x, 2) - q + x;
	df = 2.0 * q / std::pow(1.0 - x, 3) + 2.0 * (1.0 - q) / std::pow(x, 3) +
	     1.0;

	x = x - f / df;

    } while (std::fabs(f) > 1e-14);

    return x;
}

/*	Compute the Roche Radius of a point mass (other_star_mass) orbiting a star (central_star_mass)
	using one Newton-Raphson iteration iteration.


	@param central_star_mass
	@param other_star_mass
	@param old l1
	@return roche radius
*/ 
double update_l1(const double central_star_mass, const double other_star_mass,
	       double l1)
{
    const double q = central_star_mass / (central_star_mass + other_star_mass);
    double x = l1;

    // Newton Raphson, one iteration
    double f = q / std::pow(1.0 - x, 2) - (1.0 - q) / std::pow(x, 2) - q + x;
    double df =
	2.0 * q / std::pow(1.0 - x, 3) + 2.0 * (1.0 - q) / std::pow(x, 3) + 1.0;

    assert(std::fabs(f) < 1e-8);

    x = x - f / df;

    return x;
}

double eggleton_1983(const double q, const double r)
{
    const double rL =
	0.49 * std::pow(q, 2.0 / 3.0) /
	(0.6 * std::pow(q, 2.0 / 3.0) + std::log(1.0 + std::pow(q, 1.0 / 3.0)));
    return rL * r;
}
