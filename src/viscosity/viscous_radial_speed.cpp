#include "viscous_radial_speed.h"

#include "../parameters.h"
#include "../Theo.h"
#include "../global.h"
#include "../types.h"
#include "../util.h"
#include "../polargrid.h"

///
/// \brief get_nu, compute initial viscosity for the locally isothermal model
/// \param r
/// \param M
/// \return
///
static double get_nu(const double R, const double M){

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
	return nu;
}

/// does exact same as get_nu, expect it also ensures that
/// the soundspeed is inside the allowed Temperature range.
static double get_nu2(const double R, const double M, const double Sigma){


	const double v_k = std::sqrt(constants::G * M / R);
	const double h =
	parameters::ASPECTRATIO_REF * std::pow(R, parameters::FLARINGINDEX);

	double cutoff = 1.0;
	if (parameters::profile_cutoff_outer) {
		cutoff *= cutoff_outer(parameters::profile_cutoff_point_outer,
							parameters::profile_cutoff_width_outer, R);
	}
	if (parameters::profile_cutoff_inner) {
		cutoff *= 	cutoff_inner(parameters::profile_cutoff_point_inner,
						parameters::profile_cutoff_width_inner, R);
	}

	double cs_adb;
	double H;
	if(parameters::Adiabatic){
		const double gamma = parameters::ADIABATICINDEX;
		double energy = cutoff * 1.0/(gamma-1.0) * Sigma * std::pow(h*v_k, 2);

		const double temperature_floor =
			parameters::minimum_temperature *
			units::temperature.get_inverse_cgs_factor();
		const double energy_floor =
			temperature_floor *
			Sigma /
			parameters::MU * constants::R / (gamma - 1.0);

		const double temperature_ceil =
			parameters::maximum_temperature *
			units::temperature.get_inverse_cgs_factor();
		const double energy_ceil =
			temperature_ceil *
			Sigma /
			parameters::MU * constants::R / (gamma - 1.0);
		energy = std::max(energy, energy_floor);
		energy = std::min(energy, energy_ceil);

		cs_adb = std::sqrt(gamma * (gamma-1.0) * energy/Sigma);
		const double cs_iso = std::sqrt((gamma-1.0) * energy/Sigma);
		const double omega_k = v_k / R;
		H = cs_iso / omega_k;

	} else { // locally isothermal
		cs_adb = h*v_k;
		H = h*R;
	}

	const double nu = parameters::ALPHAVISCOSITY * cs_adb * H;
	return nu;
}

static double get_sigma(const double R){
	const double S = parameters::SIGMASLOPE;
	double density =
		parameters::sigma0 * std::pow(R, -S);

	const double density_floor =
		parameters::sigma_floor * parameters::sigma0;

	if (parameters::profile_cutoff_outer) {
		density *= cutoff_outer(parameters::profile_cutoff_point_outer,
							parameters::profile_cutoff_width_outer, R);
	}

	if (parameters::profile_cutoff_inner) {
		density *= 	cutoff_inner(parameters::profile_cutoff_point_inner,
						parameters::profile_cutoff_width_inner, R);
	}

	density = std::max(density, density_floor);

	return density;
}

namespace viscous_speed
{

/// \brief derive_5th_order, see https://en.wikipedia.org/wiki/Numerical_differentiation
/// \param r, position variable
/// \param mass, parameter
/// \param f, function to derive
/// \return df / dr at r
double derive(const double r, const double mass, double (*f)(double, double)){

	const double x = r;
	const double h = 8.0e-4 * x;

	const double f1 = -1.0 * f(x+2.0*h, mass);
	const double f2 =  8.0 * f(x+h, mass);
	const double f3 = -8.0 * f(x-h, mass);
	const double f4 =  1.0 * f(x-2.0*h, mass);

	const double derivative = (f1 + f2 + f3 + f4) / (12.0 * h);

	return derivative;
}

///
/// \brief get_w
/// \return pressure supported, smoothed potential angular velocity (omega)
///
double get_w(const double r, const double mass){
	return initial_locally_isothermal_smoothed_v_az(r, mass) / r;
}

///
/// \brief get_r2_w
/// \return r^2 * omega
///
double get_r2_w(const double r, const double mass){
	const double omega = get_w(r,mass);
	return std::pow(r, 2) * omega;
}

///
/// \brief get_nu_S_r3_dwdr
/// \param r
/// \param mass
/// \return nu * Sigma * r^3 * dw/dr
///
double get_nu_S_r3_dwdr(const double r, const double mass){

	const double dw_dr = derive(r, mass, get_w);
	const double Sigma = get_sigma(r);
	const double nu = get_nu2(r, mass, Sigma);
	const double nu_S_r3_dwdr = nu * Sigma * std::pow(r, 3) * dw_dr;
	return nu_S_r3_dwdr;
}

/// \brief get_vr_with_numerical_viscous_speed
/// \param Rcenter, position of returned vr
/// \param Rinterface, position of the cell interface
/// vp: Rcenter == Rmed -> Rinterface = Rsup
/// vr: Rcenter == Rinf -> Rinterface = Rmed
/// \param center_pos
/// \param mass
/// \param nr
/// \param np
/// \return vr, radial velocity with pressure support and potential smoothing
///
double get_vr_with_numerical_viscous_speed(const double r, const double mass){
	// num = 1/r d/dr (nu * Sigma * r^3 *  dw / dr)
	// den = d(r^2 w) / dr
	// vr = num / den

		// numerically compute 1/r d/dr (nu * Simga * r^3 * dw / dr)
	const double num = 1.0 / r * derive(r, mass, get_nu_S_r3_dwdr);

	// numerically compute Sigma * d(r^2 w) / dr
	const double Sigma = get_sigma(r);
	const double den = Sigma * derive(r, mass, get_r2_w);

	const double vr = num/den;

	return vr;

}

double get_vr_with_numerical_viscous_speed_wrapper(const t_radialarray &Rcenter, const t_radialarray &Rinterface, const Pair center_pos, const double mass, const int nr, const int np){


	const double r_cell = Rcenter[nr];

	const double phi = (double)np * dphi;
	const double cell_x = r_cell * std::cos(phi);
	const double cell_y = r_cell * std::sin(phi);

	const double x = cell_x - center_pos.x;
	const double y = cell_y - center_pos.y;
	const double r = std::sqrt(std::pow(x, 2) + std::pow(y, 2));

	const double vr = get_vr_with_numerical_viscous_speed(r, mass);

	return vr;
}

} // viscous_speed
