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

	(void)get_nu; // to prevent unused function warning

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

static const int N = 1000;
std::vector<double> r_table(N);
std::vector<double> vr_table(N);
static double deltaLogR;
static double min_r;
static double max_r;

void init_vr_table_outer_boundary(t_data &data){

	const auto & plsys = data.get_planetary_system();
	const unsigned int np = plsys.get_number_of_planets();
	const double mass = plsys.get_mass(np);

	if(plsys.get_number_of_planets() == 2 && parameters::n_bodies_for_hydroframe_center == 1){
		// case for binary with frame on one star
		const auto &planet = plsys.get_planet(1);
		const double e = planet.get_eccentricity();
		const double a = planet.get_semi_major_axis();

		const double q = planet.get_mass() / (planet.get_mass() + plsys.get_planet(0).get_mass());
		const double center_of_mass_max_dist = a * (1.0 + e) * q;

		min_r = Rinf[NRadial-1] * parameters::damping_outer_limit - center_of_mass_max_dist;
		max_r = Rsup[NRadial-1] + center_of_mass_max_dist;
	} else {
		min_r = Rinf[NRadial-1] * parameters::damping_outer_limit;
		max_r = Rsup[NRadial-1];
	}

	min_r = std::max(RMIN, min_r);

	deltaLogR = std::log10(max_r / min_r) / (double)N;

	for(unsigned int i = 0; i < N; ++i){
		const double r = min_r * std::pow(10.0, (deltaLogR * (double)i));
		const double vr = get_vr_with_numerical_viscous_speed(r, mass);

		r_table[i] = r;
		vr_table[i] = vr;
	}
}

double lookup_initial_vr(const double r){
	const int i = int(std::floor(std::log10(r / min_r) / deltaLogR));

	if (i > N - 1) {
	printf("Cell out of bounds outwards %.5e %.5e\n", r, max_r);
	return vr_table[N-1];
	}
	if (i < 0) {
	printf("Cell out of bounds inwards\n");
	return vr_table[0];
	}

	const double r_i = r_table[i];
	const double r_i1 = r_table[i+1];

	const double vr_i = vr_table[i];
	const double vr_i1 = vr_table[i+1];

	const double x = (r - r_i) / (r_i1 - r_i);

	// linear interpolation
	const double vr = (1.0 - x) * vr_i + x * vr_i1;

	return vr;
}

/// \brief derive_5th_order, see https://en.wikipedia.org/wiki/Numerical_differentiation
/// \param r, position variable
/// \param mass, parameter
/// \param f, function to derive
/// \return df / dr at r
static double derive(const double r, const double mass, double (*f)(double, double)){

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
static double get_w(const double r, const double mass){
	return initial_locally_isothermal_smoothed_v_az(r, mass) / r;
}

///
/// \brief get_r2_w
/// \return r^2 * omega
///
static double get_r2_w(const double r, const double mass){
	const double omega = get_w(r,mass);
	return std::pow(r, 2) * omega;
}

///
/// \brief get_nu_S_r3_dwdr
/// \param r
/// \param mass
/// \return nu * Sigma * r^3 * dw/dr
///
static double get_nu_S_r3_dwdr(const double r, const double mass){

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

} // viscous_speed
