#include "viscous_radial_speed.h"

#include "../parameters.h"
#include "../Theo.h"
#include "../global.h"
#include "../util.h"
#include "../find_cell_id.h"
#include "../parameters.h"
#include "../boundary_conditions/boundary_conditions.h"
#include "../constants.h"


///
/// \brief get_nu, compute initial viscosity for the locally isothermal model
/// \param r
/// \param M
/// \return
///
[[maybe_unused]] static double get_nu(const double R, const double M){

	double sqrt_gamma;
	if(parameters::Adiabatic){
		sqrt_gamma = std::sqrt(parameters::ADIABATICINDEX);
	} else {
		sqrt_gamma = 1.0;
	}

	const double v_k = std::sqrt(constants::G * M / R);
	const double h =
	parameters::aspectratio_ref * std::pow(R, parameters::flaring_index);
	const double cs = sqrt_gamma * h * v_k;
	const double H = h * R;
	const double nu = parameters::viscous_alpha * cs * H;
	return nu;
}

/// does exact same as get_nu, expect it also ensures that
/// the soundspeed is inside the allowed Temperature range.
[[maybe_unused]] static double get_nu2(const double R, const double M, const double Sigma){

	const double v_k = std::sqrt(constants::G * M / R);
	const double h =
	parameters::aspectratio_ref * std::pow(R, parameters::flaring_index);

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
			parameters::minimum_temperature;
		const double energy_floor =
			temperature_floor *
			Sigma /
			parameters::MU * constants::R / (gamma - 1.0);

		const double temperature_ceil =
			parameters::maximum_temperature;
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

	const double nu = parameters::viscous_alpha * cs_adb * H;
	return nu;
}

static double get_sigma(const double R){
	const double S = parameters::sigma_slope;
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
    if(parameters::v_azimuthal_with_quadropole_support){
    return initial_locally_isothermal_smoothed_v_az_with_quadropole_moment(r, mass) / r;
    } else {
    return initial_locally_isothermal_smoothed_v_az(r, mass) / r;
    }
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
	// den = Sigma d(r^2 w) / dr
	// vr = num / den

		// numerically compute 1/r d/dr (nu * Simga * r^3 * dw / dr)
	const double num = 1.0 / r * derive(r, mass, get_nu_S_r3_dwdr);

	// numerically compute Sigma * d(r^2 w) / dr
	const double Sigma = get_sigma(r);
	const double den = Sigma * derive(r, mass, get_r2_w);

	const double vr = num/den;
	return vr;
}

/**
 * @brief get_vr_outer_viscous_speed_correction_factor
 * computes vr on the actual grid and compares it with get_vr_with_numerical_viscous_speed
 * @param r
 * @param mass
 * @return
 */
double get_vr_outer_viscous_speed_correction_factor(const double r, const double mass){
	// num = 1/r d/dr (nu * Sigma * r^3 *  dw / dr)
	// den = Sigma d(r^2 w) / dr
	// vr = num / den
	const unsigned int nr = std::max((unsigned int)2, IMIN + clamp_r_id_to_rmed_grid(get_rmed_id(r), true));

	const double rinf   = Radii[nr];
	const double rmed_p = GlobalRmed[nr+1];
	const double rmed   = GlobalRmed[nr];
	const double rmed_m = GlobalRmed[nr-1];
	const double rmed_m2= GlobalRmed[nr-2];

	// numerically compute 1/r d/dr (nu * Simga * r^3 * dw / dr)
	const double w_p = get_w(rmed_p, mass);
	const double w = get_w(rmed, mass);
	const double w_m = get_w(rmed_m, mass);
	const double w_m2 = get_w(rmed_m2, mass);

	const double dw_dr = (0.5*(w_p+w) - 0.5*(w+w_m))/(Radii[nr+1] - Radii[nr]);
	const double dw_dr_m = (0.5*(w+w_m) - 0.5*(w_m+w_m2))/(Radii[nr] - Radii[nr-1]);

	const double Sigma = get_sigma(rmed);
	const double nu = get_nu2(rmed, mass, Sigma);

	const double Sigma_m = get_sigma(rmed_m);
	const double nu_m = get_nu2(rmed_m, mass, Sigma_m);

	const double nu_S_r3_dwdr = nu*Sigma*std::pow(rmed, 3) * dw_dr;
	const double nu_S_r3_dwdr_m = nu_m*Sigma_m*std::pow(rmed_m, 3) * dw_dr_m;
	const double dr = rmed - rmed_m;
	const double num = 1.0 / rinf * (nu_S_r3_dwdr - nu_S_r3_dwdr_m)/dr;

	// numerically compute Sigma * d(r^2 w) / dr
	//const double Sigma_inf = 0.5*(Sigma + Sigma_m);
	const double Sigma_inf = Sigma_m;

	const double r2w   = std::pow(rmed  , 2) * get_w(rmed  , mass);
	const double r2w_m = std::pow(rmed_m, 2) * get_w(rmed_m, mass);

	const double dr2w_dr = (r2w - r2w_m)/(rmed - rmed_m);

	const double den = Sigma_inf * dr2w_dr;

	const double vr = num/den;
	const double vr_non_grid = get_vr_with_numerical_viscous_speed(rinf, mass);
	return vr/vr_non_grid;
}

/// Lookup table generation and reading
static const int N = 1000;
std::vector<double> r_table_outer(N);
std::vector<double> vr_table_outer(N);
static double deltaLogR_outer;
static double min_r_outer;
static double max_r_outer;

std::vector<double> r_table_inner(N);
std::vector<double> vr_table_inner(N);
static double deltaLogR_inner;
static double min_r_inner;
static double max_r_inner;

void init_vr_table_boundary(t_data &data){

	const auto & plsys = data.get_planetary_system();
	const unsigned int np = plsys.get_number_of_planets();
	const double mass = plsys.get_mass(np);

	//////////////// init outer arrays //////////////////////////////////////
	if(plsys.get_number_of_planets() == 2 && parameters::n_bodies_for_hydroframe_center == 1){
		// case for binary with frame on one star
		const auto &planet = plsys.get_planet(1);
		const double e = planet.get_eccentricity();
		const double a = planet.get_semi_major_axis();

		const double q = planet.get_mass() / (planet.get_mass() + plsys.get_planet(0).get_mass());
		const double center_of_mass_max_dist = a * (1.0 + e) * q;

		const double safety_dr = Rsup[NRadial-1] - Rinf[NRadial-1];

		min_r_outer = Rinf[NRadial-1] * boundary_conditions::damping_outer_limit - center_of_mass_max_dist  - safety_dr;
		max_r_outer = Rsup[NRadial-1] + center_of_mass_max_dist + safety_dr;
	} else {
		const double safety_dr = Rsup[NRadial-1] - Rinf[NRadial-1];
		min_r_outer = Rinf[NRadial-1] * boundary_conditions::damping_outer_limit - safety_dr;
		max_r_outer = Rsup[NRadial-1] + safety_dr;
	}

	min_r_outer = std::max(RMIN, min_r_outer);

	deltaLogR_outer = std::log10(max_r_outer / min_r_outer) / (double)N;

	for(unsigned int i = 0; i < N; ++i){
		const double r = min_r_outer * std::pow(10.0, (deltaLogR_outer * (double)i));
		const double vr = get_vr_with_numerical_viscous_speed(r, mass);
		const double corr = get_vr_outer_viscous_speed_correction_factor(r, mass);

		r_table_outer[i] = r;
		vr_table_outer[i] = vr*corr;
	}


	//////////////// init inner arrays //////////////////////////////////////
	const int Nr_limit = clamp_r_id_to_rmed_grid(
	get_rmed_id(RMIN * boundary_conditions::damping_inner_limit), false) + 1;
	const double safety_dr = Rsup[Nr_limit] - Rinf[Nr_limit];

	if(plsys.get_number_of_planets() == 2 && parameters::n_bodies_for_hydroframe_center == 1){
	// case for binary with frame on one star
	const auto &planet = plsys.get_planet(1);
	const double e = planet.get_eccentricity();
	const double a = planet.get_semi_major_axis();

	const double q = planet.get_mass() / (planet.get_mass() + plsys.get_planet(0).get_mass());
	const double center_of_mass_max_dist = a * (1.0 + e) * q;

	min_r_inner = Rinf[1] - center_of_mass_max_dist  - safety_dr;
	min_r_inner = std::max(min_r_inner, Rsup[0] - Rinf[0]);
	max_r_inner = Rsup[1] * boundary_conditions::damping_inner_limit + center_of_mass_max_dist + safety_dr;
} else {
	min_r_inner = Rinf[1] - safety_dr;
	max_r_inner = Rsup[1] * boundary_conditions::damping_inner_limit + safety_dr;
}

	deltaLogR_inner = std::log10(max_r_inner / min_r_inner) / (double)N;

	const unsigned int np2 = parameters::n_bodies_for_hydroframe_center;
	const double mass2 = plsys.get_mass(np2);
for(unsigned int i = 0; i < N; ++i){
	const double r = min_r_inner * std::pow(10.0, (deltaLogR_inner * (double)i));
	const double vr = get_vr_with_numerical_viscous_speed(r, mass2);

	r_table_inner[i] = r;
	vr_table_inner[i] = vr;
}
}

/**
 * @brief lookup_initial_vr_outer, note that each mpi thread computes its own domain inside its bounds
 * thus it is not bit-wise exact when changing the number of mpi cores
 * @param r
 * @return vr at distance r from center of mass
 */
double lookup_initial_vr_outer(const double r){
	const int i = int(std::floor(std::log10(r / min_r_outer) / deltaLogR_outer));

	if (i > N - 1) {
	printf("Lookup outer: Cell out of bounds outwards %.5e %.5e\n", r, max_r_outer);
	return vr_table_outer[N-1];
	}
	if (i < 0) {
	printf("Lookup outer: Cell out of bounds inwards %.5e %.5e\n", r, min_r_outer);
	return vr_table_outer[0];
	}

	const double r_i = r_table_outer[i];
	const double r_i1 = r_table_outer[i+1];

	const double vr_i = vr_table_outer[i];
	const double vr_i1 = vr_table_outer[i+1];

	const double x = (r - r_i) / (r_i1 - r_i);

	// linear interpolation
	const double vr = (1.0 - x) * vr_i + x * vr_i1;

	return vr;
}
/**
 * @brief lookup_initial_vr_inner, note that each mpi thread computes its own domain inside its bounds
 * thus it is not bit-wise exact when changing the number of mpi cores
 * @param r
 * @return vr at distance r from center of mass
 */
double lookup_initial_vr_inner(const double r){
	const int i = int(std::floor(std::log10(r / min_r_inner) / deltaLogR_inner));

	if (i > N - 1) {
	printf("Lookup inner: Cell out of bounds outwards %.5e %.5e\n", r, max_r_inner);
	return vr_table_inner[N-1];
	}
	if (i < 0) {
	printf("Lookup inner: Cell out of bounds inwards %.5e %.5e\n", r, min_r_inner);
	return vr_table_inner[0];
	}

	const double r_i = r_table_inner[i];
	const double r_i1 = r_table_inner[i+1];

	const double vr_i = vr_table_inner[i];
	const double vr_i1 = vr_table_inner[i+1];

	const double x = (r - r_i) / (r_i1 - r_i);

	// linear interpolation
	const double vr = (1.0 - x) * vr_i + x * vr_i1;

	return vr;
}

} // viscous_speed
