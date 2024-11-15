/**
	\file SourceEuler.c

	Contains routines used by the hydrodynamical loop. More specifically, it
   contains the main loop itself and all the source term substeps (with the
   exception of the evaluation of the viscous force). The transport substep is
   treated elsewhere.
*/
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cassert>
#include <cfloat>
#include <climits>


#include <vector>

#include "LowTasks.h"
#include "SideEuler.h"
#include "SourceEuler.h"
#include "Theo.h"
#include "TransportEuler.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "opacity.h"
#include "parameters.h"
#include "pvte_law.h"
#include "quantities.h"
#include "selfgravity.h"
#include "units.h"
#include "util.h"
#include "viscosity/viscosity.h"
#include "simulation.h"
#include "frame_of_reference.h"
#include "simulation.h"
#include "fld.h"
#include "compute.h"

#include <cstring>

/**
	Checks polargrid for negative entries.

	\param array polargrid to check
	\returns >0 if there are negative entries, 0 otherwise
*/
// static int DetectCrash(t_polargrid *array)
// {
//     unsigned int result = 0;

//     for (unsigned int n_radial = 0; n_radial < array->Nrad; ++n_radial) {
// 	for (unsigned int n_azimuthal = 0; n_azimuthal < array->Nsec;
// 	     ++n_azimuthal) {
// 	    /// since nan < 0 is false and nan > 0 is false
// 	    /// we need to assure that array > 0 to catch bad values
// 	    if (!((*array)(n_radial, n_azimuthal) > 0.0)) {
// 		logging::print(LOG_WARNING "%s negative in cell: (%u,%u)=%g\n",
// 			       array->get_name(), n_radial, n_azimuthal,
// 			       (*array)(n_radial, n_azimuthal));
// 		result += 1;
// 	    }
// 	}
//     }

//     return result;
// }

// static void HandleCrash(t_data &data)
// {
//     if (DetectCrash(&data[t_data::SIGMA])) {
// 	logging::print(LOG_ERROR "DetectCrash: Density < 0\n");
// 	PersonalExit(1);
//     }

//     if (parameters::Adiabatic) {
// 	if (DetectCrash(&data[t_data::ENERGY])) {
// 	    logging::print(LOG_ERROR "DetectCrash: Energy < 0\n");
// 	    PersonalExit(1);
// 	}
//     }
// }

void SetTemperatureFloorCeilValues(t_data &data, std::string filename, int line)
{
    if (assure_temperature_range(data)) {
	logging::print(
	    LOG_DEBUG
	    "Found temperature outside the valid range of %g to %g %s in %s: %d.\n",
	    parameters::minimum_temperature, parameters::maximum_temperature,
	    units::temperature.get_cgs_symbol(), filename.c_str(), line);
    }
}

/**
	Assures miminum value in each cell.

	\param dst polar grid
	\param minimum_value minimum value
*/
bool assure_minimum_value(t_polargrid &dst, double minimum_value)
{
    bool found = false;
    const bool is_dens = strcmp(dst.get_name(), "Sigma") == 0;

	const unsigned int Nr = dst.get_size_radial();
	const unsigned int Nphi = dst.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		if (dst(nr, naz) < minimum_value) {
		if (is_dens) {
		    double mass_delta =
			(minimum_value - dst(nr, naz)) *
			Surf[nr];
		    sum_without_ghost_cells(MassDelta.FloorMassCreation, mass_delta,
						nr);
		}
		dst(nr, naz) = minimum_value;
#ifndef NDEBUG
		logging::print(LOG_DEBUG
			       "assure_minimum_value: %s(%u,%u)=%g < %g\n",
				   dst.get_name(), nr, naz,
				   dst(nr, naz), minimum_value);
#endif
		found = true;
	    }
	}
    }

    return found;
}

bool assure_temperature_range(t_data &data)
{
    bool found = false;

    t_polargrid &energy = data[t_data::ENERGY];
    t_polargrid &density = data[t_data::SIGMA];

    const double Tmin = parameters::minimum_temperature;

    const double Tmax = parameters::maximum_temperature;

	const unsigned int Nr = energy.get_size_radial();
	const unsigned int Nphi = energy.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {

		const double mu = pvte::get_mu(data, nr, naz);
		const double gamma_eff = pvte::get_gamma_eff(data, nr, naz);

	    const double minimum_energy = Tmin *
					  density(nr, naz) / mu *
					  constants::R / (gamma_eff - 1.0);

	    const double maximum_energy = Tmax *
					  density(nr, naz) / mu *
					  constants::R / (gamma_eff - 1.0);

		if (!(energy(nr, naz) > minimum_energy)) {
#ifndef NDEBUG
		logging::print(
		    LOG_DEBUG "assure_minimum_temperature: (%u,%u)=%g<%g\n",
			nr, naz,
			energy(nr, naz) *
			units::temperature.get_code_to_cgs_factor() /
			density(nr, naz) * mu / constants::R *
			(gamma_eff - 1.0),
		    Tmin * units::temperature.get_code_to_cgs_factor(), Tmin);
#endif
		energy(nr, naz) =
			Tmin * density(nr, naz) / mu * constants::R /
		    (gamma_eff - 1.0);
		found = true;
	    }

		if (!(energy(nr, naz) < maximum_energy)) {
#ifndef NDEBUG
		logging::print(
		    LOG_DEBUG "assure_maximum_temperature: (%u,%u)=%g>%g\n",
			nr, naz,
			energy(nr, naz) *
			units::temperature.get_code_to_cgs_factor() /
			density(nr, naz) * mu / constants::R *
			(gamma_eff - 1.0),
		    Tmax * units::temperature.get_code_to_cgs_factor(), Tmax);
#endif
		energy(nr, naz) =
			Tmax * density(nr, naz) / mu * constants::R /
		    (gamma_eff - 1.0);
		found = true;
	    }
	}
    }

    return found;
}


void recalculate_viscosity(t_data &data, const double current_time)
{

	if (parameters::Locally_Isothermal) {
	if (parameters::aspectratio_mode > 0) {
		compute_sound_speed(data, current_time);
		compute_scale_height(data, current_time);
	}
	}
	if (parameters::Adiabatic || parameters::Polytropic) {
	if (parameters::variableGamma) {
		pvte::compute_gamma_mu(data);
	}
	compute_sound_speed(data, current_time);
	compute_scale_height(data, current_time);
	}

	viscosity::update_viscosity(data);
}

void recalculate_derived_disk_quantities(t_data &data, const double current_time)
{

    if (parameters::Locally_Isothermal) {
	if (parameters::aspectratio_mode > 0) {
		compute_sound_speed(data, current_time);
		compute_pressure(data);
		compute_temperature(data);
		compute_scale_height(data, current_time);
	} else {
		compute_pressure(data);
	}
    }
    if (parameters::Adiabatic || parameters::Polytropic) {
	if (parameters::variableGamma) {
	    pvte::compute_gamma_mu(data);
	}
	compute_temperature(data);
	compute_sound_speed(data, current_time);
	compute_scale_height(data, current_time);
	compute_pressure(data);
    }

    viscosity::update_viscosity(data);
}

void init_euler(t_data &data, const double current_time)
{
    InitCellCenterCoordinates();
    InitTransport();

	static const unsigned int N_planets =
	data.get_planetary_system().get_number_of_planets();
	g_xpl.resize(N_planets);
	g_ypl.resize(N_planets);
	g_mpl.resize(N_planets);
	g_rpl.resize(N_planets);
	g_cubic_smoothing_radius.resize(N_planets);

    if (parameters::Locally_Isothermal) {
	compute_sound_speed(data, current_time);
	compute_pressure(data);
	compute_temperature(data);
	compute_scale_height(data, current_time);
    }

    if (parameters::Adiabatic || parameters::Polytropic) {
	if (parameters::variableGamma) {
		compute_sound_speed(data, current_time);
		compute_scale_height(data, current_time);
	    pvte::compute_gamma_mu(data);
	}
	compute_temperature(data);
	compute_sound_speed(data, current_time);
	compute_scale_height(data, current_time);
	compute_pressure(data);
    }

    viscosity::update_viscosity(data);
	compute_heating_cooling_for_CFL(data, current_time);
}


void FreeEuler()
{
    FreeTransport();
    FreeCellCenterCoordinates();
}

/**
	copy one polar grid into another

	\param dst destination polar grid
	\param src source polar grid
*/
void copy_polargrid(t_polargrid &dst, const t_polargrid &src)
{
    assert((dst.get_size_radial() == src.get_size_radial()) &&
	   (dst.get_size_azimuthal() == src.get_size_azimuthal()));

    std::memcpy(dst.Field, src.Field,
		dst.get_size_radial() * dst.get_size_azimuthal() *
		    sizeof(*dst.Field));
}

/**
	switch one polar grid with another

	\param dst destination polar grid
	\param src source polar grid
	switches polar grids
*/
void move_polargrid(t_polargrid &dst, t_polargrid &src)
{
    assert((dst.get_size_radial() == src.get_size_radial()) &&
	   (dst.get_size_azimuthal() == src.get_size_azimuthal()));

    std::swap(dst.Field, src.Field);
}

static void momentum_update_radial(t_data &data, const double dt) {

	const unsigned int Nphi = data[t_data::ENERGY].get_size_azimuthal();

	t_polargrid &Sigma = data[t_data::SIGMA];
	t_polargrid &P = data[t_data::PRESSURE];
	t_polargrid &Phi = data[t_data::POTENTIAL];
	t_polargrid &vaz = data[t_data::V_AZIMUTHAL];
	t_polargrid &vrad = data[t_data::V_RADIAL];


    // update v_radial with source terms
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = One_no_ghost_vr; nr < MaxMo_no_ghost_vr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
	    // 1/Sigma * dP/dr : Sigma is calculated as a mean value between the
	    // neightbour cells
	    double gradp = 2.0 / (Sigma(nr, naz) + Sigma(nr - 1, naz));
		gradp *= (P(nr, naz) - P(nr - 1, naz));
		gradp *= InvDiffRmed[nr];

	    // dPhi/dr
	    double gradphi;
	    if (parameters::body_force_from_potential) {
			gradphi = (Phi(nr, naz) - Phi(nr-1, naz)) * InvDiffRmed[nr];
	    } else {
			gradphi = - (data[t_data::ACCEL_RADIAL](nr, naz)
				    + data[t_data::ACCEL_RADIAL](nr-1, naz)) * 0.5;
	    }

		const unsigned int naz_next = (naz == Nphi-1 ? 0 : naz + 1);
	    // v_phi^2/r : v_phi^2 is calculated by a mean in both directions
		const double vsum = vaz(nr, naz) +
			vaz(nr, naz_next) +
			vaz(nr - 1, naz) +
			vaz(nr - 1, naz_next);
		const double vt = 0.25 * vsum + Rinf[nr] * refframe::OmegaFrame;
	    const double vt2 = vt * vt;

        const double centrifugal_accel = vt2 * InvRinf[nr];

	    // add all terms to new v_radial: v_radial_new = v_radial +
	    // dt*(source terms)
		vrad(nr, naz) += dt * (-gradp - gradphi + centrifugal_accel);
	}
    }

}


static void momentum_update_azimuthal(t_data &data, const double dt) {

	const unsigned int Nphi = data[t_data::ENERGY].get_size_azimuthal();

	double supp_torque = 0.0; // for imposed disk drift

    // update v_azimuthal with source terms
	#pragma omp parallel for
	for (unsigned int nr = Zero_no_ghost; nr < Max_no_ghost; ++nr) {

	if (parameters::IMPOSEDDISKDRIFT != 0.0) {
		supp_torque = parameters::IMPOSEDDISKDRIFT * 0.5 *
			  std::pow(Rmed[nr], -2.5 + parameters::sigma_slope);
	}
	//const double invdxtheta = 1.0 / (dphi * Rmed[n_radial]);
	const double invdxtheta = 2.0 / (dphi * (Rsup[nr]+Rinf[nr]));

	for (unsigned int naz = 0; naz < Nphi; ++naz) {

		const double naz_prev = (naz == 0 ? Nphi-1 : naz - 1);
	    // 1/Sigma 1/r dP/dphi
	    const double gradp =
		2.0 /
		(data[t_data::SIGMA](nr, naz) +
		 data[t_data::SIGMA](nr, naz_prev)) *
		(data[t_data::PRESSURE](nr, naz) -
		 data[t_data::PRESSURE](nr, naz_prev)) *
		invdxtheta;

	    // 1/r dPhi/dphi
	    double gradphi;
	    if (parameters::body_force_from_potential) {
		gradphi = (data[t_data::POTENTIAL](nr, naz) -
			   data[t_data::POTENTIAL](nr, naz_prev)) *
			  invdxtheta;
	    } else {
		gradphi = -(data[t_data::ACCEL_AZIMUTHAL](nr, naz) +
			   data[t_data::ACCEL_AZIMUTHAL](nr, naz_prev)) * 0.5;
	    }

	    // add all terms to new v_azimuthal: v_azimuthal_new = v_azimuthal +
	    // dt*(source terms)
		data[t_data::V_AZIMUTHAL](nr, naz) =
		data[t_data::V_AZIMUTHAL](nr, naz) +
		dt * (-gradp - gradphi);

	    if (parameters::IMPOSEDDISKDRIFT != 0.0) {
		// add term for imposed disk drift
		data[t_data::V_AZIMUTHAL](nr, naz) +=
		    dt * supp_torque;
	    }
	}
    }
}

/**
	In this substep we take into account the source part of Euler equations.
   We evolve velocities with pressure gradients, gravitational forces and
   curvature terms
*/
void update_with_sourceterms(t_data &data, const double dt)
{

	// We do Self gravity first, so we don't have to recompute aspect ratio
	if (parameters::self_gravity) {
	selfgravity::compute(data, dt, true);
	}

	// Momentum update due to 
	momentum_update_radial(data, dt);
	momentum_update_azimuthal(data, dt);
	compression_heating(data, dt);

	if(ECC_GROWTH_MONITOR){
		quantities::calculate_disk_delta_ecc_peri(data, delta_ecc_source, delta_peri_source);
	}

}


/**
 * Perform energy update due to compression
 * i.e. due to P div v
*/
void compression_heating(t_data &data, const double dt){
	
	const unsigned int Nr = data[t_data::QMINUS].get_max_radial(); // = Size - 1
	const unsigned int Nphi = data[t_data::QMINUS].get_size_azimuthal();

	// Calculate the energy update due to div_v
	if (parameters::Adiabatic) {
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const unsigned int naz_next = (naz == Nphi-1 ? 0 : naz + 1);
		// div(v) = 1/r d(r*v_r)/dr + 1/r d(v_phi)/dphi
		const double DIV_V =
			(data[t_data::V_RADIAL](nr+1, naz)*Ra[nr+1] -
				data[t_data::V_RADIAL](nr, naz)*Ra[nr]) *
			InvDiffRsupRb[nr] +
			(data[t_data::V_AZIMUTHAL](nr, naz_next) -
			 data[t_data::V_AZIMUTHAL](nr, naz)) *
			invdphi * InvRb[nr];

		const double gamma =
			pvte::get_gamma_eff(data, nr, naz);


		// Like D'Angelo et al. 2003 eq. 24
		const double energy_old =
			data[t_data::ENERGY](nr, naz);
		const double energy_new =
		    energy_old * std::exp(-(gamma - 1.0) * dt * DIV_V);
		data[t_data::ENERGY](nr, naz) = energy_new;
		
	    } // end azimuthal loop
	} // end radial loop
    } // end if adiabatic
}


static void viscous_heating(t_data &data) {
	const unsigned int Nr_m1 = data[t_data::QPLUS].get_max_radial();
	const unsigned int Nphi = data[t_data::QPLUS].get_size_azimuthal();

	/* We calculate the heating source term Qplus from i=1 to max-1 */
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nr_m1; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {

		if (data[t_data::VISCOSITY](nr, naz) != 0.0) {
			const unsigned int naz_next = (naz == Nphi-1 ? 0 : naz + 1);
		    // average tau_r_phi over 4 cells
            const double tau_r_phi =
			0.25 *
			(data[t_data::TAU_R_PHI](nr, naz) +
			 data[t_data::TAU_R_PHI](nr + 1, naz) +
			 data[t_data::TAU_R_PHI](nr, naz_next) +
			 data[t_data::TAU_R_PHI](nr + 1, naz_next));

		    double qplus =
			1.0 /
			(2.0 * data[t_data::VISCOSITY](nr, naz) *
			 data[t_data::SIGMA](nr, naz)) *
			(std::pow(data[t_data::TAU_R_R](nr, naz),
				  2) +
			 2 * std::pow(tau_r_phi, 2) +
			 std::pow(
				 data[t_data::TAU_PHI_PHI](nr, naz),
			     2));
		    qplus +=
			(2.0 / 9.0) *
			data[t_data::VISCOSITY](nr, naz) *
			data[t_data::SIGMA](nr, naz) *
			std::pow(data[t_data::DIV_V](nr, naz), 2);

		    qplus *= parameters::heating_viscous_factor;
			data[t_data::QPLUS](nr, naz) += qplus;
		}
	    }
	}
}

static void irradiation_single(t_data &data, const t_planet &planet) {

	const double rampup_time = planet.get_irradiation_rampuptime();

	double ramping = 1.0;
	
	if (sim::time < rampup_time) {
	    ramping =
		1.0 -
		std::pow(std::cos(sim::time * M_PI / 2.0 / rampup_time), 2);
	}

	const double x = planet.get_x();
	const double y = planet.get_y();
	const double R_star = planet.get_planet_radial_extend();
	const double T_star = planet.get_temperature();

	const double l1 = planet.get_dimensionless_roche_radius() *
			planet.get_semi_major_axis() * (1.0 - planet.get_eccentricity());
	double min_dist;
	if(x*x + y*y > 1e-10){
	    min_dist = std::max(R_star, l1 * planet.get_cubic_smoothing_factor());
	} else {
		min_dist = R_star;
	}

	const unsigned int Nrad = data[t_data::QPLUS].get_max_radial();
	const unsigned int Naz = data[t_data::QPLUS].get_max_azimuthal();
	// Star heating from Menou & Goodman 2004
	// Implementation follows D'Angelo & Marzari 2012 (doi:10.1088/0004-637X/757/1/50)

	#pragma omp parallel for collapse(2)
	for (unsigned int nrad = 1; nrad < Nrad; ++nrad) {
	for (unsigned int naz = 0; naz <= Naz; ++naz) {
		const unsigned int ncell = nrad * data[t_data::SIGMA].get_size_azimuthal() + naz;
		const double xc = CellCenterX->Field[ncell];
		const double yc = CellCenterY->Field[ncell];
		const double distance_measured = std::sqrt(std::pow(x - xc, 2) + std::pow(y - yc, 2));
		const double distance = std::max(distance_measured, min_dist);
		const double roverd = distance < R_star ? 1.0 : R_star/distance;

		const double HoverR = data[t_data::ASPECTRATIO](nrad, naz);
		const double sigma = constants::sigma.get_code_value();
		const double tau_eff = data[t_data::TAU_EFF](nrad, naz);
		const double eps = 0.5;
		// choose according to Chiang & Goldreich (1997)
		const double dlogH_dlogr = 9.0 / 7.0;
		
		// irradiation contribution near and far from the star
		const double W_G = 0.4 * roverd + HoverR * (dlogH_dlogr - 1.0);

		const double T_irrad_pow4 = (1.0 - eps) * std::pow(T_star, 4) * std::pow(roverd, 2) * W_G;

		const double qplus = 2.0 * sigma * T_irrad_pow4 / tau_eff;

		data[t_data::QPLUS](nrad, naz) += ramping * qplus;
	}
	}

}


static void irradiation(t_data &data) {

	const auto & plsys = data.get_planetary_system();
	const unsigned int Npl = plsys.get_number_of_planets();

	for (unsigned int npl=0; npl < Npl; npl ++) { 
		const auto& planet = plsys.get_planet(npl);
		if (planet.get_irradiate()) {
			irradiation_single(data, planet);
		}
	}
}


void calculate_qplus(t_data &data)
{
    data[t_data::QPLUS].clear();

    if (parameters::heating_viscous_enabled) {
			viscous_heating(data);
    }
    if (parameters::heating_star_enabled) {
		if (!(parameters::cooling_surface_enabled || parameters::cooling_scurve_enabled)) {
			compute::midplane_density(data, sim::time);
		    compute::kappa_eff(data);
		}
		irradiation(data);
	}    
}

/* Perform thermal relaxation also called beta cooling.
*/
static void thermal_relaxation(t_data &data, const double current_time) {

	t_polargrid &Qminus = data[t_data::QMINUS];
	const unsigned int Nr = Qminus.get_max_radial(); // = Size - 1
	const unsigned int Nphi = Qminus.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nr; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {
		// Q- = E Omega/beta
		const double r = Rmed[nr];
		const double omega_k = calculate_omega_kepler(r);
		const double E = data[t_data::ENERGY](nr, naz);
		const double t_ramp_up = parameters::cooling_beta_ramp_up;

		double beta_inv = 1 / parameters::cooling_beta;
		if (t_ramp_up > 0.0) {
			const double t = current_time;
		    double ramp_factor =
			1 - std::exp(-std::pow(2 * t / t_ramp_up, 2));
		    beta_inv = beta_inv * ramp_factor;
		}

		double delta_E = E;
		if (parameters::cooling_beta_reference) {
		    const double sigma =
			data[t_data::SIGMA](nr, naz);
		    const double sigma0 =
			data[t_data::SIGMA0](nr, naz);
		    const double E0 =
			data[t_data::ENERGY0](nr, naz);
		    delta_E -= E0 / sigma0 * sigma;
		}
		if (parameters::cooling_beta_model) {
		    const double sigma =
			data[t_data::SIGMA](nr, naz);
		    const double E0 =
			1.0 / (parameters::ADIABATICINDEX - 1.0) *
			std::pow(parameters::aspectratio_ref, 2) *
			std::pow(Rmed[nr], 2.0 * parameters::flaring_index - 1.0) *
			constants::G * hydro_center_mass * sigma;
		    delta_E -= E0;
		}
		if(parameters::cooling_beta_floor){
		    const double Tmin = parameters::minimum_temperature;
		    const double gamma_eff = pvte::get_gamma_eff(data, nr, naz);
		    const double mu = pvte::get_mu(data, nr, naz);
		    const double minimum_energy = Tmin *
						  data[t_data::SIGMA](nr, naz) / mu *
						  constants::R / (gamma_eff - 1.0);
		    delta_E -= minimum_energy;
		}
		const double qminus = delta_E * omega_k * beta_inv;


		Qminus(nr, naz) += qminus;
	    }
	}
}


static void thermal_cooling(t_data &data) {

	t_polargrid & Qminus = data[t_data::QMINUS];
	const unsigned int Nr = Qminus.get_max_radial(); // = Size - 1
	const unsigned int Nphi = Qminus.get_size_azimuthal();

	compute::midplane_density(data, sim::time);
	compute::kappa_eff(data);

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nr; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {

		const double temperature = data[t_data::TEMPERATURE](nr, naz);

		// Q = factor 2 sigma_sb T^4 / tau_eff
		const double factor = parameters::surface_cooling_factor;
		const double sigma_sb = constants::sigma.get_code_value();
		const double T4 = std::pow(temperature, 4);
		const double tau_eff =
			data[t_data::TAU_EFF](nr, naz);
		const double Tmin4 =
			std::pow(parameters::minimum_temperature, 4);

		const double qminus =
		    factor * 2 * sigma_sb * (T4 - Tmin4) / tau_eff;

		Qminus(nr, naz) += qminus;
	    }
	}
}


static void scurve_cooling(t_data &data) {

	const double SigmaCGS_threshold = 2.0;
	const double temperatureCGS_threshold = 1200.0;

	double muExponent;
	double F_hot_const;
	if(parameters::cooling_scurve_type){
	    F_hot_const = 23.405; // Kimura et al. 2020 (https://doi.org/10.1093/pasj/psz144)
	    muExponent = 0.31; // Sign change in Kimura
	} else {
	    F_hot_const = 25.49; // Ichikawa & Osaki 1992 (https://ui.adsabs.harvard.edu/abs/1992PASJ...44...15I/abstract)
	    muExponent = -0.31;
	}

	t_polargrid & Qminus = data[t_data::QMINUS];
	const unsigned int Nr = Qminus.get_max_radial(); // = Size - 1
	const unsigned int Nphi = Qminus.get_size_azimuthal();

    /// Scruve cooling according to Ichikawa & Osaki (1992)
    /// See page 21 & 22 in https://articles.adsabs.harvard.edu/full/1992PASJ...44...15I
    #pragma omp parallel for collapse(2)
    for (unsigned int nr = 1; nr < Nr; ++nr) {
        for (unsigned int naz = 0; naz < Nphi; ++naz) {

	const double Sigma = data[t_data::SIGMA](nr, naz);

	const double SigmaCGS = Sigma*units::surface_density.get_code_to_cgs_factor();
	const double SigmaCGS_tmp = std::max(SigmaCGS, SigmaCGS_threshold);

	const double temperatureCGS = data[t_data::TEMPERATURE](nr, naz) * units::temperature.get_code_to_cgs_factor();
	const double temperatureCGS_tmp = std::max(temperatureCGS, temperatureCGS_threshold);

	const double rCGS = Rmed[nr] * units::length.get_code_to_cgs_factor();

	const double mu = pvte::get_mu(data, nr, naz);

	const double M = hydro_center_mass*units::mass.get_code_to_cgs_factor();

	const double cgs_G = constants::global_cgs_G;

	const double omega_keplerCGS = std::sqrt(cgs_G * M / (rCGS * rCGS * rCGS));

	const double sigma_sb_cgs = constants::sigma.get_cgs_value();

	// Solve sigma T^4 = F_cool(T) for T
	const double logTA = -1.0/5.49 * (0.62 * std::log10(omega_keplerCGS)
					    + 1.62 * std::log10(SigmaCGS_tmp)
					    + muExponent * std::log10(mu) - 25.48
					    - std::log10(sigma_sb_cgs));
	const double TA = std::pow(10.0, logTA);

	const double FA = sigma_sb_cgs * std::pow(TA, 4);
	const double logFA = std::log10(FA);

	const double KCGS = 11.0 + 0.4 * std::log10(2.0e10/rCGS);
	const double logFB = std::max(KCGS, logFA);

	// Solve F_hot(T) = logFB for T
	const double logTB_aux = std::log10(omega_keplerCGS)+ 2.0 * std::log10(SigmaCGS_tmp) +
				 0.5 * std::log10(mu) + F_hot_const;
	const double logTB = (logFB + logTB_aux) / 8.0;
	const double TB = std::pow(10.0, logTB);

	double logFtot;

	if (temperatureCGS_tmp < TA) {
	// F_cold
	    logFtot = 9.49 * std::log10(temperatureCGS_tmp) + 0.62 *
							      std::log10(omega_keplerCGS) + 1.62 * std::log10(SigmaCGS_tmp) +
		      muExponent * std::log10(mu) - 25.48;
	} else if (temperatureCGS_tmp > TB) {
	// F_hot
	    logFtot = 8.0 * std::log10(temperatureCGS_tmp) -
		      std::log10(omega_keplerCGS) - 2.0 * std::log10(SigmaCGS_tmp) -
		      0.5 * std::log10(mu) - F_hot_const;
	} else {
	// F_intermediate
	    logFtot = (logFA - logFB) * std::log10(temperatureCGS_tmp / TB)
			  / std::log10(TA / TB) + logFB;
    }

	const double T4 = std::pow(data[t_data::TEMPERATURE](nr, naz), 4);
	const double sigma_sb = constants::sigma.get_code_value();
	const double factor = parameters::surface_cooling_factor;

	/// Limiting cooling to black body radiation is not described in the equations
	/// of either Ichikawa 1992 or Kimura 2020, but it is shown
	/// in Fig.3 of Kimura et. al 2020
	/// without it, the cooling by Kimura has a jump from cold to hot branch
	/// at densities below ~0.7 g/cm^2 that reduces the timestep
	double F_tot = std::pow(10.0, logFtot)*units::energy_flux.get_cgs_to_code_factor();
	// Scale with power law below theshold Temperature or Surface density
	F_tot *= std::pow((SigmaCGS / SigmaCGS_tmp), 0.5);
	F_tot *= std::pow(temperatureCGS / temperatureCGS_tmp, 2);

	const double F_Blackbody = sigma_sb * T4;
	const double qminus_scurve = 2.0 * factor * std::min(F_tot, F_Blackbody);

	Qminus(nr, naz) += qminus_scurve;

	const double tau_eff =  factor * 2 * sigma_sb * T4 / qminus_scurve;
	data[t_data::TAU_EFF](nr, naz) = tau_eff;
	} // azimuthal loop
    } // radial loop
}


void calculate_qminus(t_data &data, const double current_time)
{
    // clear up all Qminus terms
    data[t_data::QMINUS].clear();

    // beta cooling
    if (parameters::cooling_beta_enabled) {
		thermal_relaxation(data, current_time);
    }

    // local radiative cooling
    if (parameters::cooling_surface_enabled) {
		thermal_cooling(data);
    }

    //S-curve cooling
    if(parameters::cooling_scurve_enabled){
		scurve_cooling(data);
    }
}

/**
	In this substep we take into account the source part of energy equation.
   We evolve internal energy with compression/dilatation and heating terms
*/
void SubStep3(t_data &data, const double current_time, const double dt)
{
	compute_temperature(data);
	calculate_qminus(data, current_time); // first to calculate teff
	calculate_qplus(data);

	const unsigned int Nr = data[t_data::TAU_COOL].get_size_radial();
	const unsigned int Nphi = data[t_data::TAU_COOL].get_size_azimuthal();

    // calculate tau_cool if needed for output
    if (data[t_data::TAU_COOL].get_write_1D() ||
	data[t_data::TAU_COOL].get_write_2D()) {
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nr-1; ++nr) { // Q_plus / Q_minus are not computed for nr = 1 or nr = Nr-1 (check max_radial usage there)
		for (unsigned int naz = 0; naz < Nphi; ++naz) {
		data[t_data::TAU_COOL](nr, naz) =
			data[t_data::ENERGY](nr, naz) /
			data[t_data::QMINUS](nr, naz);
	    }
	}
    }

    // calculate pDV for write out
    if (data[t_data::P_DIVV].get_write_1D() ||
	data[t_data::P_DIVV].get_write_2D() ||
	fld::radiative_diffusion_enabled) {

	data.pdivv_total = 0;
	double &pdivv_tmp = data.pdivv_total;

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nr-1; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {
		double pdivv =
			(pvte::get_gamma_eff(data, nr, naz) - 1.0) *
			dt * data[t_data::DIV_V](nr, naz) *
			data[t_data::ENERGY](nr, naz);
		data[t_data::P_DIVV](nr, naz) = pdivv;

		sum_without_ghost_cells(pdivv_tmp, pdivv, nr);
	    }
	}
    }

    // Now we can update energy with source terms
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nr-1; ++nr) {// Don't update ghost cells
	for (unsigned int naz = 0; naz < Nphi; ++naz) {

	    const double sigma_sb = constants::sigma;
	    const double c = constants::c;
		const double mu = pvte::get_mu(data, nr, naz);
	    const double gamma =
		pvte::get_gamma_eff(data, nr, naz);

	    const double Rgas = constants::R;

		const double H = data[t_data::SCALE_HEIGHT](nr, naz);

		const double sigma = data[t_data::SIGMA](nr, naz);
		const double energy = data[t_data::ENERGY](nr, naz);

	    const double inv_pow4 =
		std::pow(mu * (gamma - 1.0) / (Rgas * sigma), 4);
	    double alpha = 1.0 + 2.0 * H * 4.0 * sigma_sb / c * inv_pow4 *
				     std::pow(energy, 3);

		data[t_data::QPLUS](nr, naz) /= alpha;
		data[t_data::QMINUS](nr, naz) /= alpha;
		const double Qplus = data[t_data::QPLUS](nr, naz);
		const double Qminus = data[t_data::QMINUS](nr, naz);

	    double energy_new = energy + dt * (Qplus - Qminus);

	    const double SigmaFloor =
		10.0 * parameters::sigma0 * parameters::sigma_floor;
	    // If the cell is too close to the density floor
	    // we set energy to equilibrium energy
	    if ((sigma < SigmaFloor)) {
		const double tau_eff =
			data[t_data::TAU_EFF](nr, naz);
		const double e4 = Qplus * tau_eff / (2.0 * sigma_sb);
		const double constant = (Rgas / mu * sigma / (gamma - 1.0));
		// energy, where current heating cooling rate are in equilibirum
		const double eq_energy = std::pow(e4, 1.0 / 4.0) * constant;

		data[t_data::QMINUS](nr, naz) = Qplus;
		energy_new = eq_energy;
	    }

		data[t_data::ENERGY](nr, naz) = energy_new;
	}
    }

	SetTemperatureFloorCeilValues(data, __FILE__, __LINE__);
}


static void compute_sound_speed_normal(t_data &data)
{

	const unsigned int Nr = data[t_data::SOUNDSPEED].get_size_radial();
	const unsigned int Nphi = data[t_data::SOUNDSPEED].get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
	    if (parameters::Adiabatic) {
		const double gamma_eff =
			pvte::get_gamma_eff(data, nr, naz);
		const double gamma1 =
			pvte::get_gamma1(data, nr, naz);

		data[t_data::SOUNDSPEED](nr, naz) =
		    std::sqrt(gamma1 * (gamma_eff - 1.0) *
				  data[t_data::ENERGY](nr, naz) /
				  data[t_data::SIGMA](nr, naz));

	    } else if (parameters::Polytropic) {
		const double gamma_eff =
			pvte::get_gamma_eff(data, nr, naz);
		data[t_data::SOUNDSPEED](nr, naz) =
		    std::sqrt(gamma_eff * constants::R / parameters::MU *
				  data[t_data::TEMPERATURE](nr, naz));
	    } else { // isothermal
		const double h0 = parameters::aspectratio_ref;
		const double beta = parameters::flaring_index;
		const double G = constants::G;
		const double vK = std::sqrt(G * hydro_center_mass / Rb[nr]);
		const double h = h0 * std::pow(Rb[nr], beta);
		// This follows from: cs/v_Kepler = H/r
		const double cs = h*vK;
		data[t_data::SOUNDSPEED](nr, naz) = cs;
	    }
	}
    }
}

static void compute_iso_sound_speed_center_of_mass(t_data &data)
{

    const Pair r_cm = data.get_planetary_system().get_center_of_mass();
    const double m_cm = data.get_planetary_system().get_mass();

	const unsigned int Nr = data[t_data::SOUNDSPEED].get_size_radial();
	const unsigned int Nphi = data[t_data::SOUNDSPEED].get_size_azimuthal();

    // Cs^2 = h^2 * vk * r ^ 2*Flaring
	#pragma omp parallel for collapse(2)
	for (unsigned int n_rad = 0; n_rad < Nr; ++n_rad) {
	for (unsigned int n_az = 0; n_az < Nphi; ++n_az) {

	    const int cell = get_cell_id(n_rad, n_az);
	    const double x = CellCenterX->Field[cell];
	    const double y = CellCenterY->Field[cell];

	    /// since the mass is distributed homogeniously distributed
	    /// inside the cell, we assume that the planet is always at
	    /// least cell_size / 2 plus planet radius away from the gas
	    /// this is an rough estimate without explanation
	    /// alternatively you can think about it yourself
	    const double min_dist =
		0.5 * std::max(Rsup[n_rad] - Rinf[n_rad], Rmed[n_rad] * dphi);

	    const double dx = x - r_cm.x;
	    const double dy = y - r_cm.y;

	    const double dist = std::max(
		std::sqrt(std::pow(dx, 2) + std::pow(dy, 2)), min_dist);

	    const double Cs2 = constants::G * m_cm / dist;

	    const double Cs =
		parameters::aspectratio_ref * std::pow(dist, parameters::flaring_index) * std::sqrt(Cs2);
	    data[t_data::SOUNDSPEED](n_rad, n_az) = Cs;
	}
    }
}

static void compute_iso_sound_speed_nbody(t_data &data, const double current_time)
{

    static const unsigned int N_planets =
	data.get_planetary_system().get_number_of_planets();

    // setup planet data
    for (unsigned int k = 0; k < N_planets; k++) {
	t_planet &planet = data.get_planetary_system().get_planet(k);
	g_mpl[k] = planet.get_rampup_mass(current_time);
	g_xpl[k] = planet.get_x();
	g_ypl[k] = planet.get_y();
	g_rpl[k] = planet.get_planet_radial_extend();
    }

    assert(N_planets > 1);

	const unsigned int Nr = data[t_data::SOUNDSPEED].get_size_radial();
	const unsigned int Nphi = data[t_data::SOUNDSPEED].get_size_azimuthal();

    // Cs^2 = h^2 * vk * r ^ 2*Flaring
	#pragma omp parallel for collapse(2)
	for (unsigned int n_rad = 0; n_rad < Nr; ++n_rad) {
	for (unsigned int n_az = 0; n_az < Nphi; ++n_az) {

	    const int cell = get_cell_id(n_rad, n_az);
	    const double x = CellCenterX->Field[cell];
	    const double y = CellCenterY->Field[cell];

		double Cs2 = 0.0;

		for (unsigned int k = 0; k < N_planets; k++) {

			/// since the mass is distributed homogeniously distributed
			/// inside the cell, we assume that the planet is always at
			/// least cell_size / 2 plus planet radius away from the gas
			/// this is an rough estimate without explanation
			/// alternatively you can think about it yourself
			const double min_dist =
				0.5 * std::max(Rsup[n_rad] - Rinf[n_rad],
					   Rmed[n_rad] * dphi) +
				g_rpl[k];

			const double dx = x - g_xpl[k];
			const double dy = y - g_ypl[k];

			const double dist = std::max(
				std::sqrt(std::pow(dx, 2) + std::pow(dy, 2)), min_dist);

		Cs2 += (parameters::aspectratio_ref * parameters::aspectratio_ref * 
				std::pow(dist, 2.0*parameters::flaring_index)
				* constants::G * g_mpl[k]) / dist;
	    }

		const double Cs = std::sqrt(Cs2);
	    data[t_data::SOUNDSPEED](n_rad, n_az) = Cs;
	}
    }
}

void compute_sound_speed(t_data &data, const double current_time)
{
    if (parameters::Adiabatic || parameters::Polytropic) {
	compute_sound_speed_normal(data);
    }

    if (parameters::Locally_Isothermal) {
	switch (parameters::aspectratio_mode) {
	case 0:
		compute_sound_speed_normal(data);
	    break;
	case 1:
		compute_iso_sound_speed_nbody(data, current_time);
	    break;
	case 2:
		compute_iso_sound_speed_center_of_mass(data);
	    break;
	default:
		compute_sound_speed_normal(data);
	}
    }
}

void compute_scale_height_old(t_data &data)
{

	#pragma omp parallel for
    for (unsigned int n_radial = 0;
	 n_radial <= data[t_data::SCALE_HEIGHT].get_max_radial(); ++n_radial) {
	const double inv_omega_kepler = 1.0 / calculate_omega_kepler(Rb[n_radial]);

	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::SCALE_HEIGHT].get_max_azimuthal();
	     ++n_azimuthal) {

	    if (parameters::Adiabatic || parameters::Polytropic) {
		// h = H/r = c_s,iso / v_k = c_s/sqrt(gamma) / v_k
		// H = h*r = c_s,iso / W_k = c_s/sqrt(gamma) / W_k
		const double gamma1 =
		    pvte::get_gamma1(data, n_radial, n_azimuthal);
		data[t_data::SCALE_HEIGHT](n_radial, n_azimuthal) =
		    data[t_data::SOUNDSPEED](n_radial, n_azimuthal) /
		    (std::sqrt(gamma1)) * inv_omega_kepler;
	    } else {
		// h = H/r = c_s/v_k
		// H = h*r = c_s/v_k * R = c_s / W_k
		data[t_data::SCALE_HEIGHT](n_radial, n_azimuthal) =
		    data[t_data::SOUNDSPEED](n_radial, n_azimuthal) *
		    inv_omega_kepler;
	    }
		if(parameters::heating_star_enabled || parameters::self_gravity || parameters::disk_feedback || parameters::body_force_from_potential){
			const double h = data[t_data::SCALE_HEIGHT](n_radial, n_azimuthal) / Rb[n_radial];
		data[t_data::ASPECTRATIO](n_radial, n_azimuthal) = h;
		}
	}
    }
}

void compute_scale_height_nbody(t_data &data, const double current_time)
{

    static const unsigned int N_planets =
	data.get_planetary_system().get_number_of_planets();
    // setup planet data
    for (unsigned int k = 0; k < N_planets; k++) {
	const t_planet &planet = data.get_planetary_system().get_planet(k);
	g_mpl[k] = planet.get_rampup_mass(current_time);
	g_xpl[k] = planet.get_x();
	g_ypl[k] = planet.get_y();
	g_rpl[k] = planet.get_planet_radial_extend();
    }

	const unsigned int Nr = data[t_data::SCALE_HEIGHT].get_size_radial();
	const unsigned int Nphi = data[t_data::SCALE_HEIGHT].get_size_azimuthal();

    // h = H/r
    // H = = c_s,iso / (GM/r^3) = c_s/sqrt(gamma) / / (GM/r^3)
    // for an Nbody system, H^-2 = sum_n (H_n)^-2
    // See Günter & Kley 2003 Eq. 8, but beware of wrong extra square.
    // Better see Thun et al. 2017 Eq. 8 instead.
	#pragma omp parallel for collapse(2)
	for (unsigned int n_rad = 0; n_rad < Nr; ++n_rad) {
	for (unsigned int n_az = 0; n_az < Nphi; ++n_az) {

	    const int cell = get_cell_id(n_rad, n_az);
	    const double x = CellCenterX->Field[cell];
	    const double y = CellCenterY->Field[cell];
	    const double cs2 =
		std::pow(data[t_data::SOUNDSPEED](n_rad, n_az), 2);

		double inv_H2 = 0.0; // inverse scale height squared
		double inv_h2 = 0.0; // inverse aspectratio squared

	    for (unsigned int k = 0; k < N_planets; k++) {

		/// since the mass is distributed homogeniously distributed
		/// inside the cell, we assume that the planet is always at
		/// least cell_size / 2 plus planet radius away from the gas
		/// this is an rough estimate without explanation
		/// alternatively you can think about it yourself
		const double min_dist =
		    0.5 * std::max(Rsup[n_rad] - Rinf[n_rad],
				   Rmed[n_rad] * dphi) +
			g_rpl[k];

		const double dx = x - g_xpl[k];
		const double dy = y - g_ypl[k];

		const double dist = std::max(
			std::sqrt(std::pow(dx, 2) + std::pow(dy, 2)), min_dist);
		const double dist3 = std::pow(dist, 3);

		// H^2 = (GM / dist^3 / Cs_iso^2)^-1
		if (parameters::Adiabatic || parameters::Polytropic) {
		    const double gamma1 = pvte::get_gamma1(data, n_rad, n_az);
			const double tmp_inv_H2 =
			constants::G * g_mpl[k] * gamma1 / (dist3 * cs2);
		    inv_H2 += tmp_inv_H2;

			if(parameters::heating_star_enabled || parameters::self_gravity || parameters::disk_feedback || parameters::body_force_from_potential){
			const double tmp_inv_h2 =
			constants::G * g_mpl[k] * gamma1 / (dist * cs2);
			inv_h2 += tmp_inv_h2;
			}

		} else {
		    const double tmp_inv_H2 =
			constants::G * g_mpl[k] / (dist3 * cs2);
		    inv_H2 += tmp_inv_H2;

			if(parameters::heating_star_enabled || parameters::self_gravity || parameters::disk_feedback || parameters::body_force_from_potential){
			const double tmp_inv_h2 =
			constants::G * g_mpl[k] / (dist * cs2);
			inv_h2 += tmp_inv_h2;
			}
		}
	    }

	    const double H = std::sqrt(1.0 / inv_H2);
	    data[t_data::SCALE_HEIGHT](n_rad, n_az) = H;

		if(parameters::heating_star_enabled || parameters::self_gravity || parameters::disk_feedback || parameters::body_force_from_potential){
			const double h = std::sqrt(1.0 / inv_h2);
			data[t_data::ASPECTRATIO](n_rad, n_az) = h;
		}
	}
    }
}

void compute_scale_height_center_of_mass(t_data &data)
{

    const Pair r_cm = data.get_planetary_system().get_center_of_mass();
    const double m_cm = data.get_planetary_system().get_mass();

	const unsigned int Nr = data[t_data::SCALE_HEIGHT].get_size_radial();
	const unsigned int Nphi = data[t_data::SCALE_HEIGHT].get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int n_rad = 0; n_rad < Nr; ++n_rad) {
	for (unsigned int n_az = 0; n_az < Nphi; ++n_az) {

	    const int cell = get_cell_id(n_rad, n_az);
	    const double x = CellCenterX->Field[cell];
	    const double y = CellCenterY->Field[cell];
		const double cs = data[t_data::SOUNDSPEED](n_rad, n_az);

	    // const double min_dist =
	    //	0.5 * std::max(Rsup[n_rad] - Rinf[n_rad],
	    //		   Rmed[n_rad] * dphi);

	    const double dx = x - r_cm.x;
	    const double dy = y - r_cm.y;

	    // const double dist = std::max(
	    //	std::sqrt(std::pow(dx, 2) + std::pow(dy, 2)), min_dist);
	    const double dist = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));

		// h^2 = Cs_iso / vk = (Cs_iso^2 / (GM / dist))
		// H^2 = Cs_iso / Omegak = (Cs_iso^2 / (GM / dist^3))
		// H = h * dist
	    if (parameters::Adiabatic || parameters::Polytropic) {
		/// Convert sound speed to isothermal sound speed cs,iso = cs /
		/// sqrt(gamma)
		const double gamma1 = pvte::get_gamma1(data, n_rad, n_az);
		const double h = cs * std::sqrt(dist / (constants::G * m_cm * gamma1));

		if(parameters::heating_star_enabled || parameters::self_gravity || parameters::disk_feedback || parameters::body_force_from_potential){
		data[t_data::ASPECTRATIO](n_rad, n_az) = h;
		}
		const double H = dist * h;
		data[t_data::SCALE_HEIGHT](n_rad, n_az) = H;

	    } else { // locally isothermal
		const double h = cs * std::sqrt(dist / (constants::G * m_cm));
		if(parameters::heating_star_enabled || parameters::self_gravity || parameters::disk_feedback || parameters::body_force_from_potential){
		data[t_data::ASPECTRATIO](n_rad, n_az) = h;
		}
		const double H = dist * h;
		data[t_data::SCALE_HEIGHT](n_rad, n_az) = H;
	    }
	}
    }
}

static void adjust_scale_height_for_sg(t_data &data) {

	t_polargrid &H = data[t_data::SCALE_HEIGHT];

	const unsigned int Nr = H.get_size_radial();
	const unsigned int Naz = H.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Naz; ++naz) {
		const double Q = data[t_data::TOOMRE](nr, naz);
		// Hsg = sqrt(2/pi) * H  * f(Q)
		// f(Q) = pi / (4Q) * (sqrt(1 + 8*Q^2/pi) - 1)
		const double f = M_PI * (std::sqrt(1 + 8*Q*Q/M_PI) - 1) / (4*Q);
		H(nr, naz) *= f*std::sqrt(2/M_PI);
	}
	}
}

void compute_scale_height(t_data &data, const double current_time)
{
    switch (parameters::aspectratio_mode) {
    case 0:
	compute_scale_height_old(data);
	break;
    case 1:
	compute_scale_height_nbody(data, current_time);
	break;
    case 2:
	compute_scale_height_center_of_mass(data);
	break;
    default:
	compute_scale_height_old(data);
    }
	
	if (parameters::self_gravity && 
		parameters::self_gravity_mode == parameters::t_sg::sg_BK) {
        compute::toomreQ(data);
        adjust_scale_height_for_sg(data);
    }
}

void compute_pressure(t_data &data)
{

	const unsigned int Nr = data[t_data::PRESSURE].get_size_radial();
	const unsigned int Nphi = data[t_data::PRESSURE].get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {

	    if (parameters::Adiabatic) {
		const double gamma_eff =
			pvte::get_gamma_eff(data, nr, naz);
		data[t_data::PRESSURE](nr, naz) =
		    (gamma_eff - 1.0) *
			data[t_data::ENERGY](nr, naz);
	    } else if (parameters::Polytropic) {
		data[t_data::PRESSURE](nr, naz) =
			data[t_data::SIGMA](nr, naz) *
			std::pow(data[t_data::SOUNDSPEED](nr, naz),
			     2) / parameters::ADIABATICINDEX;
	    } else { // Isothermal
		// since SoundSpeed is not update from initialization, cs
		// remains axisymmetric
		data[t_data::PRESSURE](nr, naz) =
			data[t_data::SIGMA](nr, naz) *
			std::pow(data[t_data::SOUNDSPEED](nr, naz),
			     2);
	    }
	}
    }
}

void compute_temperature(t_data &data)
{
	auto &T = data[t_data::TEMPERATURE];
	auto &Sig = data[t_data::SIGMA];
	auto &E = data[t_data::ENERGY];
	auto &P = data[t_data::PRESSURE];

	const unsigned int Nr = T.get_size_radial();
	const unsigned int Nphi = T.get_size_azimuthal();

	const double Rgas = constants::R;
	const double polyconst = parameters::POLYTROPIC_CONSTANT;

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
	    if (parameters::Adiabatic) {
			const double mu = pvte::get_mu(data, nr, naz);
			const double gamma_eff = pvte::get_gamma_eff(data, nr, naz);
			const double c_v_inv = mu / Rgas * (gamma_eff - 1.0);
			T(nr, naz) = c_v_inv * E(nr, naz) / Sig(nr, naz);
	    } else if (parameters::Polytropic) {
			const double mu = pvte::get_mu(data, nr, naz);
			const double gamma_eff = pvte::get_gamma_eff(data, nr, naz);
			T(nr, naz) = mu / Rgas * polyconst * std::pow(Sig(nr, naz),gamma_eff - 1.0);
	    } else { // Isothermal
			T(nr, naz) = parameters::MU / Rgas * P(nr, naz) / Sig(nr, naz);
	    }
	}
    }
}

void compute_heating_cooling_for_CFL(t_data &data, const double current_time)
{
    if (parameters::Adiabatic) {

	viscosity::update_viscosity(data);
	viscosity::compute_viscous_stress_tensor(data);

	calculate_qminus(data, current_time); // first to calculate teff
	calculate_qplus(data);

	const unsigned int Nr = data[t_data::ENERGY].get_max_radial();
	const unsigned int Nphi = data[t_data::ENERGY].get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nr; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {

		const double sigma_sb = constants::sigma;
		const double c = constants::c;
		const double mu = pvte::get_mu(data, nr, naz);
		const double gamma =
			pvte::get_gamma_eff(data, nr, naz);
		const double Rgas = constants::R;

		const double H =
			data[t_data::SCALE_HEIGHT](nr, naz);

		const double sigma = data[t_data::SIGMA](nr, naz);
		const double energy = data[t_data::ENERGY](nr, naz);

		const double inv_pow4 =
		    std::pow(mu * (gamma - 1.0) / (Rgas * sigma), 4);
		double alpha = 1.0 + 2.0 * H * 4.0 * sigma_sb / c * inv_pow4 *
					 std::pow(energy, 3);

		data[t_data::QPLUS](nr, naz) /= alpha;
		data[t_data::QMINUS](nr, naz) /= alpha;
	    }
	}
    }
}
