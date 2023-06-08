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
#include "Pframeforce.h"
#include "SideEuler.h"
#include "SourceEuler.h"
#include "Theo.h"
#include "TransportEuler.h"
#include "accretion.h"
#include "boundary_conditions.h"
#include "commbound.h"
#include "constants.h"
#include "gas_torques.h"
#include "global.h"
#include "logging.h"
#include "nongnu.h"
#include "opacity.h"
#include "output.h"
#include "parameters.h"
#include "particles/particles.h"
#include "pvte_law.h"
#include "quantities.h"
#include "selfgravity.h"
#include "stress.h"
#include "units.h"
#include "util.h"
#include "viscosity/viscosity.h"
#include "viscosity/artificial_viscosity.h"
#include "simulation.h"
#include "frame_of_reference.h"
#include "cfl.h"
#include "simulation.h"

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
			units::temperature.get_cgs_factor() /
			density(nr, naz) * mu / constants::R *
			(gamma_eff - 1.0),
		    Tmin * units::temperature.get_cgs_factor(), Tmin);
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
			units::temperature.get_cgs_factor() /
			density(nr, naz) * mu / constants::R *
			(gamma_eff - 1.0),
		    Tmax * units::temperature.get_cgs_factor(), Tmax);
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
	if (parameters::ASPECTRATIO_MODE > 0) {
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
	if (parameters::ASPECTRATIO_MODE > 0) {
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
	g_l1pl.resize(N_planets);

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

/**

*/
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

    double supp_torque = 0.0; // for imposed disk drift

	const unsigned int Nphi = data[t_data::ENERGY].get_size_azimuthal();

    // update v_radial with source terms
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = One_no_ghost_vr; nr < MaxMo_no_ghost_vr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
	    // 1/Sigma * dP/dr : Sigma is calculated as a mean value between the
	    // neightbour cells
	    const double gradp =
		2.0 /
		(data[t_data::SIGMA](nr, naz) +
		 data[t_data::SIGMA](nr - 1, naz)) *
		(data[t_data::PRESSURE](nr, naz) -
		 data[t_data::PRESSURE](nr - 1, naz)) *
		InvDiffRmed[nr];

	    // dPhi/dr
	    double gradphi;
	    if (parameters::body_force_from_potential) {
		gradphi = (data[t_data::POTENTIAL](nr, naz) -
			   data[t_data::POTENTIAL](nr - 1, naz)) *
			  InvDiffRmed[nr];
	    } else {
		gradphi = -data[t_data::ACCEL_RADIAL](nr, naz);
	    }

		const unsigned int naz_next = (naz == Nphi-1 ? 0 : naz + 1);

		// v_phi^2/r : v_phi^2 is calculated by a mean in both directions
		double vt2 =
		data[t_data::V_AZIMUTHAL](nr, naz) +
		data[t_data::V_AZIMUTHAL](nr, naz_next) +
		data[t_data::V_AZIMUTHAL](nr - 1, naz) +
		data[t_data::V_AZIMUTHAL](nr - 1, naz_next);
		vt2 = 0.25 * vt2 + Rinf[nr] * refframe::OmegaFrame;
		vt2 = vt2 * vt2;

		// add all terms to new v_radial: v_radial_new = v_radial +
		// dt*(source terms)
		data[t_data::V_RADIAL](nr, naz) =
		data[t_data::V_RADIAL](nr, naz) +
		dt * (-gradp - gradphi + vt2 * InvRinf[nr]);

	}
    }

    // update v_azimuthal with source terms
	#pragma omp parallel for
	for (unsigned int nr = Zero_no_ghost; nr < Max_no_ghost; ++nr) {

	if (parameters::IMPOSEDDISKDRIFT != 0.0) {
		supp_torque = parameters::IMPOSEDDISKDRIFT * 0.5 *
			  std::pow(Rmed[nr], -2.5 + parameters::SIGMASLOPE);
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
		gradphi = -data[t_data::ACCEL_AZIMUTHAL](nr, naz);
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

    if (parameters::Adiabatic) {
    #pragma omp parallel for collapse(2)
    for (unsigned int nr = Zero_no_ghost; nr < Max_no_ghost; ++nr) {
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

        /*
        // Like D'Angelo et al. 2003 eq. 25
        const double P = (gamma - 1.0) * data[t_data::ENERGY](n_radial,
        n_azimuthal); const double dE = dt * (-P*DIV_V + 0.5*(gamma
        - 1.0) * P * dt * std::pow(DIV_V, 2)); const double energy_old =
        data[t_data::ENERGY](n_radial, n_azimuthal); const double
        energy_new = energy_old + dE; data[t_data::ENERGY](n_radial,
        n_azimuthal) = energy_new;
        */

        // Like D'Angelo et al. 2003 eq. 24
        const double energy_old =
            data[t_data::ENERGY](nr, naz);
        const double energy_new =
            energy_old * std::exp(-(gamma - 1.0) * dt * DIV_V);
        data[t_data::ENERGY](nr, naz) = energy_new;

        /*
        // Zeus2D like, see Stone & Norman 1992
        // produces poor results with shock tube test
        const double P = (gamma - 1.0);
        const double energy_old = data[t_data::ENERGY](n_radial,
        n_azimuthal); const double energy_new = energy_old*(1.0 -
        0.5*dt*P*DIV_V)/(1.0 + 0.5*dt*P*DIV_V);
        data[t_data::ENERGY](n_radial, n_azimuthal) = energy_new;
        */
        }
    }
    }

	if(ECC_GROWTH_MONITOR){
		quantities::calculate_disk_delta_ecc_peri(data, delta_ecc_source, delta_peri_source);
	}

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
	
	if (sim::PhysicalTime < rampup_time) {
	    ramping =
		1.0 -
		std::pow(std::cos(sim::PhysicalTime * M_PI / 2.0 / rampup_time), 2);
	}

	const double x = planet.get_x();
	const double y = planet.get_y();
	const double radius = planet.get_planet_radial_extend();
	const double temperature = planet.get_temperature();

	const double l1 = planet.get_dimensionless_roche_radius() *
			planet.get_distance_to_primary();
	double min_dist;
	if(x*x + y*y > 1e-10){
		min_dist = l1 * parameters::klahr_smoothing_radius;
	} else {
		min_dist = 0.0;
	}

	const unsigned int Nrad = data[t_data::QPLUS].get_max_radial();
	const unsigned int Naz = data[t_data::QPLUS].get_max_azimuthal();
	// Simple star heating (see Masterthesis Alexandros Ziampras)

	#pragma omp parallel for collapse(2)
	for (unsigned int nrad = 1; nrad < Nrad; ++nrad) {
	for (unsigned int naz = 0; naz <= Naz; ++naz) {
		const unsigned int ncell = nrad * data[t_data::SIGMA].get_size_azimuthal() + naz;
		const double xc = CellCenterX->Field[ncell];
		const double yc = CellCenterY->Field[ncell];
		const double distance_measured = std::sqrt(std::pow(x - xc, 2) + std::pow(y - yc, 2));
		const double distance = std::max(distance_measured, min_dist);

		const double HoverR = data[t_data::ASPECTRATIO](nrad, naz);
		const double sigma = constants::sigma.get_code_value();
		const double tau_eff = data[t_data::TAU_EFF](nrad, naz);
		const double eps = 0.5; // TODO: add a parameter
		// choose according to Chiang & Goldreich (1997)
		const double dlogH_dlogr = 9.0 / 7.0;
		
		// irradiation contribution near and far from the star
		// see D'Angelo & Marzari 2012 (doi:10.1088/0004-637X/757/1/5), 
		const double roverd = distance < radius ? 1.0 : radius/distance;
		const double W_G = 0.4 * roverd + HoverR * (dlogH_dlogr - 1.0);

		// use eq. 7 from Menou & Goodman (2004) (rearranged), Qirr
		// = 2*(1-eps)*L_star/(4 pi r^2)*(dlogH/dlogr - 1) * H/r *
		// 1/Tau_eff here we use (1-eps) =
		// parameters::heating_star_factor L_star = 4 pi R_star^2
		// sigma_sb T_star^4
		// with modifications from D'Angelo & Marzari 2012
		// for near and far field
		double qplus = 2.0 * (1.0 - eps);
		qplus *= sigma * std::pow(temperature, 4) * std::pow(roverd, 2); // *L_star/(4 pi r^2)
		qplus *= W_G;
		qplus /= tau_eff;			// * 1/Tau_eff
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
		if (!parameters::cooling_radiative_enabled) {
			// TODO: make it properly!
			die("Need to calulate Tau_eff first!\n"); 
		}
		irradiation(data);
	}    
}

void calculate_qminus(t_data &data, const double current_time)
{
    // clear up all Qminus terms
    data[t_data::QMINUS].clear();

	const unsigned int Nr = data[t_data::QMINUS].get_max_radial(); // = Size - 1
	const unsigned int Nphi = data[t_data::QMINUS].get_size_azimuthal();

    // beta cooling
    if (parameters::cooling_beta_enabled) {

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
		if (parameters::cooling_beta_initial) {
		    const double sigma =
			data[t_data::SIGMA](nr, naz);
		    const double sigma0 =
			data[t_data::SIGMA0](nr, naz);
		    const double E0 =
			data[t_data::ENERGY0](nr, naz);
		    delta_E -= E0 / sigma0 * sigma;
		}
		if (parameters::cooling_beta_aspect_ratio) {
		    const double sigma =
			data[t_data::SIGMA](nr, naz);
		    const double E0 =
			1.0 / (parameters::ADIABATICINDEX - 1.0) *
			std::pow(parameters::ASPECTRATIO_REF, 2) *
			std::pow(Rmed[nr], 2.0 * parameters::FLARINGINDEX - 1.0) *
			constants::G * hydro_center_mass * sigma;
		    delta_E -= E0;
		}
		const double qminus = delta_E * omega_k * beta_inv;


		data[t_data::QMINUS](nr, naz) += qminus;
	    }
	}
    }

    // local radiative cooling
    if (parameters::cooling_radiative_enabled) {
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nr; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {
		// calculate Rosseland mean opacity kappa. opaclin needs values
		// in cgs units
		const double temperatureCGS =
			data[t_data::TEMPERATURE](nr, naz) *
		    units::temperature;

		const double H =
			data[t_data::SCALE_HEIGHT](nr, naz);

		const double densityCGS =
			data[t_data::SIGMA](nr, naz) /
		    (parameters::density_factor * H) * units::density;

		const double kappaCGS =
		    opacity::opacity(densityCGS, temperatureCGS);

		data[t_data::KAPPA](nr, naz) =
		    parameters::kappa_factor * kappaCGS *
		    units::opacity.get_inverse_cgs_factor();

		// mean vertical optical depth: tau = 1/2 kappa Sigma
		data[t_data::TAU](nr, naz) =
		    parameters::tau_factor *
		    (1.0 / parameters::density_factor) *
			data[t_data::KAPPA](nr, naz) *
			data[t_data::SIGMA](nr, naz);


		if(parameters::heating_star_enabled){
		//  irradiated disk tau_eff = 3/8 tau + 1/2 + 1/(4*tau+tau_min)
		//  compare D'Angelo & Marzari 2012
		data[t_data::TAU_EFF](nr, naz) =
			3.0 / 8.0 * data[t_data::TAU](nr, naz) +
			0.5 +
		    1.0 /
			(4.0 * data[t_data::TAU](nr, naz) + parameters::tau_min);
		} else {
			//  non irradiated disk tau_eff = 3/8 tau + sqrt(3)/4 + 1/(4*tau+tau_min)
			data[t_data::TAU_EFF](nr, naz) =
				3.0 / 8.0 * data[t_data::TAU](nr, naz) +
				std::sqrt(3.0) / 4.0 +
				1.0 /
				(4.0 * data[t_data::TAU](nr, naz) + parameters::tau_min);
		}

		if (parameters::opacity ==
		    parameters::opacity_simple) { // Compare D'Angelo et. al
						  // 2003 eq.(28)
			data[t_data::TAU_EFF](nr, naz) =
			3.0 / 8.0 * data[t_data::TAU](nr, naz);
		}
		// Q = factor 2 sigma_sb T^4 / tau_eff

		const double factor = parameters::cooling_radiative_factor;
		const double sigma_sb = constants::sigma.get_code_value();
		const double T4 = std::pow(
			data[t_data::TEMPERATURE](nr, naz), 4);
		const double tau_eff =
			data[t_data::TAU_EFF](nr, naz);
		const double Tmin4 =
			std::pow(parameters::minimum_temperature, 4);

		const double qminus =
		    factor * 2 * sigma_sb * (T4 - Tmin4) / tau_eff;

		data[t_data::QMINUS](nr, naz) += qminus;
	    }
	}
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
	parameters::radiative_diffusion_enabled) {

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

static inline double flux_limiter(const double R)
{
    // flux limiter
    if (R <= 2) {
	return 2.0 / (3 + std::sqrt(9 + 10 * std::pow(R, 2)));
    } else {
	return 10.0 / (10 * R + 9 + std::sqrt(180 * R + 81));
    }
}

void radiative_diffusion(t_data &data, const double current_time, const double dt)
{
    static bool grids_allocated = false;
    static t_polargrid Ka, Kb;
    static t_polargrid A, B, C, D, E;
    static t_polargrid Told;
    static double *SendInnerBoundary, *SendOuterBoundary, *RecvInnerBoundary,
	*RecvOuterBoundary;

    if (!grids_allocated) {
	Ka.set_vector(true);
	Ka.set_size(data.get_n_radial(), data.get_n_azimuthal());
	Kb.set_scalar(true);
	Kb.set_size(data.get_n_radial(), data.get_n_azimuthal());

	A.set_scalar(true);
	A.set_size(data.get_n_radial(), data.get_n_azimuthal());
	B.set_scalar(true);
	B.set_size(data.get_n_radial(), data.get_n_azimuthal());
	C.set_scalar(true);
	C.set_size(data.get_n_radial(), data.get_n_azimuthal());
	D.set_scalar(true);
	D.set_size(data.get_n_radial(), data.get_n_azimuthal());
	E.set_scalar(true);
	E.set_size(data.get_n_radial(), data.get_n_azimuthal());

	Told.set_scalar(true);
	Told.set_size(data.get_n_radial(), data.get_n_azimuthal());

	// create arrays for communcation
	SendInnerBoundary =
	    (double *)malloc(NAzimuthal * CPUOVERLAP * sizeof(double));
	SendOuterBoundary =
	    (double *)malloc(NAzimuthal * CPUOVERLAP * sizeof(double));
	RecvInnerBoundary =
	    (double *)malloc(NAzimuthal * CPUOVERLAP * sizeof(double));
	RecvOuterBoundary =
	    (double *)malloc(NAzimuthal * CPUOVERLAP * sizeof(double));

	grids_allocated = true;
    }

    auto &Temperature = data[t_data::TEMPERATURE];
    auto &Sigma = data[t_data::SIGMA];
    auto &Energy = data[t_data::ENERGY];
    auto &Scale_height = data[t_data::SCALE_HEIGHT];

	const unsigned int Nphi = Kb.get_size_azimuthal();

	// We set minimum Temperature here such that H (thru Cs) is also computed with the minimum Temperature
	#pragma omp parallel for
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const unsigned int nr_max = Energy.get_max_radial();

	if (CPU_Rank == 0 &&
	    parameters::boundary_inner == parameters::boundary_condition_open) {
	    const double Tmin = parameters::minimum_temperature;

	    const double mu = pvte::get_mu(data, 1, naz);
	    const double gamma_eff = pvte::get_gamma_eff(data, 1, naz);
	    Sigma(0, naz) = Sigma(1, naz);

	    const double minimum_energy =
		Tmin * Sigma(1, naz) / mu * constants::R / (gamma_eff - 1.0);

	    Energy(0, naz) = minimum_energy;
	}

	if (CPU_Rank == CPU_Highest &&
	    parameters::boundary_outer == parameters::boundary_condition_open) {
	    const double Tmin = parameters::minimum_temperature; 

	    const double mu = pvte::get_mu(data, nr_max - 1, naz);
	    const double gamma_eff = pvte::get_gamma_eff(data, nr_max - 1, naz);
	    Sigma(nr_max, naz) = Sigma(nr_max - 1, naz);

	    const double minimum_energy = Tmin * Sigma(nr_max - 1, naz) / mu *
					  constants::R / (gamma_eff - 1.0);

	    Energy(nr_max, naz) = minimum_energy;
	}
    }

    // update temperature, soundspeed and aspect ratio
	compute_temperature(data);
	compute_sound_speed(data, current_time);
	compute_scale_height(data, current_time);

	const unsigned int NrKa = Ka.get_size_radial();
    // calcuate Ka for K(i/2,j)
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < NrKa - 1; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
	    const unsigned int n_azimuthal_plus =
		(naz == Ka.get_max_azimuthal() ? 0 : naz + 1);
	    const unsigned int n_azimuthal_minus =
		(naz == 0 ? Ka.get_max_azimuthal() : naz - 1);

	    // average temperature radially
	    const double temperature =
		0.5 * (Temperature(nr - 1, naz) + Temperature(nr, naz));
	    const double density = 0.5 * (Sigma(nr - 1, naz) + Sigma(nr, naz));
	    const double scale_height =
		0.5 * (Scale_height(nr - 1, naz) + Scale_height(nr, naz));

	    const double temperatureCGS = temperature * units::temperature;
	    const double H = scale_height;
	    const double densityCGS =
		density / (parameters::density_factor * H) * units::density;

	    const double kappaCGS =
		opacity::opacity(densityCGS, temperatureCGS);
	    const double kappa = parameters::kappa_factor * kappaCGS *
				 units::opacity.get_inverse_cgs_factor();

	    const double denom = 1.0 / (density * kappa);

	    // Levermore & Pomraning 1981
	    // R = 4 |nabla T\/T * 1/(rho kappa)
	    const double dT_dr =
		(Temperature(nr, naz) - Temperature(nr - 1, naz)) *
		InvDiffRmed[nr];
	    const double dT_dphi =
		InvRinf[nr] *
		(0.5 * (Temperature(nr - 1, n_azimuthal_plus) +
			Temperature(nr, n_azimuthal_plus)) -
		 0.5 * (Temperature(nr - 1, n_azimuthal_minus) +
			Temperature(nr, n_azimuthal_minus))) /
		(2.0 * dphi);

	    const double nabla_T =
		std::sqrt(std::pow(dT_dr, 2) + std::pow(dT_dphi, 2));

	    const double R = 4.0 * nabla_T / temperature * denom * H *
			     parameters::density_factor;

	    const double lambda = flux_limiter(R);

	    Ka(nr, naz) = 8.0 * 4.0 * constants::sigma.get_code_value() *
			  lambda * H * H * std::pow(temperature, 3) * denom;
	}
    }

	#pragma omp parallel for
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const unsigned int nr_max = Ka.get_max_radial();
		if(CPU_Rank == CPU_Highest && parameters::boundary_outer == parameters::boundary_condition_reflecting){
		Ka(nr_max-1, naz) = 0.0;
		}


		if(CPU_Rank == 0 && parameters::boundary_inner == parameters::boundary_condition_reflecting){
		Ka(1, naz) = 0.0;
		}
	}

	#pragma omp parallel for
	// Similar to Tobi's original implementation
	for (unsigned int naz = 0; naz < Ka.get_size_azimuthal(); ++naz) {
		const unsigned int nr_max = Ka.get_max_radial();
		if(CPU_Rank == CPU_Highest && !(parameters::boundary_outer == parameters::boundary_condition_reflecting || parameters::boundary_outer == parameters::boundary_condition_open)){
		Ka(nr_max-1, naz) = Ka(nr_max-2, naz);
		}


		if(CPU_Rank == 0 && !(parameters::boundary_inner == parameters::boundary_condition_reflecting || parameters::boundary_inner == parameters::boundary_condition_open)){
		Ka(1, naz) = Ka(2, naz);
		}
	}


    // Similar to Tobi's original implementation
    for (unsigned int naz = 0; naz < Ka.get_size_azimuthal(); ++naz) {
	const unsigned int nr_max = Ka.get_max_radial();
	if (CPU_Rank == CPU_Highest &&
	    !(parameters::boundary_outer ==
		  parameters::boundary_condition_reflecting ||
	      parameters::boundary_outer ==
		  parameters::boundary_condition_open)) {
	    Ka(nr_max - 1, naz) = Ka(nr_max - 2, naz);
	}

	if (CPU_Rank == 0 && !(parameters::boundary_inner ==
				   parameters::boundary_condition_reflecting ||
			       parameters::boundary_inner ==
				   parameters::boundary_condition_open)) {
	    Ka(1, naz) = Ka(2, naz);
	}
    }

    // calcuate Kb for K(i,j/2)
	const unsigned int NrKb = Kb.get_size_radial();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < NrKb-1; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
	    // unsigned int n_azimuthal_plus = (n_azimuthal ==
	    // Kb.get_max_azimuthal() ? 0 : n_azimuthal + 1);
		const unsigned int naz_prev =
		(naz == 0 ? Nphi-1 : naz - 1);

	    // average temperature azimuthally
	    const double temperature =
		0.5 * (Temperature(nr, naz_prev) + Temperature(nr, naz));
		const double density = 0.5 * (Sigma(nr, naz_prev) + Sigma(nr, naz));
	    const double scale_height =
		0.5 * (Scale_height(nr, naz_prev) + Scale_height(nr, naz));

	    const double temperatureCGS = temperature * units::temperature;
	    const double H = scale_height;
	    const double densityCGS =
		density / (parameters::density_factor * H) * units::density;

	    const double kappaCGS =
		opacity::opacity(densityCGS, temperatureCGS);
	    const double kappa = parameters::kappa_factor * kappaCGS *
				 units::opacity.get_inverse_cgs_factor();

	    const double denom = 1.0 / (density * kappa);

	    // Levermore & Pomraning 1981
	    // R = 4 |nabla T\/T * 1/(rho kappa)
	    const double dT_dr =
		(0.5 * (Temperature(nr - 1, naz_prev) + Temperature(nr - 1, naz)) -
		 0.5 * (Temperature(nr + 1, naz_prev) + Temperature(nr + 1, naz))) /
		(Ra[nr - 1] - Ra[nr + 1]);
	    const double dT_dphi =
		InvRmed[nr] * (Temperature(nr, naz) - Temperature(nr, naz_prev)) /
		dphi;

	    const double nabla_T =
		std::sqrt(std::pow(dT_dr, 2) + std::pow(dT_dphi, 2));

	    const double R = 4.0 * nabla_T / temperature * denom * H *
			     parameters::density_factor;

	    const double lambda = flux_limiter(R);
	    /*if (n_radial == 4) {
		    printf("kb:
	    phi=%lg\tR=%lg\tlambda=%lg\tdphi=%lg\tdr=%lg\tnabla=%lg\tT=%lg\tH=%lg\n",
	    dphi*n_azimuthal, R, lambda,dT_dphi,dT_dr,nabla_T,temperature,H);
	    }*/

	    Kb(nr, naz) = 8.0 * 4.0 * constants::sigma.get_code_value() *
			  lambda * H * H * std::pow(temperature, 3) * denom;
	    // Kb(n_radial, n_azimuthal)
	    // = 16.0*parameters::density_factor*constants::sigma.get_code_value()*lambda*H*pow3(temperature)*denom;
	}
    }

    const double c_v = constants::R / (parameters::MU * (parameters::ADIABATICINDEX - 1.0));

    // calculate A,B,C,D,E
	const unsigned int Nr = Temperature.get_size_radial();
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nr - 1; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const double Sig = Sigma(nr, naz);
	    const double common_factor =
		-dt * parameters::density_factor / (Sig * c_v);

	    // 2/(dR^2)
	    const double common_AC =
		common_factor * 2.0 /
		(std::pow(Ra[nr + 1], 2) - std::pow(Ra[nr], 2));
	    A(nr, naz) = common_AC * Ka(nr, naz) * Ra[nr] * InvDiffRmed[nr];
	    C(nr, naz) =
		common_AC * Ka(nr + 1, naz) * Ra[nr + 1] * InvDiffRmed[nr + 1];

	    // 1/(r^2 dphi^2)
	    const double common_DE =
		common_factor / (std::pow(Rb[nr], 2) * std::pow(dphi, 2));
	    D(nr, naz) = common_DE * Kb(nr, naz);
	    E(nr, naz) =
		common_DE * Kb(nr, naz == Kb.get_max_azimuthal() ? 0 : naz + 1);

	    B(nr, naz) =
		-A(nr, naz) - C(nr, naz) - D(nr, naz) - E(nr, naz) + 1.0;

	    Told(nr, naz) = Temperature(nr, naz);

	    /*double energy_change = dt*data[t_data::QPLUS](n_radial,
	    n_azimuthal)
		- dt*data[t_data::QMINUS](n_radial, n_azimuthal)
		- dt*data[t_data::P_DIVV](n_radial, n_azimuthal);

	    double temperature_change =
		MU/R*(parameters::ADIABATICINDEX-1.0)*energy_change/Sigma(n_radial,n_azimuthal);
	    Told(n_radial, n_azimuthal) += temperature_change;

	    if (Told(n_radial, n_azimuthal) <
	    parameters::minimum_temperature*units::temperature.get_inverse_cgs_factor())
		{ Temperature(n_radial, n_azimuthal) =
	    parameters::minimum_temperature*units::temperature.get_inverse_cgs_factor();
	    }
	    */
	}
    }

    static unsigned int old_iterations =
	parameters::radiative_diffusion_max_iterations;
    static int direction = 1;
    static double omega = parameters::radiative_diffusion_omega;

    unsigned int iterations = 0;
	double absolute_norm = std::numeric_limits<double>::max();
	double norm_change = std::numeric_limits<double>::max();

    const int l = CPUOVERLAP * NAzimuthal;
    const int oo = (Temperature.Nrad - CPUOVERLAP) * NAzimuthal;
    const int o = (Temperature.Nrad - 2 * CPUOVERLAP) * NAzimuthal;

    // do SOR
    while ((norm_change > 1e-12) &&
	   (parameters::radiative_diffusion_max_iterations > iterations)) {
	// if ((CPU_Rank == CPU_Highest) && parameters::boundary_outer ==
	// parameters::boundary_condition_open) {
	// 	// set temperature to T_min in outermost cells
	// 	for (unsigned int n_azimuthal = 0; n_azimuthal <=
	// Temperature.get_max_azimuthal(); ++n_azimuthal) {
	// 		Temperature(Temperature.get_max_radial(),
	// n_azimuthal) =
	// parameters::minimum_temperature*units::temperature.get_inverse_cgs_factor();
	// 	}
	// }

	// if ((CPU_Rank == 0) && parameters::boundary_inner ==
	// parameters::boundary_condition_open) {
	// 	// set temperature to T_min in innermost cells
	// 	for (unsigned int n_azimuthal = 0; n_azimuthal <=
	// Temperature.get_max_azimuthal(); ++n_azimuthal) {
	// 		Temperature(0, n_azimuthal) =
	// parameters::minimum_temperature*units::temperature.get_inverse_cgs_factor();
	// 	}
	// }
	boundary_conditions::apply_boundary_condition(data, current_time, 0.0, false);

	norm_change = absolute_norm;
	absolute_norm = 0.0;
    /// TODO: Cannot be OpenMP parallelized due to Temperature being iteratively computed ??
#pragma omp parallel for collapse(2) reduction(+ : absolute_norm)
	for (unsigned int nr = 1; nr < Nr - 1; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {

		const double old_value = Temperature(nr, naz);
		const unsigned int naz_next =
		    (naz == Temperature.get_max_azimuthal() ? 0 : naz + 1);
		const unsigned int naz_prev =
		    (naz == 0 ? Temperature.get_max_azimuthal() : naz - 1);

		Temperature(nr, naz) =
		    (1.0 - omega) * Temperature(nr, naz) -
		    omega / B(nr, naz) *
			(A(nr, naz) * Temperature(nr - 1, naz) +
			 C(nr, naz) * Temperature(nr + 1, naz) +
			 D(nr, naz) * Temperature(nr, naz_prev) +
			 E(nr, naz) * Temperature(nr, naz_next) - Told(nr, naz));

		// only non ghostcells to norm and don't count overlap cell's
		// twice
		const bool isnot_ghostcell_rank_0 =
		    nr > ((CPU_Rank == 0) ? GHOSTCELLS_B : CPUOVERLAP);
		const bool isnot_ghostcell_rank_highest =
		    (nr <
		     (Temperature.get_max_radial() -
		      ((CPU_Rank == CPU_Highest) ? GHOSTCELLS_B : CPUOVERLAP)));

		if (isnot_ghostcell_rank_0 && isnot_ghostcell_rank_highest) {
		    absolute_norm +=
			std::pow(old_value - Temperature(nr, naz), 2);
		}
	    }
	}

    const double tmp = absolute_norm;
	MPI_Allreduce(&tmp, &absolute_norm, 1, MPI_DOUBLE, MPI_SUM,
		      MPI_COMM_WORLD);
	absolute_norm = std::sqrt(absolute_norm) / (GlobalNRadial * NAzimuthal);

	norm_change = fabs(absolute_norm - norm_change);
	iterations++;

	// communicate with other nodes
	memcpy(SendInnerBoundary, Temperature.Field + l, l * sizeof(double));
	memcpy(SendOuterBoundary, Temperature.Field + o, l * sizeof(double));

	MPI_Request req1, req2, req3, req4;

	if (CPU_Rank % 2 == 0) {
	    if (CPU_Rank != 0) {
		MPI_Isend(SendInnerBoundary, NAzimuthal * CPUOVERLAP,
			  MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
		MPI_Irecv(RecvInnerBoundary, NAzimuthal * CPUOVERLAP,
			  MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req2);
	    }
	    if (CPU_Rank != CPU_Highest) {
		MPI_Isend(SendOuterBoundary, NAzimuthal * CPUOVERLAP,
			  MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req3);
		MPI_Irecv(RecvOuterBoundary, NAzimuthal * CPUOVERLAP,
			  MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req4);
	    }
	} else {
	    if (CPU_Rank != CPU_Highest) {
		MPI_Irecv(RecvOuterBoundary, NAzimuthal * CPUOVERLAP,
			  MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req3);
		MPI_Isend(SendOuterBoundary, NAzimuthal * CPUOVERLAP,
			  MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req4);
	    }
	    if (CPU_Rank != 0) {
		MPI_Irecv(RecvInnerBoundary, NAzimuthal * CPUOVERLAP,
			  MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
		MPI_Isend(SendInnerBoundary, NAzimuthal * CPUOVERLAP,
			  MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req2);
	    }
	}

	if (CPU_Rank != 0) {
	    MPI_Wait(&req1, &global_MPI_Status);
	    MPI_Wait(&req2, &global_MPI_Status);
	    memcpy(Temperature.Field, RecvInnerBoundary, l * sizeof(double));
	}

	if (CPU_Rank != CPU_Highest) {
	    MPI_Wait(&req3, &global_MPI_Status);
	    MPI_Wait(&req4, &global_MPI_Status);
	    memcpy(Temperature.Field + oo, RecvOuterBoundary,
		   l * sizeof(double));
	}
    }

    if (iterations == parameters::radiative_diffusion_max_iterations) {
	logging::print_master(
	    LOG_WARNING
	    "Maximum iterations (%u) reached in radiative_diffusion (omega = %lg). Norm is %lg with a last change of %lg.\n",
	    parameters::radiative_diffusion_max_iterations, omega,
	    absolute_norm, norm_change);
    }

    // adapt omega
    if (old_iterations < iterations) {
	direction *= -1;
    }

    if (parameters::radiative_diffusion_omega_auto_enabled) {
	omega += direction * 0.01;
    }

    if (omega >= 2.0) {
	omega = 1.99;
	direction = -1;
    }

    if (omega <= 1.0) {
	omega = 1.0;
	direction = 1;
    }

    old_iterations = iterations;

    logging::print_master(LOG_VERBOSE "%u iterations, omega=%lf\n", iterations,
			  omega);

    // compute energy from temperature
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nr - 1; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		Energy(nr, naz) = Temperature(nr, naz) *
						Sigma(nr, naz) / (parameters::ADIABATICINDEX - 1.0) /
					    parameters::MU * constants::R;
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
		const double h0 = parameters::ASPECTRATIO_REF;
		const double beta = parameters::FLARINGINDEX;
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
		parameters::ASPECTRATIO_REF * std::pow(dist, parameters::FLARINGINDEX) * std::sqrt(Cs2);
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

		Cs2 += (parameters::ASPECTRATIO_REF * parameters::ASPECTRATIO_REF * 
				std::pow(dist, 2.0*parameters::FLARINGINDEX)
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
	switch (parameters::ASPECTRATIO_MODE) {
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

/**
	computes aspect ratio
*/
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

/**
	computes aspect ratio for an entire Nbody system
*/
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

/**
	computes aspect ratio with respect to the center of mass
*/
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

void compute_scale_height(t_data &data, const double current_time)
{
    switch (parameters::ASPECTRATIO_MODE) {
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
}

/**
	computes pressure
*/
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

/**
	computes temperature
*/
void compute_temperature(t_data &data)
{

	const unsigned int Nr = data[t_data::TEMPERATURE].get_size_radial();
	const unsigned int Nphi = data[t_data::TEMPERATURE].get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {
	    if (parameters::Adiabatic) {
		const double mu = pvte::get_mu(data, nr, naz);
		const double gamma_eff =
			pvte::get_gamma_eff(data, nr, naz);

		data[t_data::TEMPERATURE](nr, naz) =
		    mu / constants::R * (gamma_eff - 1.0) *
			data[t_data::ENERGY](nr, naz) /
			data[t_data::SIGMA](nr, naz);
	    } else if (parameters::Polytropic) {
		const double mu = pvte::get_mu(data, nr, naz);
		const double gamma_eff =
			pvte::get_gamma_eff(data, nr, naz);
		data[t_data::TEMPERATURE](nr, naz) =
		    mu / constants::R * parameters::POLYTROPIC_CONSTANT *
			std::pow(data[t_data::SIGMA](nr, naz),
			     gamma_eff - 1.0);
	    } else { // Isothermal
		data[t_data::TEMPERATURE](nr, naz) =
		    parameters::MU / constants::R *
			data[t_data::PRESSURE](nr, naz) /
			data[t_data::SIGMA](nr, naz);
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
