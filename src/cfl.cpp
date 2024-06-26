#include "cfl.h"
#include "global.h"
#include "logging.h"
#include "parameters.h"
#include "simulation.h"
#include "constants.h"
#include "pvte_law.h"

namespace cfl {

static std::vector<double> v_mean;
static std::vector<double> v_residual;

void init(t_data &data) {
	const t_polargrid &v_azimuthal = data[t_data::V_AZIMUTHAL];
	v_mean.resize(v_azimuthal.get_size_radial());
	v_residual.resize(v_azimuthal.get_size_radial()*v_azimuthal.get_size_azimuthal());
}


static void print_info()
{
	logging::print_master(LOG_INFO
			  "\nInteractive status requested with SIGUSR1\n");
	logging::print_master(LOG_INFO "hydro dt = %g\n", sim::last_dt);
	logging::print_master(LOG_INFO "output number = %d\n", sim::N_snapshot);
	logging::print_master(LOG_INFO "outer loop iteration = %d\n", sim::N_monitor);
	logging::print_master(LOG_INFO "N hydro step = %d\n", sim::N_hydro_iter);
	logging::print_master(LOG_INFO "time = %g\n", sim::time);
}

static inline unsigned int cell_number(const unsigned int nrad, 
								const unsigned int naz, 
								const unsigned int Naz_tot) {
    return nrad * Naz_tot + naz;
}

static void timestep_debug_report(t_data &data,
								  const unsigned int n_radial_debug,
								  const unsigned int n_azimuthal_debug,
								  const double dt_global,
								  const double dt_cell,
								  const double invdt1,
								  const double invdt2,
								  const double invdt3,
								  const double invdt4,
								  const double invdt5,
								  const double invdt6,
								  const double dt_shear,
								  const double dt_stable_visc
								  ){

	PRINT_SIG_INFO = false;

	const unsigned int n_radial = n_radial_debug;
	const unsigned int n_azimuthal = n_azimuthal_debug;

	// debugging variables
	double viscRadial = 0.0, viscAzimuthal = 0.0;
	double itdbg1 = std::numeric_limits<double>::max(), itdbg2 = std::numeric_limits<double>::max(), itdbg3 = std::numeric_limits<double>::max(),
	   itdbg4 = std::numeric_limits<double>::max(), itdbg5 = std::numeric_limits<double>::max(), itdbg6 = std::numeric_limits<double>::max();

	if (invdt1 != 0) {
		itdbg1 = 1.0 / invdt1 * parameters::CFL;
	}
	if (invdt2 != 0) {
		itdbg2 = std::fabs(1.0 / invdt2) * parameters::CFL;
	}
	if (invdt3 != 0) {
		itdbg3 = std::fabs(1.0 / invdt3) * parameters::CFL;
	}
	if (invdt4 != 0) {
		itdbg4 = 1.0 / invdt4 * parameters::CFL;
	}
	if (invdt5 != 0) {
		itdbg5 = 1.0 / invdt5 * parameters::CFL;
	}
	if (invdt6 != 0) {
		itdbg6 = 1.0 / invdt6 * parameters::CFL;
	}
	if ((parameters::artificial_viscosity ==
		 parameters::artificial_viscosity_SN) &&
			(parameters::artificial_viscosity_factor > 0)) {


		const t_polargrid &v_radial = data[t_data::V_RADIAL];
		const t_polargrid &v_azimuthal = data[t_data::V_AZIMUTHAL];

		// velocity differences in radial & azimuthal direction
		double dvRadial = v_radial(n_radial + 1, n_azimuthal) -
				v_radial(n_radial, n_azimuthal);
		double dvAzimuthal =
				v_azimuthal(n_radial,
							n_azimuthal == v_radial.get_max_azimuthal()
							? 0
							: n_azimuthal + 1) - v_azimuthal(n_radial, n_azimuthal);

		if (dvRadial >= 0.0) {
			dvRadial = std::numeric_limits<double>::min();
		} else {
			dvRadial = -dvRadial;
		}

		if (dvAzimuthal >= 0.0) {
			dvAzimuthal = std::numeric_limits<double>::min();
		} else {
			dvAzimuthal = -dvAzimuthal;
		}

		// cell sizes in radial & azimuthal direction
		double dxRadial = Rsup[n_radial] - Rinf[n_radial];
		double dxAzimuthal = Rmed[n_radial] * 2.0 * M_PI /
					 (double)(v_radial.get_size_azimuthal());

		viscRadial =
				dxRadial / dvRadial / 4.0 /
				std::pow(parameters::artificial_viscosity_factor,
						 2) * parameters::CFL;
		viscAzimuthal =
				dxAzimuthal / dvAzimuthal / 4.0 /
				std::pow(parameters::artificial_viscosity_factor,
						 2) * parameters::CFL;
	}

	print_info();

	logging::print(LOG_INFO
				   "Timestep control information for CPU %d: \n",
				   CPU_Rank);
	logging::print(
				LOG_INFO
				"Most restrictive cell at nRadial=%d and nAzimuthal=%d\n",
				IMIN + n_radial_debug, n_azimuthal_debug);
	logging::print(LOG_INFO "located at radius Rmed         : %g\n",
				   Rmed[n_radial_debug]);
	const double SigmaFloor =
			parameters::sigma0 * parameters::sigma_floor;
	logging::print(LOG_INFO "Cell has a surface denstity of : %g\t%g g/cm^2\t%g 1/SigmaFloor\n",
				   data[t_data::SIGMA](n_radial_debug, n_azimuthal_debug), data[t_data::SIGMA](n_radial_debug, n_azimuthal_debug)*units::surface_density,
			data[t_data::SIGMA](n_radial_debug, n_azimuthal_debug)/SigmaFloor);
	if (parameters::Adiabatic) {
		const double mu = pvte::get_mu(data, n_radial_debug, n_azimuthal_debug);
		const double gamma_eff =
				pvte::get_gamma_eff(data, n_radial_debug, n_azimuthal_debug);
		const double cell_temperature =
				mu / constants::R * (gamma_eff - 1.0) *
				data[t_data::ENERGY](n_radial_debug, n_azimuthal_debug) /
				data[t_data::SIGMA](n_radial_debug, n_azimuthal_debug) *
				units::temperature.get_code_to_cgs_factor();
		logging::print(LOG_INFO "Cell has a Temperature of      : %g K\n",
					   cell_temperature);
	}
	logging::print(LOG_INFO "Sound speed limit              : %g\n",
				   itdbg1);
	logging::print(LOG_INFO "Radial motion limit            : %g\n",
				   itdbg2);
	logging::print(LOG_INFO "Residual circular motion limit : %g\n",
				   itdbg3);

	logging::print(LOG_INFO "Articifial Viscosity limit     : %g\n",
					   itdbg4);
	if (parameters::artificial_viscosity == parameters::artificial_viscosity_SN) {
		logging::print(LOG_INFO "   Arise from r with limit     : %g\n",
					   viscRadial);
		logging::print(LOG_INFO "   and from theta with limit   : %g\n",
					   viscAzimuthal);
	}

	if (StabilizeViscosity == 2) {
		logging::print(LOG_INFO "Viscosity stability limit %g\n", dt_stable_visc);
	}

	logging::print(LOG_INFO "Fargo shear limit		   : %g\n", dt_shear);
	logging::print(LOG_INFO "Kinematic viscosity limit      : %g\n",
				   itdbg5);
	logging::print(LOG_INFO "Heating cooling limit		: %g Qm=%.3g Qp=%.3g\n",
				   itdbg6, data[t_data::QMINUS](n_radial, n_azimuthal), data[t_data::QPLUS](n_radial, n_azimuthal));
	logging::print(LOG_INFO "Limit time step for this cell  : %g\n",
				   dt_cell);
	logging::print(LOG_INFO "Limit time step adopted        : %g\n",
				   dt_global);
}


double condition_cfl(t_data &data, const double dt_global_input)
{

	const t_polargrid &v_radial = data[t_data::V_RADIAL];
	const t_polargrid &v_azimuthal = data[t_data::V_AZIMUTHAL];
	const t_polargrid &soundspeed = data[t_data::SOUNDSPEED];


	dt_parabolic_local = std::numeric_limits<double>::max();

	// Calculate and fill VMean array
	#pragma omp parallel for
	for (unsigned int n_radial = 0; n_radial < v_azimuthal.get_size_radial();
	 ++n_radial) {
	v_mean[n_radial] = 0.0;
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < v_azimuthal.get_size_azimuthal(); ++n_azimuthal) {
		v_mean[n_radial] += v_azimuthal(n_radial, n_azimuthal);
	}
	v_mean[n_radial] /= (double)(v_azimuthal.get_size_azimuthal());
	}

	const double denom = fabs(v_mean[0]*InvRmed[0] - v_mean[1]*InvRmed[1]) + 1.0e-100;
	double dt_core = parameters::CFL * dphi / denom;

	#pragma omp parallel for reduction(min : dt_core, dt_parabolic_local)
	for (unsigned int nr = radial_first_active;
	 nr < radial_active_size; ++nr) {
		// FARGO algorithm timestep criterion. See Masset 2000 Sect. 3.3.
		const double denom = fabs(v_mean[nr] * InvRmed[nr] -
							 v_mean[nr + 1] * InvRmed[nr + 1]) + 1.0e-100;
		double shear_dt = parameters::CFL * dphi / denom;

		if (shear_dt < dt_core){
			dt_core = shear_dt;
		}


	// cell sizes in radial & azimuthal direction
	const double dxRadial = Rsup[nr] - Rinf[nr];
	const double dxAzimuthal = Rmed[nr] * dphi;
	const double cell_size = std::min(dxRadial, dxAzimuthal);

	for (unsigned int naz = 0;
		 naz < v_radial.get_size_azimuthal(); ++naz) {
		if (parameters::fast_transport) {
		// FARGO algorithm
		v_residual[cell_number(nr, naz, v_azimuthal.get_size_azimuthal())] =
			v_azimuthal(nr, naz) - v_mean[nr];
		} else {
		// Standard algorithm
		v_residual[cell_number(nr, naz, v_azimuthal.get_size_azimuthal())] = v_azimuthal(nr, naz);
		}
	}

	for (unsigned int naz = 0;
		 naz < v_radial.get_size_azimuthal(); ++naz) {

		// sound speed limit
		const double invdt1 = soundspeed(nr, naz) / cell_size;

		// radial motion limit
		// we do not need abs() because only square of it is used later
		const double invdt2 = v_radial(nr, naz) / dxRadial;

		// residual circular motion limit
		// we do not need abs() because only square of it is used later
		const double invdt3 = v_residual[cell_number(nr, naz, v_azimuthal.get_size_azimuthal())] / dxAzimuthal;

		double leapfrog_cfl_factor;
		if (parameters::hydro_integrator == LEAPFROG_INTEGRATOR){
			// since we are applying source terms twice with half the timestep,
			// we can consider that by allowing double the timestep size for source terms cfl
			leapfrog_cfl_factor = 0.6;
		} else {
			leapfrog_cfl_factor = 1.0;
		}

		// artificial viscosity limit
		double invdt4;
		if (parameters::artificial_viscosity ==
		parameters::artificial_viscosity_SN) {
			const unsigned int naz_next = naz == v_radial.get_max_azimuthal()
					? 0	: naz + 1;

			// velocity differences in radial & azimuthal direction
			double dvRadial = v_radial(nr + 1, naz) -
					v_radial(nr, naz);
			double dvAzimuthal =
					v_azimuthal(nr, naz_next) -
					v_azimuthal(nr, naz);

		if (dvRadial > 0.0) {
			dvRadial = 0.0;
		} else {
			dvRadial = -dvRadial;
		}

		if (dvAzimuthal > 0.0) {
			dvAzimuthal = 0.0;
		} else {
			dvAzimuthal = -dvAzimuthal;
		}

		invdt4 =
			4.0 * std::pow(parameters::artificial_viscosity_factor, 2) *
			std::max(dvRadial / dxRadial, dvAzimuthal / dxAzimuthal) * leapfrog_cfl_factor;
		} else { // TW artificial viscosity
			const unsigned int naz_next = naz == v_radial.get_max_azimuthal()
					? 0	: naz + 1;
			// div(v) = 1/r d(r v_r)/dr + 1/r d(v_phi)/dphi
			//  	 == d(v_r)/dr + 1/r [ d(v_phi)/dphi + v_r]
			const double eps_rr = (v_radial(nr+1, naz) - v_radial(nr, naz)) * InvDiffRsup[nr];
			const double eps_pp =  InvRb[nr] * ((v_azimuthal(nr, naz_next) - v_azimuthal(nr, naz)) * invdphi + 0.5*(v_radial(nr + 1, naz) + v_radial(nr, naz)));
			const double mdiv_V =  -std::min(eps_rr + eps_pp, 0.0);
			invdt4 = 4.0 * std::pow(parameters::artificial_viscosity_factor, 2)  * mdiv_V * leapfrog_cfl_factor;
		}

		// kinematic viscosity limit
		// for diffusion stability: dt < dx^2 / (4 * nu)
		const double invdt5 = 4.0 * data[t_data::VISCOSITY](nr, naz) / std::pow(cell_size, 2)
				* leapfrog_cfl_factor;

		// heating / cooling limit
		double invdt6;
		if (parameters::Adiabatic) {
		// Limit energy update from heating / cooling to given fraction
		// per dt
		const double inv_limit =
			1.0 / parameters::HEATING_COOLING_CFL_LIMIT;
		const double Qp = data[t_data::QPLUS](nr, naz);
		const double Qm = data[t_data::QMINUS](nr, naz);
		const double E = data[t_data::ENERGY](nr, naz);
		invdt6 = inv_limit * std::fabs((Qp - Qm) / E) * leapfrog_cfl_factor;
		} else {
		invdt6 = 0.0;
		}

		double dt_cell;
		// calculate new dt based on different limits
		dt_cell = parameters::CFL /
			  std::sqrt(std::pow(invdt1, 2) + std::pow(invdt2, 2) +
					std::pow(invdt3, 2) + std::pow(invdt4, 2) +
					std::pow(invdt5, 2) + std::pow(invdt6, 2));


		double dt_stable_visc = std::numeric_limits<double>::max();
		if (StabilizeViscosity == 2) {
		const double cphi =
			data[t_data::VISCOSITY_CORRECTION_FACTOR_PHI](nr,
								  naz);
		const double cr = data[t_data::VISCOSITY_CORRECTION_FACTOR_R](
			nr, naz);
		const double c =
			std::min(cphi, cr); // c < 0.0 is negative, so take min to
					// get 'larger' negative number

		if(c*dt_cell*parameters::CFL < -1.0){
		logging::print(LOG_INFO "Viscosity unstable with c=%.5e\n", c);
		}

		if (c != 0.0) {
			const double dtStable = -parameters::CFL / c;
			dt_stable_visc = dtStable;
			dt_cell = std::min(dt_cell, dtStable);
		}
		}

		if (dt_cell < dt_core) {
		dt_core = dt_cell;
		}

		if(PRINT_SIG_INFO){
		if(dt_core == dt_global_input){
			timestep_debug_report(data,
								  nr,
								  naz,
								  dt_global_input,
								  dt_cell,
								  invdt1,
								  invdt2,
								  invdt3,
								  invdt4,
								  invdt5,
								  invdt6,
								  shear_dt,
								  dt_stable_visc);
		}
		}

	}
	}

	double dt_global;
	MPI_Allreduce(&dt_core, &dt_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	return dt_global;
}


} // close namespace cfl
