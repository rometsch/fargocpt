#include "cfl.h"
#include "global.h"
#include "logging.h"
#include "parameters.h"
#include "simulation.h"
#include "constants.h"
#include "pvte_law.h"

namespace cfl {

/**
	\param VRadial radial velocity polar grid
	\param VAzimuthal azimuthal velocity polar grid
	\param SoundSpeed sound speed polar grid
*/
static void print_info()
{
	logging::print_master(LOG_INFO
			  "\nInteractive status requested with SIGUSR1\n");
	logging::print_master(LOG_INFO "hydro dt = %g\n", sim::last_dt);
	logging::print_master(LOG_INFO "output number = %d\n", sim::N_output);
	logging::print_master(LOG_INFO "outer loop iteration = %d\n", sim::N_outer_loop);
	logging::print_master(LOG_INFO "N hydro step = %d\n", sim::N_hydro_iter);
	logging::print_master(LOG_INFO "PhysicalTime = %g\n", sim::PhysicalTime);
}

double condition_cfl(t_data &data)
{
    t_polargrid &v_radial = data[t_data::V_RADIAL];
    t_polargrid &v_azimuthal = data[t_data::V_AZIMUTHAL];
	t_polargrid &soundspeed = data[t_data::SOUNDSPEED];
    
	dt_parabolic_local = std::numeric_limits<double>::max();
	std::vector<double> v_mean(v_radial.get_size_radial());
	std::vector<double> v_residual(v_radial.get_size_azimuthal());
	double dt_core = std::numeric_limits<double>::max();
    double dt_cell;

	// debugging variables
	double viscRadial = 0.0, viscAzimuthal = 0.0;
	unsigned int n_azimuthal_debug = 0, n_radial_debug = 0;
	double itdbg1 = std::numeric_limits<double>::max();
    double itdbg2 = std::numeric_limits<double>::max();
    double itdbg3 = std::numeric_limits<double>::max();
    double itdbg4 = std::numeric_limits<double>::max();
    double itdbg5 = std::numeric_limits<double>::max();
    double itdbg6 = std::numeric_limits<double>::max();

	// Calculate and fill VMean array
	for (unsigned int n_radial = 0; n_radial < v_azimuthal.get_size_radial();
	 ++n_radial) {
	v_mean[n_radial] = 0.0;
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < v_azimuthal.get_size_azimuthal(); ++n_azimuthal) {
		v_mean[n_radial] += v_azimuthal(n_radial, n_azimuthal);
	}
	v_mean[n_radial] /= (double)(v_azimuthal.get_size_azimuthal());
	}

	for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size; ++n_radial) {
	// cell sizes in radial & azimuthal direction
	double dxRadial = Rsup[n_radial] - Rinf[n_radial];
	double dxAzimuthal = Rmed[n_radial] * 2.0 * M_PI /
				 (double)(v_radial.get_size_azimuthal());

	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < v_radial.get_size_azimuthal(); ++n_azimuthal) {
		if (parameters::fast_transport) {
		// FARGO algorithm
		v_residual[n_azimuthal] =
			v_azimuthal(n_radial, n_azimuthal) - v_mean[n_radial];
		} else {
		// Standard algorithm
		v_residual[n_azimuthal] = v_azimuthal(n_radial, n_azimuthal);
		}
	}

	// there is no v_residual[v_radial.Nsec]
	// v_residual[v_radial.Nsec]=v_residual[0];

	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < v_radial.get_size_azimuthal(); ++n_azimuthal) {
		double invdt1, invdt2, invdt3, invdt4, invdt5, invdt6;

		// velocity differences in radial & azimuthal direction
		double dvRadial = v_radial(n_radial + 1, n_azimuthal) -
				  v_radial(n_radial, n_azimuthal);
		double dvAzimuthal =
		v_azimuthal(n_radial,
				n_azimuthal == v_radial.get_max_azimuthal()
				? 0
				: n_azimuthal + 1) -
		v_azimuthal(n_radial, n_azimuthal);

		// sound speed limit
		invdt1 = soundspeed(n_radial, n_azimuthal) /
			 (std::min(dxRadial, dxAzimuthal));

		// radial motion limit
		invdt2 = fabs(v_radial(n_radial, n_azimuthal)) / dxRadial;

		// residual circular motion limit
		invdt3 = fabs(v_residual[n_azimuthal]) / dxAzimuthal;

		// artificial viscosity limit
		if (parameters::artificial_viscosity ==
		parameters::artificial_viscosity_SN) {
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

		invdt4 =
			4.0 * std::pow(parameters::artificial_viscosity_factor, 2) *
			std::max(dvRadial / dxRadial, dvAzimuthal / dxAzimuthal);
		} else {
		invdt4 = 0.0;
		}

		// kinematic viscosity limit
		invdt5 = 4.0 * data[t_data::VISCOSITY](n_radial, n_azimuthal) *
			 std::max(1 / std::pow(dxRadial, 2),
				  1 / std::pow(dxAzimuthal, 2));

		// heating / cooling limit
		if (parameters::Adiabatic) {
		// Limit energy update from heating / cooling to given fraction
		// per dt
		const double inv_limit =
			1.0 / parameters::HEATING_COOLING_CFL_LIMIT;
		const double Qp = data[t_data::QPLUS](n_radial, n_azimuthal);
		const double Qm = data[t_data::QMINUS](n_radial, n_azimuthal);
		const double E = data[t_data::ENERGY](n_radial, n_azimuthal);
		invdt6 = inv_limit * std::fabs((Qp - Qm) / E);
		} else {
		invdt6 = 0.0;
		}

	    if (parameters::EXPLICIT_VISCOSITY) {
		// calculate new dt based on different limits
		dt_cell = parameters::CFL /
			  std::sqrt(std::pow(invdt1, 2) + std::pow(invdt2, 2) +
					std::pow(invdt3, 2) + std::pow(invdt4, 2) +
					std::pow(invdt5, 2) + std::pow(invdt6, 2));
		} else {
		// viscous timestep
		if (invdt4 > 0.0 && invdt5 > 0.0) {
			dt_parabolic_local = std::min(
			dt_parabolic_local,
			parameters::CFL / std::sqrt(std::pow(invdt4, 2) +
							std::pow(invdt5, 2)));
		}

		// calculate new dt based on different limits
		dt_cell = parameters::CFL /
			  std::sqrt(std::pow(invdt1, 2) + std::pow(invdt2, 2) +
					std::pow(invdt3, 2) + std::pow(invdt6, 2));

		dt_cell = std::min(dt_cell, 3.0 * dt_parabolic_local);
		}

		if (StabilizeViscosity == 2) {
		const double cphi =
			data[t_data::VISCOSITY_CORRECTION_FACTOR_PHI](n_radial,
								  n_azimuthal);
		const double cr = data[t_data::VISCOSITY_CORRECTION_FACTOR_R](
			n_radial, n_azimuthal);
		const double c =
			std::min(cphi, cr); // c < 0.0 is negative, so take min to
					// get 'larger' negative number

		if (c != 0.0) {
			const double dtStable = -parameters::CFL / c;
			dt_cell = std::min(dt_cell, dtStable);
		}
		}

		if (dt_cell < dt_core) {
		dt_core = dt_cell;

		if (PRINT_SIG_INFO) {
			n_radial_debug = n_radial;
			n_azimuthal_debug = n_azimuthal;
			if (invdt1 != 0) {
			itdbg1 = 1.0 / invdt1;
			}
			if (invdt2 != 0) {
			itdbg2 = 1.0 / invdt2;
			}
			if (invdt3 != 0) {
			itdbg3 = 1.0 / invdt3;
			}
			if (invdt4 != 0) {
			itdbg4 = 1.0 / invdt4;
			}
			if (invdt5 != 0) {
			itdbg5 = 1.0 / invdt5;
			}
			if (invdt6 != 0) {
			itdbg6 = 1.0 / invdt6;
			}
			if ((parameters::artificial_viscosity ==
			 parameters::artificial_viscosity_SN) &&
			(parameters::artificial_viscosity_factor > 0)) {
			viscRadial =
				dxRadial / dvRadial / 4.0 /
				std::pow(parameters::artificial_viscosity_factor,
					 2);
			viscAzimuthal =
				dxAzimuthal / dvAzimuthal / 4.0 /
				std::pow(parameters::artificial_viscosity_factor,
					 2);
			}
		}
		}
	}
	}

	// FARGO algorithm timestep criterion. See Masset 2000 Sect. 3.3.
	for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size - 1; ++n_radial) {
	const double azimuthal_cell_size = 2.0 * M_PI / (double)NAzimuthal;
	const double shear_dt =
		parameters::CFL * azimuthal_cell_size /
		fabs(v_mean[n_radial] * InvRmed[n_radial] -
		 v_mean[n_radial + 1] * InvRmed[n_radial + 1]);

	if (shear_dt < dt_core)
		dt_core = shear_dt;
	}

	double dt_global;
	MPI_Allreduce(&dt_core, &dt_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

	if (PRINT_SIG_INFO) {
	print_info();

	if (dt_core == dt_global) {

		logging::print(LOG_INFO
			   "Timestep control information for CPU %d: \n",
			   CPU_Rank);
		logging::print(
		LOG_INFO
		"Most restrictive cell at nRadial=%d and nAzimuthal=%d\n",
		n_radial_debug, n_azimuthal_debug);
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
			units::temperature.get_cgs_factor();
		logging::print(LOG_INFO "Cell has a Temperature of      : %g K\n",
			   cell_temperature);
		}
		logging::print(LOG_INFO "Sound speed limit              : %g\n",
			   itdbg1);
		logging::print(LOG_INFO "Radial motion limit            : %g\n",
			   itdbg2);
		logging::print(LOG_INFO "Residual circular motion limit : %g\n",
			   itdbg3);

		if (parameters::artificial_viscosity_factor > 0) {
		logging::print(LOG_INFO "Articifial Viscosity limit     : %g\n",
				   itdbg4);
		logging::print(LOG_INFO "   Arise from r with limit     : %g\n",
				   viscRadial);
		logging::print(LOG_INFO "   and from theta with limit   : %g\n",
				   viscAzimuthal);
		} else {
		logging::print(LOG_INFO
				   "Articifial Viscosity limit     : disabled\n");
		}
		logging::print(LOG_INFO "Kinematic viscosity limit      : %g\n",
			   itdbg5);
		logging::print(LOG_INFO "Heating cooling limit      : %g\n",
			   itdbg6);
		logging::print(LOG_INFO "Limit time step for this cell  : %g\n",
			   dt_core);
		logging::print(LOG_INFO "Limit time step adopted        : %g\n",
			   dt_global);
	}

	PRINT_SIG_INFO = 0;
	}
    return dt_global;
}


} // close namespace cfl