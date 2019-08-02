/**
	\file SourceEuler.c

	Contains routines used by the hydrodynamical loop. More specifically, it contains the main loop itself and all the source term substeps (with the exception of the evaluation of the viscous force). The transport substep is treated elsewhere.
*/

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <limits.h>

#include "output.h"
#include "viscosity.h"
#include "SourceEuler.h"
#include "Pframeforce.h"
#include "stress.h"
#include "SideEuler.h"
#include "TransportEuler.h"
#include "LowTasks.h"
#include "selfgravity.h"
#include "commbound.h"
#include "constants.h"
#include "Pframeforce.h"
#include "logging.h"
#include "global.h"
#include "nongnu.h"
#include "util.h"
#include "Theo.h"
#include "opacity.h"
#include "units.h"
#include "parameters.h"
#include "boundary_conditions.h"
#include "quantities.h"
#include "particles.h"
#include <cstring>
extern boolean Corotating;

extern boolean FastTransport;
Pair DiskOnPrimaryAcceleration;
double dtemp;

/**
	Checks polargrid for negative entries.

	\param array polargrid to check
	\returns >0 if there are negative entries, 0 otherwise
*/
int DetectCrash(t_polargrid* array)
{
	unsigned int result = 0;

	for (unsigned int n_radial = 0; n_radial < array->Nrad; ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal < array->Nsec; ++n_azimuthal) {
			if ((*array)(n_radial,n_azimuthal) < 0.0) {
				logging::print(LOG_WARNING "%s negative in cell: (%u,%u)=%g\n",array->get_name(),n_radial,n_azimuthal,(*array)(n_radial,n_azimuthal));
				result += 1;
			}
		}
	}

	return result;
}

/**
	Assures miminum value in each cell.

	\param dst polar grid
	\param minimum_value minimum value
*/
bool assure_minimum_value(t_polargrid &dst, double minimum_value)
{
	bool found = false;
	bool is_dens = strcmp(dst.get_name(), "dens") == 0;

	for (unsigned int n_radial = 0; n_radial <= dst.get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= dst.get_max_azimuthal(); ++n_azimuthal) {
			if (dst(n_radial, n_azimuthal) < minimum_value) {
				if (is_dens) {
					double mass_delta = (minimum_value - dst(n_radial, n_azimuthal))*Surf[n_radial];
					sum_without_ghost_cells(MassDelta.FloorPositive, mass_delta, n_radial);
				}
				dst(n_radial, n_azimuthal) = minimum_value;
#ifndef NDEBUG
				logging::print(LOG_DEBUG "assure_minimum_value: %s(%u,%u)=%g < %g\n",dst.get_name(),n_radial,n_azimuthal,dst(n_radial,n_azimuthal),minimum_value);
#endif
				found = true;
			}
		}
	}

	return found;
}

/**
	Assures a miminum temperature in each cell.

	\param energy energy polar grid
	\param minimum_value minimum temperature
*/
bool assure_minimum_temperature(t_polargrid &energy, t_polargrid &density, double minimum_value)
{
	bool found = false;

	for (unsigned int n_radial = 0; n_radial <= energy.get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= energy.get_max_azimuthal(); ++n_azimuthal) {
			if (energy(n_radial, n_azimuthal) < minimum_value*density(n_radial, n_azimuthal)/parameters::MU*constants::R/(ADIABATICINDEX-1.0)) {
#ifndef NDEBUG
				logging::print(LOG_DEBUG "assure_minimum_temperature: (%u,%u)=%g<%g\n",n_radial, n_azimuthal, energy(n_radial, n_azimuthal)*units::temperature.get_cgs_factor()/density(n_radial,n_azimuthal)*parameters::MU/constants::R*(ADIABATICINDEX-1.0),minimum_value*units::temperature.get_cgs_factor(),minimum_value);
#endif
				energy(n_radial, n_azimuthal) = minimum_value*density(n_radial, n_azimuthal)/parameters::MU*constants::R/(ADIABATICINDEX-1.0);
				found = true;
			}
		}
	}

	return found;
}

/**
	Assures a miminum temperature in each cell.

	\param energy energy polar grid
	\param minimum_value minimum temperature
*/
bool assure_maximum_temperature(t_polargrid &energy, t_polargrid &density, double maximum_value)
{
	if (isnan(maximum_value))
		return false;

	bool found = false;

	for (unsigned int n_radial = 0; n_radial <= energy.get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= energy.get_max_azimuthal(); ++n_azimuthal) {
			if (energy(n_radial, n_azimuthal) > maximum_value*density(n_radial, n_azimuthal)/parameters::MU*constants::R/(ADIABATICINDEX-1.0)) {
#ifndef NDEBUG
				logging::print(LOG_DEBUG "assure_maximum_temperature: (%u,%u)=%g>%g\n",n_radial, n_azimuthal, energy(n_radial, n_azimuthal)*units::temperature.get_cgs_factor()/density(n_radial,n_azimuthal)*parameters::MU/constants::R*(ADIABATICINDEX-1.0),maximum_value*units::temperature.get_cgs_factor(),maximum_value);
#endif
				energy(n_radial, n_azimuthal) = maximum_value*density(n_radial, n_azimuthal)/parameters::MU*constants::R/(ADIABATICINDEX-1.0);
				found = true;
			}
		}
	}

	return found;
}

void recalculate_derived_disk_quantities(t_data &data, bool force_update)
{

	if(parameters::Locally_Isothermal)
	{
		compute_pressure(data, force_update);
	}


	if (parameters::Adiabatic || parameters::Polytropic) {
		compute_temperature(data, force_update);
		compute_sound_speed(data, force_update);
		compute_aspect_ratio(data, force_update);
		compute_pressure(data, force_update);
	}

	viscosity::update_viscosity(data);
}

void init_euler(t_data &data)
{
	InitTransport();
	InitComputeAccel();

	if(parameters::Locally_Isothermal)
	{
		compute_sound_speed(data, false);
		compute_pressure(data, false);
		compute_temperature(data, false);
		compute_aspect_ratio(data, false);
	}

	if (parameters::Adiabatic || parameters::Polytropic) {
		compute_temperature(data, false);
		compute_sound_speed(data, false);
		compute_aspect_ratio(data, false);
		compute_pressure(data, false);
	}

	viscosity::update_viscosity(data);
}


/**

*/
void FreeEuler()
{
		FreeTransport();
}

/**
	copy one polar grid into another

	\param dst destination polar grid
	\param src source polar grid
*/
void copy_polargrid(t_polargrid &dst, t_polargrid &src)
{
	assert((dst.get_size_radial() == src.get_size_radial()) && (dst.get_size_azimuthal() == src.get_size_azimuthal()));

	for (unsigned int n_radial = 0; n_radial <= dst.get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= dst.get_max_azimuthal(); ++n_azimuthal) {
			dst(n_radial, n_azimuthal) = src(n_radial, n_azimuthal);
		}
	}
}

/**
	switches polar grids

	\param dst destination polar grid
	\param src source polar grid
*/
void SwitchPolarGrid(t_polargrid *dst, t_polargrid *src)
{
	double *tmp;

	assert(dst->Nsec = src->Nsec);
	assert(dst->Nrad = dst->Nrad);

	tmp = dst->Field;
	dst->Field = src->Field;
	src->Field = tmp;
}

/**

	\param force
	\param data
	\param sys
*/
void AlgoGas(unsigned int nTimeStep, Force* force, t_data &data)
{
	double local_gas_time_step_cfl = 1.0;
	double global_gas_time_step_cfl;

	double dt;
	double OmegaNew, domega;
	double planet0_old_x = 0.0, planet0_old_y = 0.0;

    dtemp=0.0;

	if (parameters::calculate_disk) {
		CommunicateBoundaries(&data[t_data::DENSITY], &data[t_data::V_RADIAL], &data[t_data::V_AZIMUTHAL], &data[t_data::ENERGY]);
		local_gas_time_step_cfl = condition_cfl(data, data[t_data::V_RADIAL], data[t_data::V_AZIMUTHAL], data[t_data::SOUNDSPEED], DT-dtemp);
	}

	MPI_Allreduce(&local_gas_time_step_cfl, &global_gas_time_step_cfl, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	dt = DT/global_gas_time_step_cfl;
	boundary_conditions::apply_boundary_condition(data, dt, false);

	if (data[t_data::ALPHA_GRAV_MEAN].get_write()) {
		quantities::calculate_alpha_grav_mean_reset(data);
	}
	if (data[t_data::ALPHA_GRAV_MEAN_1D].get_write()) {
		quantities::calculate_radial_alpha_grav_mean_reset(data);
	}
	if (data[t_data::ALPHA_REYNOLDS_MEAN].get_write()) {
		quantities::calculate_alpha_reynolds_mean_reset(data);
	}
	if (data[t_data::ALPHA_REYNOLDS_MEAN_1D].get_write()) {
		quantities::calculate_radial_alpha_reynolds_mean_reset(data);
	}

	while (dtemp < DT) {
		logging::print_master(LOG_VERBOSE "AlgoGas: Total: %*i/%i (%5.2f %%) - Timestep: %#7f/%#7f (%5.2f %%)\n",(int)ceil(log10(NTOT)),nTimeStep,NTOT,(double)nTimeStep/(double)NTOT*100.0,dtemp,DT,dtemp/DT*100.0);

        dtemp += dt;

		DiskOnPrimaryAcceleration.x = 0.0;
		DiskOnPrimaryAcceleration.y = 0.0;

		if (Corotating == YES) {
			// save old planet positions
			planet0_old_x = data.get_planetary_system().get_planet(0).get_x();
			planet0_old_y = data.get_planetary_system().get_planet(0).get_y();
		}

		if (parameters::calculate_disk) {
			/* Indirect term star's potential computed here */
			if (parameters::disk_feedback)
				DiskOnPrimaryAcceleration = ComputeAccel(force, data, 0.0, 0.0, 0.0);

			/* Gravitational potential from star and planet(s) is computed and stored here*/
			FillForcesArrays(data);
			/* Planets' velocities are updated here from gravitationnal interaction with disk */
			AdvanceSystemFromDisk(force, data, dt);
		}

		/* Planets' positions and velocities are updated from gravitational interaction with star and other planets */
		if (parameters::integrate_planets) {
			AdvanceSystemRK5(data, dt);
		}

		if (parameters::integrate_particles) {
			particles::integrate(data, dt);
		}


		/* Below we correct v_azimuthal, planet's position and velocities if we work in a frame non-centered on the star */
		if (Corotating == YES) {
			double distance_new = sqrt( pow2(data.get_planetary_system().get_planet(0).get_x()) + pow2(data.get_planetary_system().get_planet(0).get_y()) );
			double distance_old = sqrt( pow2(planet0_old_x) + pow2(planet0_old_y) );
			double cross = planet0_old_x*data.get_planetary_system().get_planet(0).get_y()-data.get_planetary_system().get_planet(0).get_x()*planet0_old_y;

			// new = r_new x r_old = distance_new * distance_old * sin(alpha*dt)
			OmegaNew = asin(cross/(distance_new*distance_old))/dt;

			domega = (OmegaNew-OmegaFrame);
			if (parameters::calculate_disk)
				correct_v_azimuthal(data[t_data::V_AZIMUTHAL], domega);
			OmegaFrame = OmegaNew;
		}

		if (parameters::integrate_planets) {
			data.get_planetary_system().rotate(OmegaFrame*dt);
		}

		FrameAngle += OmegaFrame*dt;

		/* Now we update gas */
		if (parameters::calculate_disk) {

			if (DetectCrash(&data[t_data::DENSITY])) {
				logging::print(LOG_ERROR "DetectCrash: Density < 0\n");
				PersonalExit(1);
			}

			if (parameters::Adiabatic) {
				if (DetectCrash(&data[t_data::ENERGY])) {
					logging::print(LOG_ERROR "DetectCrash: Energy < 0\n");
					PersonalExit(1);
				}
			}

			SubStep1(data, dt);
			SubStep2(data, dt);

			// set new velocities
			copy_polargrid(data[t_data::V_RADIAL], data[t_data::V_RADIAL_NEW]);
			copy_polargrid(data[t_data::V_AZIMUTHAL], data[t_data::V_AZIMUTHAL_NEW]);

			boundary_conditions::apply_boundary_condition(data, dt, false);

			if (parameters::Adiabatic) {
				if ((parameters::artificial_viscosity == parameters::artificial_viscosity_TW) && (parameters::artificial_viscosity_dissipation)) {
					viscosity::compute_viscous_terms(data, true);
				} else {
					viscosity::compute_viscous_terms(data, false);
				}

				SubStep3(data, dt);

				copy_polargrid(data[t_data::ENERGY], data[t_data::ENERGY_NEW]);

				if (assure_minimum_temperature(data[t_data::ENERGY], data[t_data::DENSITY], parameters::minimum_temperature*units::temperature.get_inverse_cgs_factor())) {
					logging::print(LOG_DEBUG "Found temperature < %g %s after SubStep3.\n",parameters::minimum_temperature, units::temperature.get_cgs_symbol());
				}

				if (assure_maximum_temperature(data[t_data::ENERGY], data[t_data::DENSITY], parameters::maximum_temperature*units::temperature.get_inverse_cgs_factor())) {
					logging::print(LOG_DEBUG "Found temperature > %g %s after SubStep3.\n",parameters::maximum_temperature, units::temperature.get_cgs_symbol());
				}

				if (parameters::radiative_diffusion_enabled) {
					radiative_diffusion(data, dt);

					if (assure_minimum_temperature(data[t_data::ENERGY], data[t_data::DENSITY], parameters::minimum_temperature*units::temperature.get_inverse_cgs_factor())) {
						logging::print(LOG_DEBUG "Found temperature < %g %s after radiative_diffusion.\n",parameters::minimum_temperature, units::temperature.get_cgs_symbol());
					}

					if (assure_maximum_temperature(data[t_data::ENERGY], data[t_data::DENSITY], parameters::maximum_temperature*units::temperature.get_inverse_cgs_factor())) {
						logging::print(LOG_DEBUG "Found temperature > %g %s after radiative_diffusion.\n",parameters::maximum_temperature, units::temperature.get_cgs_symbol());
					}
				}


			}

			Transport(data, &data[t_data::DENSITY], &data[t_data::V_RADIAL], &data[t_data::V_AZIMUTHAL], &data[t_data::ENERGY], dt);

			// assure minimum density after all substeps & transport
			assure_minimum_value(data[t_data::DENSITY],parameters::sigma_floor*parameters::sigma0);

			if (parameters::Adiabatic) {
				// assure minimum temperature after all substeps & transport. it is crucial the check minimum density before!
				if (assure_minimum_temperature(data[t_data::ENERGY], data[t_data::DENSITY], parameters::minimum_temperature*units::temperature.get_inverse_cgs_factor())) {
					logging::print(LOG_DEBUG "Found temperature < %g %s after Transport.\n",parameters::minimum_temperature, units::temperature.get_cgs_symbol());
				}

				if (assure_maximum_temperature(data[t_data::ENERGY], data[t_data::DENSITY], parameters::maximum_temperature*units::temperature.get_inverse_cgs_factor())) {
					logging::print(LOG_DEBUG "Found temperature > %g %s after Transport.\n",parameters::maximum_temperature, units::temperature.get_cgs_symbol());
				}
			}
			boundary_conditions::apply_boundary_condition(data, dt, true);

			mdcp = CircumPlanetaryMass(data);
			exces_mdcp = mdcp - mdcp0;

			if (data[t_data::ALPHA_GRAV_MEAN].get_write()) {
				quantities::calculate_alpha_grav_mean_sumup(data, nTimeStep, dt);
			}
			if (data[t_data::ALPHA_GRAV_MEAN_1D].get_write()) {
				quantities::calculate_radial_alpha_grav_mean_sumup(data, nTimeStep, dt);
			}
			if (data[t_data::ALPHA_REYNOLDS_MEAN].get_write()) {
				quantities::calculate_alpha_reynolds_mean_sumup(data, nTimeStep, dt);
			}
			if (data[t_data::ALPHA_REYNOLDS_MEAN_1D].get_write()) {
				quantities::calculate_radial_alpha_reynolds_mean_sumup(data, nTimeStep, dt);
			}
		}

		PhysicalTime += dt;
		N_iter = N_iter +1;
		logging::print_runtime_info(nTimeStep/NINTERM, nTimeStep, dt);

		if (parameters::calculate_disk) {
			CommunicateBoundaries(&data[t_data::DENSITY], &data[t_data::V_RADIAL], &data[t_data::V_AZIMUTHAL], &data[t_data::ENERGY]);

			recalculate_derived_disk_quantities(data, true);

			if (!SloppyCFL) {
				local_gas_time_step_cfl = 1.0;
				local_gas_time_step_cfl = condition_cfl(data, data[t_data::V_RADIAL], data[t_data::V_AZIMUTHAL], data[t_data::SOUNDSPEED], DT-dtemp);
				MPI_Allreduce(&local_gas_time_step_cfl, &global_gas_time_step_cfl, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
				dt = (DT-dtemp)/global_gas_time_step_cfl;
			}
			AccreteOntoPlanets(data, dt);

		}
	}

	if (data[t_data::ALPHA_GRAV_MEAN].get_write()) {
		quantities::calculate_alpha_grav_mean_finalize(data, DT);
	}
	if (data[t_data::ALPHA_GRAV_MEAN_1D].get_write()) {
		quantities::calculate_radial_alpha_grav_mean_finalize(data, DT);
	}
	if (data[t_data::ALPHA_REYNOLDS_MEAN].get_write()) {
		quantities::calculate_alpha_reynolds_mean_finalize(data, DT);
	}
	if (data[t_data::ALPHA_REYNOLDS_MEAN_1D].get_write()) {
		quantities::calculate_radial_alpha_reynolds_mean_finalize(data, DT);
	}
}

/**
	In this substep we take into account the source part of Euler equations. We evolve velocities with pressure gradients, gravitational forces and curvature terms
*/
void SubStep1(t_data &data, double dt)
{
	double gradp, gradphi, vt2;
	double invdxtheta;
	double supp_torque=0.0; // for imposed disk drift

	// update v_radial with source terms
	for (unsigned int n_radial = 1; n_radial <= data[t_data::V_RADIAL_SOURCETERMS].get_max_radial()-1; ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::V_RADIAL_SOURCETERMS].get_max_azimuthal(); ++n_azimuthal) {
			// 1/Sigma * dP/dr : Sigma is calculated as a mean value between the neightbour cells
			gradp = 2.0/(data[t_data::DENSITY](n_radial,n_azimuthal)+data[t_data::DENSITY](n_radial-1,n_azimuthal))*(data[t_data::PRESSURE](n_radial,n_azimuthal)-data[t_data::PRESSURE](n_radial-1,n_azimuthal))*InvDiffRmed[n_radial];

			// dPhi/dr
			gradphi = (data[t_data::POTENTIAL](n_radial,n_azimuthal)-data[t_data::POTENTIAL](n_radial-1,n_azimuthal))*InvDiffRmed[n_radial];

			// v_phi^2/r : v_phi^2 is calculated by a mean in both directions
			vt2 = data[t_data::V_AZIMUTHAL](n_radial,n_azimuthal)+data[t_data::V_AZIMUTHAL](n_radial,n_azimuthal == data[t_data::V_AZIMUTHAL].get_max_azimuthal() ? 0 : n_azimuthal+1)+data[t_data::V_AZIMUTHAL](n_radial-1,n_azimuthal)+data[t_data::V_AZIMUTHAL](n_radial-1,n_azimuthal == data[t_data::V_AZIMUTHAL].get_max_azimuthal() ? 0 : n_azimuthal+1);
			vt2 = 0.25*vt2+Rinf[n_radial]*OmegaFrame;
			vt2 = vt2*vt2;

			// add all terms to new v_radial: v_radial_new = v_radial + dt*(source terms)
			data[t_data::V_RADIAL_SOURCETERMS](n_radial,n_azimuthal) = data[t_data::V_RADIAL](n_radial,n_azimuthal)+dt*(-gradp-gradphi+vt2*InvRinf[n_radial]);
		}
	}


	// update v_azimuthal with source terms
	for (unsigned int n_radial = 0; n_radial <= data[t_data::V_AZIMUTHAL_SOURCETERMS].get_max_radial(); ++n_radial) {
		supp_torque = IMPOSEDDISKDRIFT*0.5*pow(Rmed[n_radial],-2.5+SIGMASLOPE);
		invdxtheta = 1.0/(2.0*PI/(double)data[t_data::V_AZIMUTHAL].get_size_azimuthal()*Rmed[n_radial]);

		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::V_AZIMUTHAL_SOURCETERMS].get_max_azimuthal(); ++n_azimuthal) {
			// 1/Sigma 1/r dP/dphi
			gradp = 2.0/(data[t_data::DENSITY](n_radial,n_azimuthal)+data[t_data::DENSITY](n_radial,n_azimuthal == 0 ? data[t_data::DENSITY].get_max_azimuthal() : n_azimuthal-1)) * (data[t_data::PRESSURE](n_radial,n_azimuthal)-data[t_data::PRESSURE](n_radial,n_azimuthal == 0 ? data[t_data::PRESSURE].get_max_azimuthal() : n_azimuthal-1))*invdxtheta;

			// 1/r dPhi/dphi
			gradphi = (data[t_data::POTENTIAL](n_radial,n_azimuthal)-data[t_data::POTENTIAL](n_radial,n_azimuthal == 0 ? data[t_data::POTENTIAL].get_max_azimuthal() : n_azimuthal-1))*invdxtheta;

			// add all terms to new v_azimuthal: v_azimuthal_new = v_azimuthal + dt*(source terms)
			data[t_data::V_AZIMUTHAL_SOURCETERMS](n_radial,n_azimuthal) = data[t_data::V_AZIMUTHAL](n_radial,n_azimuthal)+dt*(-gradp-gradphi);

			// add term for imposed disk drift
			data[t_data::V_AZIMUTHAL_SOURCETERMS](n_radial,n_azimuthal) += dt*supp_torque;
		}
	}

	if (parameters::self_gravity) {
		selfgravity::compute(data[t_data::DENSITY], data[t_data::V_RADIAL_SOURCETERMS], data[t_data::V_AZIMUTHAL_SOURCETERMS], dt, true);
	}

	// compute and add acceleartions due to disc viscosity as a source term
	if (parameters::artificial_viscosity == parameters::artificial_viscosity_TW) {
		viscosity::compute_viscous_terms(data, true);
	} else {
		viscosity::compute_viscous_terms(data, false);
	}

	viscosity::update_velocities_with_viscosity(data, data[t_data::V_RADIAL_SOURCETERMS], data[t_data::V_AZIMUTHAL_SOURCETERMS], dt);


	if ( (parameters::boundary_inner != parameters::boundary_condition_evanescent) || (parameters::boundary_outer != parameters::boundary_condition_evanescent) || (parameters::boundary_inner != parameters::boundary_condition_boundary_layer) || (parameters::boundary_outer != parameters::boundary_condition_boundary_layer) )
		ApplySubKeplerianBoundary(data[t_data::V_AZIMUTHAL_SOURCETERMS]);

	/* if ( !Evanescent )
		ApplySubKeplerianBoundary(&v_azimuthal_sourceterms); */
}

/**
	In this substep we add the articifial viscous pressure source terms.
	Shocks are spread over CVNR zones: von Neumann-Richtmyer viscosity constant; Beware of misprint in Stone and Norman's paper : use C2^2 instead of C2
*/
void SubStep2(t_data &data, double dt)
{
	if (parameters::artificial_viscosity == parameters::artificial_viscosity_SN) {
		double dxtheta, invdxtheta;

		// calculate q_r and q_phi
		for (unsigned int n_radial = 0; n_radial <= data[t_data::Q_R].get_max_radial(); ++n_radial) {
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::Q_R].get_max_azimuthal(); ++n_azimuthal) {
				double dv_r = data[t_data::V_RADIAL_SOURCETERMS](n_radial+1,n_azimuthal)-data[t_data::V_RADIAL_SOURCETERMS](n_radial,n_azimuthal);
				if (dv_r < 0.0) {
					data[t_data::Q_R](n_radial,n_azimuthal) = pow2(parameters::artificial_viscosity_factor)*data[t_data::DENSITY](n_radial,n_azimuthal)*pow2(dv_r);
				} else {
					data[t_data::Q_R](n_radial,n_azimuthal) = 0.0;
				}

				double dv_phi = data[t_data::V_AZIMUTHAL_SOURCETERMS](n_radial,n_azimuthal == data[t_data::V_AZIMUTHAL_SOURCETERMS].get_max_azimuthal() ? 0 : n_azimuthal+1)-data[t_data::V_AZIMUTHAL_SOURCETERMS](n_radial,n_azimuthal);
				if (dv_phi < 0.0) {
					data[t_data::Q_PHI](n_radial,n_azimuthal) = pow2(parameters::artificial_viscosity_factor)*data[t_data::DENSITY](n_radial,n_azimuthal)*pow2(dv_phi);
				} else {
					data[t_data::Q_PHI](n_radial,n_azimuthal) = 0.0;
				}
			}
		}

		// add artificial viscous pressure source term to v_radial
		for (unsigned int n_radial = 1; n_radial <= data[t_data::V_RADIAL_NEW].get_max_radial()-1; ++n_radial) {
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::V_RADIAL_NEW].get_max_azimuthal(); ++n_azimuthal) {
				// 1/Sigma dq_r/dr : Sigma is calculated as a mean value between the neightbour cells
				data[t_data::V_RADIAL_NEW](n_radial,n_azimuthal) = data[t_data::V_RADIAL_SOURCETERMS](n_radial,n_azimuthal)
								- dt*2.0/(data[t_data::DENSITY](n_radial,n_azimuthal)+ data[t_data::DENSITY](n_radial-1,n_azimuthal))*(data[t_data::Q_R](n_radial,n_azimuthal) - data[t_data::Q_R](n_radial-1,n_azimuthal))*InvDiffRmed[n_radial];
			}
		}


		// add artificial viscous pressure source term to v_azimuthal
		for (unsigned int n_radial = 0; n_radial <= data[t_data::V_AZIMUTHAL_NEW].get_max_radial(); ++n_radial) {
			dxtheta = 2.0*PI/(double)data[t_data::DENSITY].get_size_azimuthal()*Rmed[n_radial];
			invdxtheta = 1.0/dxtheta;
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::V_AZIMUTHAL_NEW].get_max_azimuthal(); ++n_azimuthal) {
				// 1/Sigma 1/r dq_phi/dphi : Sigma is calculated as a mean value between the neightbour cells
				data[t_data::V_AZIMUTHAL_NEW](n_radial,n_azimuthal) = data[t_data::V_AZIMUTHAL_SOURCETERMS](n_radial,n_azimuthal)
								- dt*2.0/(data[t_data::DENSITY](n_radial,n_azimuthal) + data[t_data::DENSITY](n_radial,n_azimuthal == 0 ? data[t_data::DENSITY].get_max_azimuthal() : n_azimuthal-1))*(data[t_data::Q_PHI](n_radial,n_azimuthal) - data[t_data::Q_PHI](n_radial,n_azimuthal == 0 ? data[t_data::Q_PHI].get_max_azimuthal() : n_azimuthal-1))*invdxtheta;

			}
		}


		// If gas disk is adiabatic, we add artificial viscosity as a source term for advection of thermal energy polargrid
		if (parameters::Adiabatic) {
			if (parameters::artificial_viscosity_dissipation) {
				for (unsigned int n_radial = 0; n_radial <= data[t_data::ENERGY_INT].get_max_radial(); ++n_radial) {
					dxtheta = 2.0*PI/(double)data[t_data::DENSITY].get_size_azimuthal()*Rmed[n_radial];
					invdxtheta = 1.0/dxtheta;
					for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ENERGY_INT].get_max_azimuthal(); ++n_azimuthal) {
						data[t_data::ENERGY_INT](n_radial,n_azimuthal) = data[t_data::ENERGY](n_radial,n_azimuthal)
						                           - dt*data[t_data::Q_R](n_radial,n_azimuthal)*(data[t_data::V_RADIAL_SOURCETERMS](n_radial+1,n_azimuthal)-data[t_data::V_RADIAL_SOURCETERMS](n_radial,n_azimuthal))*InvDiffRsup[n_radial]
						                           - dt*data[t_data::Q_PHI](n_radial,n_azimuthal)*(data[t_data::V_AZIMUTHAL_SOURCETERMS](n_radial,n_azimuthal == data[t_data::V_AZIMUTHAL_SOURCETERMS].get_max_azimuthal() ? 0 : n_azimuthal+1)-data[t_data::V_AZIMUTHAL_SOURCETERMS](n_radial, n_azimuthal))*invdxtheta;
					}
				}
			} else {
				copy_polargrid(data[t_data::ENERGY_INT],data[t_data::ENERGY]);
			}
		}
	} else {
		copy_polargrid(data[t_data::V_RADIAL_NEW],data[t_data::V_RADIAL_SOURCETERMS]);
		copy_polargrid(data[t_data::V_AZIMUTHAL_NEW],data[t_data::V_AZIMUTHAL_SOURCETERMS]);
		copy_polargrid(data[t_data::ENERGY_INT],data[t_data::ENERGY]);
	}
}

void calculate_qplus(t_data &data) {
	// clear up all Qplus terms
	data[t_data::QPLUS].clear();
	data.qplus_total = 0;

	if (parameters::heating_viscous_enabled) {
		/* We calculate the heating source term Qplus from i=1 to max-1 */
		for (unsigned int n_radial = 1; n_radial <= data[t_data::QPLUS].get_max_radial()-1; ++n_radial) {
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::QPLUS].get_max_azimuthal(); ++n_azimuthal) {
				if (data[t_data::VISCOSITY](n_radial,n_azimuthal) != 0.0) {
					// average tau_r_phi over 4 cells
					double tau_r_phi = 0.25*(
						 data[t_data::TAU_R_PHI](n_radial,n_azimuthal)
						+data[t_data::TAU_R_PHI](n_radial+1,n_azimuthal)
						+data[t_data::TAU_R_PHI](n_radial,n_azimuthal == data[t_data::TAU_R_PHI].get_max_azimuthal() ? 0 : n_azimuthal+1)
						+data[t_data::TAU_R_PHI](n_radial+1,n_azimuthal == data[t_data::TAU_R_PHI].get_max_azimuthal() ? 0 : n_azimuthal+1));

					double qplus = 1.0/(2.0*data[t_data::VISCOSITY](n_radial,n_azimuthal)*data[t_data::DENSITY](n_radial,n_azimuthal))*(pow2(data[t_data::TAU_R_R](n_radial,n_azimuthal))+2*pow2(tau_r_phi)+pow2(data[t_data::TAU_PHI_PHI](n_radial,n_azimuthal)));
					qplus += (2.0/9.0)*data[t_data::VISCOSITY](n_radial,n_azimuthal)*data[t_data::DENSITY](n_radial,n_azimuthal)*pow2(data[t_data::DIV_V](n_radial,n_azimuthal));

					qplus *= parameters::heating_viscous_factor;
					data[t_data::QPLUS](n_radial,n_azimuthal) += qplus;
					sum_without_ghost_cells(data.qplus_total, qplus, n_radial);
				}
			}
		}

		/* We calculate the heating source term Qplus for i=max */
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::QPLUS].get_max_azimuthal(); ++n_azimuthal) {
			if (data[t_data::VISCOSITY](data[t_data::QPLUS].get_max_radial(),n_azimuthal) != 0.0) {
				// power-law extrapolation
				double qplus = data[t_data::QPLUS](data[t_data::QPLUS].get_max_radial()-1, n_azimuthal)*exp( log(data[t_data::QPLUS](data[t_data::QPLUS].get_max_radial()-1, n_azimuthal)/data[t_data::QPLUS](data[t_data::QPLUS].get_max_radial()-2, n_azimuthal)) * log(Rmed[data[t_data::QPLUS].get_max_radial()]/Rmed[data[t_data::QPLUS].get_max_radial()-1]) / log(Rmed[data[t_data::QPLUS].get_max_radial()-1]/Rmed[data[t_data::QPLUS].get_max_radial()-2]) );

				data[t_data::QPLUS](data[t_data::QPLUS].get_max_radial(),n_azimuthal) += qplus;
				sum_without_ghost_cells(data.qplus_total, qplus, data[t_data::QPLUS].get_max_radial());
			}
		}
		/* We calculate the heating source term Qplus for i=0 */
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::QPLUS].get_max_azimuthal(); ++n_azimuthal) {
			if (data[t_data::VISCOSITY](0,n_azimuthal) != 0.0) {
				// power-law extrapolation
				double qplus = data[t_data::QPLUS](1, n_azimuthal)*exp( log(data[t_data::QPLUS](1, n_azimuthal)/data[t_data::QPLUS](2, n_azimuthal)) * log(Rmed[0]/Rmed[1]) / log(Rmed[1]/Rmed[2]) );

				data[t_data::QPLUS](0,n_azimuthal) += qplus;
				sum_without_ghost_cells(data.qplus_total, qplus, 0);
			}
		}
	}

	if (parameters::heating_star_enabled) {
		double ramping = 1.0;
		if (PhysicalTime < parameters::heating_star_ramping_time*DT) {
			ramping = 1.0-pow2(cos(PhysicalTime*PI/2.0/(parameters::heating_star_ramping_time*DT)));
		}

		if (parameters::heating_star_simple) {
			if (!parameters::cooling_radiative_enabled) {
				die("Need to calulate Tau_eff first!\n"); // TODO: make it properly!
			}
			// Simple star heating (see Masterthesis Alexandros Ziampras)
			for (unsigned int n_radial = 1; n_radial <= data[t_data::QPLUS].get_max_radial()-1; ++n_radial) {
				for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::QPLUS].get_max_azimuthal(); ++n_azimuthal) {

					const double distance = Rmed[n_radial];
					const double HoverR = data[t_data::ASPECTRATIO](n_radial, n_azimuthal);
					const double sigma = constants::sigma.get_code_value();
					const double T_star = parameters::star_temperature;
					const double R_star = parameters::star_radius;
					const double tau_eff = data[t_data::TAU_EFF](n_radial, n_azimuthal);
					const double eps = 0.5;  // TODO: add a parameter
					// choose according to Chiang & Goldreich (1997)
					const double dlogH_dlogr = 9.0/7.0;
					// use eq. 7 from Menou & Goodman (2004) (rearranged), Qirr = 2*(1-eps)*L_star/(4 pi r^2)*(dlogH/dlogr - 1) * H/r * 1/Tau_eff
					// here we use (1-eps) = parameters::heating_star_factor
					// L_star = 4 pi R_star^2 sigma_sb T_star^4
					double qplus = 2*(1-eps); // 2*(1-eps)
					qplus *= sigma*pow4(T_star)*pow2(R_star/distance); // *L_star/(4 pi r^2)
					qplus *= dlogH_dlogr - 1; // *(dlogH/dlogr - 1)
					qplus *= HoverR; // * H/r
					qplus /= tau_eff; // * 1/Tau_eff
					data[t_data::QPLUS](n_radial,n_azimuthal) += ramping*qplus;
				}
			}
		} else {
			unsigned int *zbuffer = (unsigned int *)malloc(parameters::zbuffer_size*data[t_data::QPLUS].get_size_azimuthal()*sizeof(unsigned int));

			double dtheta = parameters::zbuffer_maxangle/(double)(parameters::zbuffer_size-1);

			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::QPLUS].get_max_azimuthal(); ++n_azimuthal) {
				// init z-buffer
				for (unsigned int n_theta = 0; n_theta < parameters::zbuffer_size; ++n_theta) {
					zbuffer[n_azimuthal*parameters::zbuffer_size+n_theta] = UINT_MAX;
				}

				for (int n_radial = data[t_data::QPLUS].get_max_radial() - ( CPU_Rank == CPU_Highest ? GHOSTCELLS_B : CPUOVERLAP ); n_radial >= (int)((CPU_Rank == 0) ? GHOSTCELLS_B : CPUOVERLAP); n_radial--) {
					for (unsigned int n_theta = 0; n_theta < parameters::zbuffer_size; ++n_theta) {
						double theta = (double)n_theta*dtheta;

						if (tan(theta) < data[t_data::ASPECTRATIO](n_radial, n_azimuthal)) {
							if (zbuffer[n_azimuthal*parameters::zbuffer_size+n_theta] > IMIN+n_radial) {
								zbuffer[n_azimuthal*parameters::zbuffer_size+n_theta] = IMIN+n_radial;
							} else {
								break;
							}
						}
					}
				}
			}

			// sync
			unsigned int *zbuffer_global = (unsigned int *)malloc(parameters::zbuffer_size*data[t_data::QPLUS].get_size_azimuthal()*sizeof(unsigned int));
			MPI_Allreduce(zbuffer, zbuffer_global, parameters::zbuffer_size*data[t_data::QPLUS].get_size_azimuthal(), MPI_UNSIGNED, MPI_MIN, MPI_COMM_WORLD);
			free(zbuffer);
			zbuffer = zbuffer_global;

			// calculate visiblity
			for (unsigned int n_radial = 1; n_radial <= data[t_data::QPLUS].get_max_radial()-1; ++n_radial) {
				for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::QPLUS].get_max_azimuthal(); ++n_azimuthal) {
					// check for self-shadowing
					unsigned int n_theta = 0;

					// get next nt
					while ((tan(n_theta*dtheta) <= data[t_data::ASPECTRATIO](n_radial, n_azimuthal)) && (n_theta < parameters::zbuffer_size))
						n_theta++;

					if (zbuffer[n_azimuthal*parameters::zbuffer_size+n_theta] >= IMIN+n_radial) {
						data[t_data::VISIBILITY](n_radial, n_azimuthal) = 1.0;
					} else {
						data[t_data::VISIBILITY](n_radial, n_azimuthal) = 0.0;
					}
				}
			}

			for (unsigned int n_radial = 1; n_radial <= data[t_data::QPLUS].get_max_radial()-1; ++n_radial) {
				for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::QPLUS].get_max_azimuthal(); ++n_azimuthal) {
					// calculate "mean" visibility
					double visibility = 0;

					visibility += data[t_data::VISIBILITY](n_radial-1, n_azimuthal == 0 ? data[t_data::VISIBILITY].get_max_azimuthal() : n_azimuthal-1);
					visibility += data[t_data::VISIBILITY](n_radial-1, n_azimuthal);
					visibility += data[t_data::VISIBILITY](n_radial-1, n_azimuthal == data[t_data::VISIBILITY].get_max_azimuthal() ? 0 : n_azimuthal+1);

					visibility += data[t_data::VISIBILITY](n_radial, n_azimuthal == 0 ? data[t_data::VISIBILITY].get_max_azimuthal() : n_azimuthal-1);
					visibility += data[t_data::VISIBILITY](n_radial, n_azimuthal);
					visibility += data[t_data::VISIBILITY](n_radial, n_azimuthal == data[t_data::VISIBILITY].get_max_azimuthal() ? 0 : n_azimuthal+1);

					visibility += data[t_data::VISIBILITY](n_radial+1, n_azimuthal == 0 ? data[t_data::VISIBILITY].get_max_azimuthal() : n_azimuthal-1);
					visibility += data[t_data::VISIBILITY](n_radial+1, n_azimuthal);
					visibility += data[t_data::VISIBILITY](n_radial+1, n_azimuthal == data[t_data::VISIBILITY].get_max_azimuthal() ? 0 : n_azimuthal+1);

					visibility /= 9.0;

					// see Günther et. al (2004) or Phil Armitage "Astrophysics of Planet Format" p. 46
					double alpha = (data[t_data::ASPECTRATIO](n_radial, n_azimuthal)*Rmed[n_radial]-data[t_data::ASPECTRATIO](n_radial-1, n_azimuthal)*Rmed[n_radial-1])*InvDiffRmed[n_radial] - data[t_data::ASPECTRATIO](n_radial, n_azimuthal);

					if (alpha < 0.0) {
						alpha = 0.0;
					}

					// primary star
					double qplus = ramping*visibility*parameters::heating_star_factor*2.0*alpha*constants::sigma.get_code_value()*pow4(parameters::star_temperature)*pow2(parameters::star_radius/Rmed[n_radial]);

					data[t_data::QPLUS](n_radial,n_azimuthal) += qplus;
					sum_without_ghost_cells(data.qplus_total, qplus, n_radial);

					/* // secondary star/plantes
					for (unsigned int planet = 1; planet > data.get_planetary_system().get_number_of_planets(); ++planet) {
						double planet_x = data.get_planetary_system().get_planet(planet).get_x();
						double planet_y = data.get_planetary_system().get_planet(planet).get_y();

						double cell_x = Rmed[n_radial]*cos((double)n_azimuthal/(double)data[t_data::QPLUS].get_size_azimuthal()*2.0*PI);
						double cell_y = Rmed[n_radial]*sin((double)n_azimuthal/(double)data[t_data::QPLUS].get_size_azimuthal()*2.0*PI);

						double distance = sqrt(pow2(planet_x-cell_x)+pow2(planet_y-cell_y));

						data[t_data::QPLUS](n_radial,n_azimuthal) += parameters::heating_star_factor*alpha*constants::sigma.get_code_value()*pow4(data.get_planetary_system().get_planet(planet).get_temperature())*pow2(data.get_planetary_system().get_planet(planet).get_radius()/distance);
					}
					*/
				}
			}

			free(zbuffer);
		}
	}
}

void calculate_qminus(t_data &data) {
	// clear up all Qminus terms
	data[t_data::QMINUS].clear();
	data.qminus_total = 0;

	// beta cooling
	if (parameters::cooling_beta_enabled) {
		for (unsigned int n_radial = 0; n_radial <= data[t_data::QPLUS].get_max_radial(); ++n_radial) {
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::QPLUS].get_max_azimuthal(); ++n_azimuthal) {
				// Q- = E Omega/beta
				double qminus = data[t_data::ENERGY](n_radial, n_azimuthal)*(0.5*(data[t_data::V_AZIMUTHAL](n_radial,n_azimuthal)+data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal == data[t_data::V_AZIMUTHAL].get_max_azimuthal() ? 0 : n_azimuthal+1))/Rmed[n_radial])/parameters::cooling_beta;

				data[t_data::QMINUS](n_radial, n_azimuthal) += qminus;
				sum_without_ghost_cells(data.qminus_total, qminus, n_radial);
			}
		}
	}

	// local radiative cooling
	if (parameters::cooling_radiative_enabled) {
		double kappaCGS;
		double temperatureCGS;
		double densityCGS;

		for (unsigned int n_radial = 0; n_radial <= data[t_data::QMINUS].get_max_radial(); ++n_radial) {
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::QMINUS].get_max_azimuthal(); ++n_azimuthal) {
				// calculate Rosseland mean opacity kappa. opaclin needs values in cgs units
				temperatureCGS = data[t_data::TEMPERATURE](n_radial, n_azimuthal)*units::temperature;

				// TODO: user aspect ratio
				densityCGS = data[t_data::DENSITY](n_radial, n_azimuthal)/(parameters::density_factor*data[t_data::SOUNDSPEED](n_radial, n_azimuthal)/sqrt(ADIABATICINDEX)/omega_kepler(Rmed[n_radial]) )*units::density;

				kappaCGS = opacity::opacity(densityCGS,temperatureCGS);

				data[t_data::KAPPA](n_radial, n_azimuthal) = parameters::kappa_factor*kappaCGS*units::opacity.get_inverse_cgs_factor();

				// mean vertical optical depth: tau = 1/2 kappa Sigma
				data[t_data::TAU](n_radial, n_azimuthal) = parameters::tau_factor*(1.0/parameters::density_factor)*data[t_data::KAPPA](n_radial, n_azimuthal)*data[t_data::DENSITY](n_radial, n_azimuthal);

				//  tau_eff = 3/8 tau + sqrt(3)/4 + 1/(4*tau+tau_min)
				data[t_data::TAU_EFF](n_radial, n_azimuthal) = 3.0/8.0*data[t_data::TAU](n_radial, n_azimuthal) + sqrt(3.0)/4.0 + 1.0/(4.0*data[t_data::TAU](n_radial, n_azimuthal)+0.01);

				// effective temperature: T_eff^4 tau_eff = T^4
				//temperatureEff = (*Temperature)(n_radial, n_azimuthal)/(pow(tauEff,1.0/4.0));

				// Q = 2 sigma_R T_eff^4
				//(*Qminus)(n_radial, n_azimuthal) = 2*(constants::sigma.get_code_value())*pow(temperatureEff,4);
				double qminus = parameters::cooling_radiative_factor*2*(constants::sigma.get_code_value())*pow4(data[t_data::TEMPERATURE](n_radial, n_azimuthal))/data[t_data::TAU_EFF](n_radial, n_azimuthal);

				data[t_data::QMINUS](n_radial, n_azimuthal) += qminus;
				sum_without_ghost_cells(data.qminus_total, qminus, n_radial);
			}
		}
	}
}

/**
	In this substep we take into account the source part of energy equation. We evolve internal energy with compression/dilatation and heating terms
*/
void SubStep3(t_data &data, double dt)
{
	double num, den;

	calculate_qminus(data); // first to calculate teff
    calculate_qplus(data);


	// calculate tau_cool if needed for output
	if (data[t_data::TAU_COOL].get_write_1D() || data[t_data::TAU_COOL].get_write_2D()) {
		for (unsigned int n_radial = 0; n_radial <= data[t_data::TAU_COOL].get_max_radial(); ++n_radial) {
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::TAU_COOL].get_max_azimuthal(); ++n_azimuthal) {
					data[t_data::TAU_COOL](n_radial, n_azimuthal) = data[t_data::ENERGY](n_radial, n_azimuthal)/data[t_data::QMINUS](n_radial, n_azimuthal);
			}
		}
	}

	// calculate pDV for write out
	if (data[t_data::P_DIVV].get_write_1D() || data[t_data::P_DIVV].get_write_2D() || parameters::radiative_diffusion_enabled) {
		data.pdivv_total = 0;
		for (unsigned int n_radial = 0; n_radial <= data[t_data::ENERGY_NEW].get_max_radial(); ++n_radial) {
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ENERGY_NEW].get_max_azimuthal(); ++n_azimuthal) {
				double pdivv = (ADIABATICINDEX-1.0)*dt*data[t_data::DIV_V](n_radial, n_azimuthal)*data[t_data::ENERGY](n_radial, n_azimuthal);
				data[t_data::P_DIVV](n_radial, n_azimuthal) = pdivv;

				sum_without_ghost_cells(data.pdivv_total, pdivv, n_radial);
			}
		}
	}

	// Now we can update energy with source terms
	for (unsigned int n_radial = 0; n_radial <= data[t_data::ENERGY_NEW].get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ENERGY_NEW].get_max_azimuthal(); ++n_azimuthal) {
			// original
			/*num = dt*data[t_data::QPLUS](n_radial, n_azimuthal)
			    - dt*data[t_data::QMINUS](n_radial, n_azimuthal)
			    + data[t_data::ENERGY_INT](n_radial, n_azimuthal);
			den = 1.0+(ADIABATICINDEX-1.0)*dt*data[t_data::DIV_V](n_radial, n_azimuthal);
			*/
			// ZEUS2D like
			/*num = (1.0-(dt/2.0)*(ADIABATICINDEX-1.0)*data[t_data::DIV_V](n_radial, n_azimuthal))*data[t_data::ENERGY_INT](n_radial, n_azimuthal)+dt*data[t_data::QPLUS](n_radial, n_azimuthal)-dt*data[t_data::QMINUS](n_radial, n_azimuthal);
			den = (1.0+(dt/2.0)*(ADIABATICINDEX-1.0)*data[t_data::DIV_V](n_radial, n_azimuthal)); 
			*/

			// TWAM original with H/R = initial:
			// const double alpha = 1.0 + 2.0 * ASPECTRATIO*pow(Rmed[n_radial],1+FLARINGINDEX)*4.0
			// 	*constants::sigma/constants::c
			// 	/pow4((constants::R/parameters::MU)/(ADIABATICINDEX-1.0)*data[t_data::DENSITY](n_radial, n_azimuthal))
			// 	*pow3(data[t_data::ENERGY_INT](n_radial, n_azimuthal));

			const double R = Rmed[n_radial];
		    const double H = data[t_data::ASPECTRATIO](n_radial, n_azimuthal) * R ;
			const double inv_pow4 = pow4( (parameters::MU * (ADIABATICINDEX-1.0)) / (constants::R * data[t_data::DENSITY](n_radial, n_azimuthal)) );
			const double pow3_eint = pow3(data[t_data::ENERGY_INT](n_radial, n_azimuthal));
			double alpha = 1.0 + 2.0 * H *4.0*constants::sigma/constants::c * inv_pow4 * pow3_eint;

			num = dt*data[t_data::QPLUS](n_radial, n_azimuthal)
			    - dt*data[t_data::QMINUS](n_radial, n_azimuthal)
			    + alpha*data[t_data::ENERGY_INT](n_radial, n_azimuthal);
			den = alpha+(ADIABATICINDEX-1.0)*dt*data[t_data::DIV_V](n_radial, n_azimuthal);

			data[t_data::ENERGY_NEW](n_radial, n_azimuthal) = num/den;
		}
	}
}

static inline double flux_limiter(double R) {
	// flux limiter
	if (R <= 2) {
		return 2.0/(3+sqrt(9+10*pow2(R)));
	} else {
		return 10.0/(10*R+9+sqrt(180*R+81));
	}
}

void radiative_diffusion(t_data &data, double dt) {
	static bool grids_allocated = false;
	static t_polargrid Ka, Kb;
	static t_polargrid A, B, C, D, E;
	static t_polargrid Told;
	static double *SendInnerBoundary, *SendOuterBoundary, *RecvInnerBoundary, *RecvOuterBoundary;

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
		SendInnerBoundary = (double*)malloc(NAzimuthal*CPUOVERLAP*sizeof(double));
		SendOuterBoundary = (double*)malloc(NAzimuthal*CPUOVERLAP*sizeof(double));
		RecvInnerBoundary = (double*)malloc(NAzimuthal*CPUOVERLAP*sizeof(double));
		RecvOuterBoundary = (double*)malloc(NAzimuthal*CPUOVERLAP*sizeof(double));

		grids_allocated = true;
	}

	// update temperature, soundspeed and aspect ratio
	compute_temperature(data, true);
	compute_sound_speed(data, true);
	compute_aspect_ratio(data, true);

	double dphi = 2.0*PI/(double)(data[t_data::TEMPERATURE].get_size_azimuthal());

	// calcuate Ka for K(i/2,j)
	for (unsigned int n_radial = 1; n_radial <= Ka.get_max_radial()-1; ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= Ka.get_max_azimuthal(); ++n_azimuthal) {
			unsigned int n_azimuthal_plus = (n_azimuthal == Ka.get_max_azimuthal() ? 0 : n_azimuthal + 1);
			unsigned int n_azimuthal_minus = (n_azimuthal == 0 ? Ka.get_max_azimuthal() : n_azimuthal - 1);

			// average temperature radially
			double temperature = 0.5*(data[t_data::TEMPERATURE](n_radial-1, n_azimuthal)+data[t_data::TEMPERATURE](n_radial, n_azimuthal));
			double density     = 0.5*(data[t_data::DENSITY](n_radial-1, n_azimuthal)    +data[t_data::DENSITY](n_radial, n_azimuthal));
			double aspectratio = 0.5*(data[t_data::ASPECTRATIO](n_radial-1, n_azimuthal)+data[t_data::ASPECTRATIO](n_radial, n_azimuthal));

			double temperatureCGS = temperature*units::temperature;
			double H = aspectratio*Ra[n_radial];
			double densityCGS = density/(parameters::density_factor*H)*units::density;

			double kappaCGS = opacity::opacity(densityCGS,temperatureCGS);
			double kappa = parameters::kappa_factor*kappaCGS*units::opacity.get_inverse_cgs_factor();

			double denom = 1.0/(density*kappa);

			// Levermore & Pomraning 1981
			// R = 4 |nabla T\/T * 1/(rho kappa)
			double dT_dr = (data[t_data::TEMPERATURE](n_radial, n_azimuthal)-data[t_data::TEMPERATURE](n_radial-1, n_azimuthal))*InvDiffRmed[n_radial];
			double dT_dphi = InvRinf[n_radial]*(
			 0.5*(data[t_data::TEMPERATURE](n_radial-1, n_azimuthal_plus)+data[t_data::TEMPERATURE](n_radial, n_azimuthal_plus))
			-0.5*(data[t_data::TEMPERATURE](n_radial-1, n_azimuthal_minus)+data[t_data::TEMPERATURE](n_radial, n_azimuthal_minus))
			)/(2*dphi);

			double nabla_T = sqrt(pow2(dT_dr) + pow2(dT_dphi));

			double R;

			if ((n_radial > 1) && (n_radial < Ka.get_max_radial()-1)) {
				R = 4.0*nabla_T/temperature*denom*H*parameters::density_factor;
			} else {
				unsigned int n_radial_adjusted;

				if (n_radial == 1) {
					n_radial_adjusted = n_radial+1;
				} else {
					n_radial_adjusted = n_radial-1;
				}
				double temperature = 0.5*(data[t_data::TEMPERATURE](n_radial_adjusted-1, n_azimuthal)+data[t_data::TEMPERATURE](n_radial_adjusted, n_azimuthal));
				double density     = 0.5*(data[t_data::DENSITY](n_radial_adjusted-1, n_azimuthal)    +data[t_data::DENSITY](n_radial_adjusted, n_azimuthal));
				double aspectratio = 0.5*(data[t_data::ASPECTRATIO](n_radial_adjusted-1, n_azimuthal)+data[t_data::ASPECTRATIO](n_radial_adjusted, n_azimuthal));

				double temperatureCGS = temperature*units::temperature;
				double H = aspectratio*Ra[n_radial_adjusted];
				double densityCGS = density/(parameters::density_factor*H)*units::density;

				double kappaCGS = opacity::opacity(densityCGS,temperatureCGS);
				double kappa = parameters::kappa_factor*kappaCGS*units::opacity.get_inverse_cgs_factor();

				double denom = 1.0/(density*kappa);

				// Levermore & Pomraning 1981
				// R = 4 |nabla T\/T * 1/(rho kappa)
				double dT_dr = (data[t_data::TEMPERATURE](n_radial_adjusted, n_azimuthal)-data[t_data::TEMPERATURE](n_radial_adjusted-1, n_azimuthal))*InvDiffRmed[n_radial];
				double dT_dphi = InvRinf[n_radial]*(
				 0.5*(data[t_data::TEMPERATURE](n_radial_adjusted-1, n_azimuthal_plus)+data[t_data::TEMPERATURE](n_radial_adjusted, n_azimuthal_plus))
				-0.5*(data[t_data::TEMPERATURE](n_radial_adjusted-1, n_azimuthal_minus)+data[t_data::TEMPERATURE](n_radial_adjusted, n_azimuthal_minus))
				)/(2*dphi);

				double nabla_T = sqrt(pow2(dT_dr) + pow2(dT_dphi));

				R = 4.0*nabla_T/temperature*denom*H*parameters::density_factor;			
			}

			double lambda = flux_limiter(R);

			Ka(n_radial, n_azimuthal) = 8.0*4.0*constants::sigma.get_code_value()*lambda*H*pow3(temperature)*denom;
			//Ka(n_radial, n_azimuthal) = 16.0*parameters::density_factor*constants::sigma.get_code_value()*lambda*H*pow3(temperature)*denom;
		}
	}

	// calcuate Kb for K(i,j/2)
	for (unsigned int n_radial = 1; n_radial <= Kb.get_max_radial()-1; ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= Kb.get_max_azimuthal(); ++n_azimuthal) {
			//unsigned int n_azimuthal_plus = (n_azimuthal == Kb.get_max_azimuthal() ? 0 : n_azimuthal + 1);
			unsigned int n_azimuthal_minus = (n_azimuthal == 0 ? Kb.get_max_azimuthal() : n_azimuthal - 1);

			// average temperature azimuthally
			double temperature = 0.5*(data[t_data::TEMPERATURE](n_radial, n_azimuthal_minus)+data[t_data::TEMPERATURE](n_radial, n_azimuthal));
			double density     = 0.5*(data[t_data::DENSITY](n_radial, n_azimuthal_minus)    +data[t_data::DENSITY](n_radial, n_azimuthal));
			double aspectratio = 0.5*(data[t_data::ASPECTRATIO](n_radial, n_azimuthal_minus)+data[t_data::ASPECTRATIO](n_radial, n_azimuthal));

			double temperatureCGS = temperature*units::temperature;
			double H = aspectratio*Rb[n_radial];
			double densityCGS = density/(parameters::density_factor*H)*units::density;

			double kappaCGS = opacity::opacity(densityCGS,temperatureCGS);
			double kappa = parameters::kappa_factor*kappaCGS*units::opacity.get_inverse_cgs_factor();

			double denom = 1.0/(density*kappa);

			// Levermore & Pomraning 1981
			// R = 4 |nabla T\/T * 1/(rho kappa)
			double dT_dr = (
			 0.5*(data[t_data::TEMPERATURE](n_radial-1, n_azimuthal_minus)+data[t_data::TEMPERATURE](n_radial-1, n_azimuthal))
			-0.5*(data[t_data::TEMPERATURE](n_radial+1, n_azimuthal_minus)+data[t_data::TEMPERATURE](n_radial+1, n_azimuthal))
			)/(Ra[n_radial-1]-Ra[n_radial+1]);
			double dT_dphi = InvRmed[n_radial]*(data[t_data::TEMPERATURE](n_radial, n_azimuthal)-data[t_data::TEMPERATURE](n_radial, n_azimuthal_minus))/dphi;

			double nabla_T = sqrt(pow2(dT_dr) + pow2(dT_dphi));

			double R = 4.0*nabla_T/temperature*denom*H*parameters::density_factor;

			double lambda = flux_limiter(R);
			/*if (n_radial == 4) {
				printf("kb: phi=%lg\tR=%lg\tlambda=%lg\tdphi=%lg\tdr=%lg\tnabla=%lg\tT=%lg\tH=%lg\n", dphi*n_azimuthal, R, lambda,dT_dphi,dT_dr,nabla_T,temperature,H);
			}*/

			Kb(n_radial, n_azimuthal) = 8*4*constants::sigma.get_code_value()*lambda*H*pow3(temperature)*denom;
			//Kb(n_radial, n_azimuthal) = 16.0*parameters::density_factor*constants::sigma.get_code_value()*lambda*H*pow3(temperature)*denom;
		}
	}

	double c_v  = constants::R/(parameters::MU*(ADIABATICINDEX-1.0));

	// calculate A,B,C,D,E
	for (unsigned int n_radial = 1; n_radial <= data[t_data::TEMPERATURE].get_max_radial()-1; ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::TEMPERATURE].get_max_azimuthal(); ++n_azimuthal) {
			double H = data[t_data::ASPECTRATIO](n_radial, n_azimuthal)*Rb[n_radial];
			// -dt H /(Sigma * c_v)
			double common_factor = -dt*parameters::density_factor*H/(data[t_data::DENSITY](n_radial, n_azimuthal)*c_v);

			// 2/(dR^2)
			double common_AC= common_factor*2.0/(pow2(Ra[n_radial+1])-pow2(Ra[n_radial]));
			A(n_radial, n_azimuthal) = common_AC*Ka(n_radial, n_azimuthal)*Ra[n_radial]*InvDiffRmed[n_radial];
			C(n_radial, n_azimuthal) = common_AC*Ka(n_radial+1, n_azimuthal)*Ra[n_radial+1]*InvDiffRmed[n_radial+1];

			// 1/(r^2 dphi^2)
			double common_DE = common_factor/(pow2(Rb[n_radial])*pow2(dphi));
			D(n_radial, n_azimuthal) = common_DE*Kb(n_radial, n_azimuthal);
			E(n_radial, n_azimuthal) = common_DE*Kb(n_radial, n_azimuthal == Kb.get_max_azimuthal() ? 0 : n_azimuthal+1);

			B(n_radial, n_azimuthal) = - A(n_radial, n_azimuthal) - C(n_radial, n_azimuthal) - D(n_radial, n_azimuthal) - E(n_radial, n_azimuthal) + 1.0;

			Told(n_radial, n_azimuthal) = data[t_data::TEMPERATURE](n_radial, n_azimuthal);

			/*double energy_change = dt*data[t_data::QPLUS](n_radial, n_azimuthal)
			    - dt*data[t_data::QMINUS](n_radial, n_azimuthal)
			    - dt*data[t_data::P_DIVV](n_radial, n_azimuthal);

			double temperature_change = MU/R*(ADIABATICINDEX-1.0)*energy_change/data[t_data::DENSITY](n_radial,n_azimuthal);
			Told(n_radial, n_azimuthal) += temperature_change;

			if (Told(n_radial, n_azimuthal) < parameters::minimum_temperature*units::temperature.get_inverse_cgs_factor()) {
				data[t_data::TEMPERATURE](n_radial, n_azimuthal) = parameters::minimum_temperature*units::temperature.get_inverse_cgs_factor();
			}
			*/
		}
	}

	static unsigned int old_iterations = parameters::radiative_diffusion_max_iterations;
	static int direction = 1;
	static double omega = parameters::radiative_diffusion_omega;

	unsigned int iterations = 0;
	double absolute_norm = DBL_MAX;
	double norm_change = DBL_MAX;

	int l = CPUOVERLAP*NAzimuthal;
	int oo = (data[t_data::TEMPERATURE].Nrad-CPUOVERLAP)*NAzimuthal;
	int o = (data[t_data::TEMPERATURE].Nrad-2*CPUOVERLAP)*NAzimuthal;

	// do SOR
	while ((norm_change > 1e-12) && (parameters::radiative_diffusion_max_iterations > iterations)) {
		// if ((CPU_Rank == CPU_Highest) && parameters::boundary_outer == parameters::boundary_condition_open) {
		// 	// set temperature to T_min in outermost cells
		// 	for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::TEMPERATURE].get_max_azimuthal(); ++n_azimuthal) {
		// 		data[t_data::TEMPERATURE](data[t_data::TEMPERATURE].get_max_radial(), n_azimuthal) = parameters::minimum_temperature*units::temperature.get_inverse_cgs_factor();
		// 	}
		// }

		// if ((CPU_Rank == 0) && parameters::boundary_inner == parameters::boundary_condition_open) {
		// 	// set temperature to T_min in innermost cells
		// 	for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::TEMPERATURE].get_max_azimuthal(); ++n_azimuthal) {
		// 		data[t_data::TEMPERATURE](0, n_azimuthal) = parameters::minimum_temperature*units::temperature.get_inverse_cgs_factor();
		// 	}
		// }
                boundary_conditions::apply_boundary_condition(data, dt, false);


		norm_change = absolute_norm;
		absolute_norm = 0.0;

		for (unsigned int n_radial = 1; n_radial <= data[t_data::TEMPERATURE].get_max_radial()-1; ++n_radial) {
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::TEMPERATURE].get_max_azimuthal(); ++n_azimuthal) {
				double old_value = data[t_data::TEMPERATURE](n_radial, n_azimuthal);
				unsigned int n_azimuthal_plus = (n_azimuthal == data[t_data::TEMPERATURE].get_max_azimuthal() ? 0 : n_azimuthal + 1);
				unsigned int n_azimuthal_minus = (n_azimuthal == 0 ? data[t_data::TEMPERATURE].get_max_azimuthal() : n_azimuthal - 1);

				data[t_data::TEMPERATURE](n_radial, n_azimuthal) = (1.0 - omega) * data[t_data::TEMPERATURE](n_radial, n_azimuthal) - omega/B(n_radial, n_azimuthal) * (
					  A(n_radial, n_azimuthal) * data[t_data::TEMPERATURE](n_radial-1, n_azimuthal)
					+ C(n_radial, n_azimuthal) * data[t_data::TEMPERATURE](n_radial+1, n_azimuthal)
					+ D(n_radial, n_azimuthal) * data[t_data::TEMPERATURE](n_radial, n_azimuthal_minus)
					+ E(n_radial, n_azimuthal) * data[t_data::TEMPERATURE](n_radial, n_azimuthal_plus)
					- Told(n_radial, n_azimuthal));

        // only non ghostcells to norm and don't count overlap cell's twice
        bool isnot_ghostcell_rank_0 = n_radial > ((CPU_Rank == 0) ? GHOSTCELLS_B : CPUOVERLAP);
        bool isnot_ghostcell_rank_highest = (n_radial < (data[t_data::TEMPERATURE].get_max_radial()
                                                         - ((CPU_Rank == CPU_Highest) ? GHOSTCELLS_B : CPUOVERLAP)));
        
        if (isnot_ghostcell_rank_0 && isnot_ghostcell_rank_highest) {
					absolute_norm += pow2(old_value - data[t_data::TEMPERATURE](n_radial, n_azimuthal));
				}
			}
		}

		double tmp = absolute_norm;
		MPI_Allreduce(&tmp, &absolute_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		absolute_norm = sqrt(absolute_norm)/(GlobalNRadial * NAzimuthal);

		norm_change = fabs(absolute_norm - norm_change);
		iterations++;

		// communicate with other nodes
		memcpy(SendInnerBoundary, data[t_data::TEMPERATURE].Field+l, l*sizeof(double));
		memcpy(SendOuterBoundary, data[t_data::TEMPERATURE].Field+o, l*sizeof(double));

		MPI_Request req1, req2, req3, req4;

		if (CPU_Rank%2 == 0) {
			if (CPU_Rank != 0) {
				MPI_Isend(SendInnerBoundary, NAzimuthal * CPUOVERLAP, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
				MPI_Irecv(RecvInnerBoundary, NAzimuthal * CPUOVERLAP, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req2);
			}
			if (CPU_Rank != CPU_Highest) {
				MPI_Isend(SendOuterBoundary, NAzimuthal * CPUOVERLAP, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req3);
				MPI_Irecv(RecvOuterBoundary, NAzimuthal * CPUOVERLAP, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req4);
			}
		} else {
			if (CPU_Rank != CPU_Highest) {
				MPI_Irecv(RecvOuterBoundary, NAzimuthal * CPUOVERLAP, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req3);
				MPI_Isend(SendOuterBoundary, NAzimuthal * CPUOVERLAP, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req4);
			}
			if (CPU_Rank != 0) {
				MPI_Irecv(RecvInnerBoundary, NAzimuthal * CPUOVERLAP, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
				MPI_Isend(SendInnerBoundary, NAzimuthal * CPUOVERLAP, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req2);
			}
		}

		if (CPU_Rank != 0) {
            MPI_Wait(&req1, &global_MPI_Status);
            MPI_Wait(&req2, &global_MPI_Status);
			memcpy(data[t_data::TEMPERATURE].Field, RecvInnerBoundary, l*sizeof(double));
		}

		if (CPU_Rank != CPU_Highest) {
            MPI_Wait(&req3, &global_MPI_Status);
            MPI_Wait(&req4, &global_MPI_Status);
			memcpy(data[t_data::TEMPERATURE].Field+oo, RecvOuterBoundary, l*sizeof(double));
		}
	}

	if (iterations == parameters::radiative_diffusion_max_iterations) {
		logging::print_master(LOG_WARNING "Maximum iterations (%u) reached in radiative_diffusion (omega = %lg). Norm is %lg with a last change of %lg.\n", parameters::radiative_diffusion_max_iterations, omega, absolute_norm, norm_change);
	}

	// adapt omega
	if (old_iterations < iterations) {
		direction *= -1;
	}

	if (parameters::radiative_diffusion_omega_auto_enabled) {
		omega += direction*0.01;
	}

	if (omega >= 2.0) {
		omega = 2.0;
		direction = -1;
	}

	if (omega <= 1.0) {
		omega = 1.0;
		direction = 1;
	}

	old_iterations = iterations;

	logging::print_master(LOG_VERBOSE "%u iterations, omega=%lf\n", iterations, omega);

	// compute energy from temperature
	for (unsigned int n_radial = 0; n_radial <= data[t_data::ENERGY].get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal(); ++n_azimuthal) {
			data[t_data::ENERGY](n_radial,n_azimuthal) = data[t_data::TEMPERATURE](n_radial,n_azimuthal)*data[t_data::DENSITY](n_radial,n_azimuthal)/(ADIABATICINDEX-1.0)/parameters::MU*constants::R;
		}
	}
}

/**
	\param VRadial radial velocity polar grid
	\param VAzimuthal azimuthal velocity polar grid
	\param SoundSpeed sound speed polar grid
	\param deltaT
*/
double condition_cfl(t_data& data, t_polargrid &v_radial, t_polargrid &v_azimuthal, t_polargrid &soundspeed, double deltaT)
{
	double v_mean[v_radial.get_size_radial()], v_residual[v_radial.get_size_azimuthal()];
	double dtGlobal, dtLocal;

	// debugging variables
	double viscRadial = 0.0, viscAzimuthal = 0.0;
	unsigned int n_azimuthal_debug = 0, n_radial_debug = 0;
	double itdbg1 = DBL_MAX, itdbg2 = DBL_MAX, itdbg3 = DBL_MAX, itdbg4 = DBL_MAX, itdbg5 = DBL_MAX, mdtdbg = DBL_MAX;

	dtGlobal = DBL_MAX;

	// Calculate and fill VMean array
	for (unsigned int n_radial = 0; n_radial <= v_azimuthal.get_max_radial(); ++n_radial) {
		v_mean[n_radial] = 0.0;
		for (unsigned int n_azimuthal = 0; n_azimuthal <= v_azimuthal.get_max_azimuthal(); ++n_azimuthal) {
			v_mean[n_radial] += v_azimuthal(n_radial, n_azimuthal);
		}
		v_mean[n_radial] /= (double)(v_azimuthal.get_size_azimuthal());
	}

	for (unsigned int n_radial = One_or_active; n_radial < Max_or_active; ++n_radial) {
		// cell sizes in radial & azimuthal direction
		double dxRadial = Rsup[n_radial]-Rinf[n_radial];
		double dxAzimuthal = Rmed[n_radial]*2.0*PI/(double)(v_radial.get_size_azimuthal());

		for (unsigned int n_azimuthal = 0; n_azimuthal <= v_radial.get_max_azimuthal(); ++n_azimuthal) {
			if (FastTransport) {
				// FARGO algorithm
				v_residual[n_azimuthal] = v_azimuthal(n_radial,n_azimuthal)-v_mean[n_radial];
			} else {
				// Standard algorithm
				v_residual[n_azimuthal] = v_azimuthal(n_radial,n_azimuthal);
			}
		}

		// there is no v_residual[v_radial.Nsec]
		//v_residual[v_radial.Nsec]=v_residual[0];

		for (unsigned int n_azimuthal = 0; n_azimuthal <= v_radial.get_max_azimuthal(); ++n_azimuthal) {
			double invdt1, invdt2, invdt3, invdt4, invdt5;

			// velocity differences in radial & azimuthal direction
			double dvRadial = v_radial(n_radial+1, n_azimuthal)-v_radial(n_radial, n_azimuthal);
			double dvAzimuthal = v_azimuthal(n_radial, n_azimuthal == v_radial.get_max_azimuthal() ? 0 : n_azimuthal+1)-v_azimuthal(n_radial, n_azimuthal);

			// sound speed limit
			invdt1 = soundspeed(n_radial,n_azimuthal)/(min(dxRadial,dxAzimuthal));

			// radial motion limit
			invdt2 = fabs(v_radial(n_radial,n_azimuthal))/dxRadial;

			// residual circular motion limit
			invdt3 = fabs(v_residual[n_azimuthal])/dxAzimuthal;

			// artificial viscosity limit
			if (parameters::artificial_viscosity == parameters::artificial_viscosity_SN) {
				// TODO: Change to sizes defined by constants of compiler
				if (dvRadial >= 0.0) {
					dvRadial = 1e-10;
				} else {
					dvRadial = -dvRadial;
				}

				if (dvAzimuthal >= 0.0) {
					dvAzimuthal = 1e-10;
				} else {
					dvAzimuthal = -dvAzimuthal;
				}

				invdt4 = 4.0*pow2(parameters::artificial_viscosity_factor)*max(dvRadial/dxRadial,dvAzimuthal/dxAzimuthal);
			} else {
				invdt4 = 0.0;
			}

			// kinematic viscosity limit
			// TODO: Factor 4 on errors!
			invdt5 = 4.0*data[t_data::VISCOSITY](n_radial, n_azimuthal)*max(1/pow2(dxRadial),1/pow2(dxAzimuthal));

			// calculate new dt based on different limits
			dtLocal = parameters::CFL/sqrt(pow2(invdt1)+pow2(invdt2)+pow2(invdt3)+pow2(invdt4)+pow2(invdt5));

			if (dtLocal < dtGlobal) {
				dtGlobal = dtLocal;
				if (debug) {
					n_radial_debug = n_radial;
					n_azimuthal_debug = n_azimuthal;
					if (invdt1 != 0)
						itdbg1=1.0/invdt1;
					if (invdt2 != 0)
						itdbg2=1.0/invdt2;
					if (invdt3 != 0)
						itdbg3=1.0/invdt3;
					if (invdt4 != 0)
						itdbg4=1.0/invdt4;
					if (invdt5 != 0)
						itdbg5=1.0/invdt5;
					mdtdbg = dtGlobal;
					if ((parameters::artificial_viscosity == parameters::artificial_viscosity_SN) && (parameters::artificial_viscosity_factor > 0)) {
						viscRadial = dxRadial/dvRadial/4.0/pow2(parameters::artificial_viscosity_factor);
						viscAzimuthal = dxAzimuthal/dvAzimuthal/4.0/pow2(parameters::artificial_viscosity_factor);
					}
				}
			}
		}
	}

	//for (unsigned int n_radial = Zero_or_active; n_radial < MaxMO_or_active; ++n_radial) {
	for (unsigned int n_radial = 1+(CPUOVERLAP-1) * (CPU_Rank > 0 ? 1 : 0); n_radial < NRadial-2-(CPUOVERLAP-1)*(CPU_Rank < CPU_Number-1 ? 1 : 0); ++n_radial) {
        dtLocal = 2.0*PI*parameters::CFL/(double)NAzimuthal/fabs(v_mean[n_radial]*InvRmed[n_radial]-v_mean[n_radial+1]*InvRmed[n_radial+1]);

		if (dtLocal < dtGlobal)
			dtGlobal = dtLocal;
	}

	if (debug) {
		logging::print(LOG_DEBUG "Timestep control information for CPU %d: \n", CPU_Rank);
		logging::print(LOG_DEBUG "Most restrictive cell at nRadial=%d and nAzimuthal=%d\n", n_radial_debug, n_azimuthal_debug);
		logging::print(LOG_DEBUG "located at radius Rmed         : %g\n", Rmed[n_radial_debug]);
		logging::print(LOG_DEBUG "Sound speed limit              : %g\n", itdbg1);
		logging::print(LOG_DEBUG "Radial motion limit            : %g\n", itdbg2);
		logging::print(LOG_DEBUG "Residual circular motion limit : %g\n", itdbg3);

		if (parameters::artificial_viscosity_factor > 0) {
			logging::print(LOG_DEBUG "Articifial Viscosity limit     : %g\n", itdbg4);
			logging::print(LOG_DEBUG "   Arise from r with limit     : %g\n", viscRadial);
			logging::print(LOG_DEBUG "   and from theta with limit   : %g\n", viscAzimuthal);
		} else {
			logging::print(LOG_DEBUG "Articifial Viscosity limit     : disabled\n");
		}
		logging::print(LOG_DEBUG "Kinematic viscosity limit      : %g\n", itdbg5);
		logging::print(LOG_DEBUG "Limit time step for this cell  : %g\n", mdtdbg);
		logging::print(LOG_DEBUG "Limit time step adopted        : %g\n", dtGlobal);
		if (dtGlobal < mdtdbg) {
			logging::print(LOG_DEBUG "Discrepancy arise either from shear.\n");
			logging::print(LOG_DEBUG "or from the imposed DT interval.\n");
		}
	}

	return max(deltaT/dtGlobal,1.0);
}

/**
	computes soundspeed
*/
void compute_sound_speed(t_data &data, bool force_update)
{
	static double last_physicaltime_calculated = -1;

	if ((!force_update) && (last_physicaltime_calculated == PhysicalTime)) {
			return;
	}
	last_physicaltime_calculated = PhysicalTime;

	for (unsigned int n_radial = 0; n_radial <= data[t_data::SOUNDSPEED].get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::SOUNDSPEED].get_max_azimuthal(); ++n_azimuthal) {
			if (parameters::Adiabatic) {
                /*if ((ADIABATICINDEX*(ADIABATICINDEX-1.0)*data[t_data::ENERGY](n_radial, n_azimuthal)/data[t_data::DENSITY](n_radial, n_azimuthal)) < 0) {
                    logging::print("SoundSpeed calc: %g < 0! density=%g energy=%g in cell (%u,%u)\n",(ADIABATICINDEX*(ADIABATICINDEX-1.0)*data[t_data::ENERGY](n_radial, n_azimuthal)/data[t_data::DENSITY](n_radial, n_azimuthal)), data[t_data::DENSITY](n_radial, n_azimuthal), data[t_data::ENERGY](n_radial, n_azimuthal), n_radial, n_azimuthal);
                }*/
                data[t_data::SOUNDSPEED](n_radial, n_azimuthal) = sqrt( ADIABATICINDEX*(ADIABATICINDEX-1.0)*data[t_data::ENERGY](n_radial, n_azimuthal)/data[t_data::DENSITY](n_radial, n_azimuthal) );
			} else if(parameters::Polytropic)
			{
				data[t_data::SOUNDSPEED](n_radial, n_azimuthal) = sqrt(ADIABATICINDEX * constants::R / parameters::MU * data[t_data::TEMPERATURE](n_radial, n_azimuthal));
			} else { // isothermal
				// This follows from: cs/v_Kepler = H/r
				data[t_data::SOUNDSPEED](n_radial, n_azimuthal) = ASPECTRATIO_REF * sqrt(constants::G*M/Rb[n_radial])*pow(Rb[n_radial], FLARINGINDEX);
			}
		}
	}
}

/**
	computes aspect ratio
*/
void compute_aspect_ratio(t_data &data, bool force_update)
{
	static double last_physicaltime_calculated = -1;

	if ((!force_update) && (last_physicaltime_calculated == PhysicalTime)) {
			return;
	}
	last_physicaltime_calculated = PhysicalTime;

	for (unsigned int n_radial = 0; n_radial <= data[t_data::ASPECTRATIO].get_max_radial(); ++n_radial) {
		double inv_v_kepler = 1.0/(omega_kepler(Rb[n_radial])*Rb[n_radial]);

		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ASPECTRATIO].get_max_azimuthal(); ++n_azimuthal) {
			if (parameters::Adiabatic || parameters::Polytropic) {
                // h = H/r = c_s,iso / v_k = c_s/sqrt(gamma) / v_k
                data[t_data::ASPECTRATIO](n_radial, n_azimuthal) = data[t_data::SOUNDSPEED](n_radial, n_azimuthal)/(sqrt(ADIABATICINDEX))*inv_v_kepler;
			} else {
                // h = H/r = c_s/v_k
                data[t_data::ASPECTRATIO](n_radial, n_azimuthal) = data[t_data::SOUNDSPEED](n_radial, n_azimuthal)*inv_v_kepler;
			}
		}
	}
}

/**
	computes pressure
*/
void compute_pressure(t_data &data, bool force_update)
{
	static double last_physicaltime_calculated = -1;

	if ((!force_update) && (last_physicaltime_calculated == PhysicalTime)) {
			return;
	}
	last_physicaltime_calculated = PhysicalTime;


	for (unsigned int n_radial = 0; n_radial <= data[t_data::PRESSURE].get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::PRESSURE].get_max_azimuthal(); ++n_azimuthal) {
			if (parameters::Adiabatic) {
                data[t_data::PRESSURE](n_radial,n_azimuthal) = (ADIABATICINDEX-1.0)*data[t_data::ENERGY](n_radial,n_azimuthal);
			} else if(parameters::Polytropic)
			{
				data[t_data::PRESSURE](n_radial,n_azimuthal) = data[t_data::DENSITY](n_radial,n_azimuthal)*pow2(data[t_data::SOUNDSPEED](n_radial,n_azimuthal)) / ADIABATICINDEX;
			} else { // Isothermal
				// since SoundSpeed is not update from initialization, cs remains axisymmetric
				data[t_data::PRESSURE](n_radial,n_azimuthal) = data[t_data::DENSITY](n_radial,n_azimuthal)*pow2(data[t_data::SOUNDSPEED](n_radial,n_azimuthal));
			}
		}
	}
}

/**
	computes temperature
*/
void compute_temperature(t_data &data, bool force_update)
{
	static double last_physicaltime_calculated = -1;

	if ((!force_update) && (last_physicaltime_calculated == PhysicalTime)) {
			return;
	}
	last_physicaltime_calculated = PhysicalTime;


	for (unsigned int n_radial = 0; n_radial <= data[t_data::TEMPERATURE].get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::TEMPERATURE].get_max_azimuthal(); ++n_azimuthal) {
			if (parameters::Adiabatic) {
                data[t_data::TEMPERATURE](n_radial,n_azimuthal) = parameters::MU/constants::R*(ADIABATICINDEX-1.0)*data[t_data::ENERGY](n_radial,n_azimuthal)/data[t_data::DENSITY](n_radial,n_azimuthal);
			} else if(parameters::Polytropic){
				data[t_data::TEMPERATURE](n_radial,n_azimuthal) = parameters::MU/constants::R*POLYTROPIC_CONSTANT*pow(data[t_data::DENSITY](n_radial,n_azimuthal), ADIABATICINDEX-1.0);
			} else { // Isothermal
				data[t_data::TEMPERATURE](n_radial,n_azimuthal) = parameters::MU/constants::R*data[t_data::PRESSURE](n_radial,n_azimuthal)/data[t_data::DENSITY](n_radial,n_azimuthal);
			}
		}
	}
}

/**
	computes density rho
*/
void compute_rho(t_data &data, bool force_update)
{
	static double last_physicaltime_calculated = -1;

	if ((!force_update) && (last_physicaltime_calculated == PhysicalTime)) {
			return;
	}
	last_physicaltime_calculated = PhysicalTime;

	double H;

	compute_aspect_ratio(data, force_update);

	for (unsigned int n_radial = 0; n_radial <= data[t_data::RHO].get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::RHO].get_max_azimuthal(); ++n_azimuthal) {
			if (parameters::Adiabatic) {
				H = data[t_data::ASPECTRATIO](n_radial,n_azimuthal)*Rb[n_radial];
			} else {
				H = ASPECTRATIO_REF * pow(Rb[n_radial], 1.0+FLARINGINDEX);
			}
			data[t_data::RHO](n_radial,n_azimuthal) = data[t_data::DENSITY](n_radial,n_azimuthal)/(parameters::density_factor*H);
		}
	}
}
/**
	Calculates the gas mass inside the planet's Roche lobe
*/
double CircumPlanetaryMass(t_data &data)
{
	double xpl, ypl;
	double dist, mdcplocal, mdcptotal;
	double *abs, *ord;

	/* if there's no planet, there is no mass inside its Roche lobe ;) */
	if (data.get_planetary_system().get_number_of_planets() == 0)
		return 0;

	// TODO: non global
	abs = CellAbscissa->Field;
	ord = CellOrdinate->Field;

	xpl = data.get_planetary_system().get_planet(0).get_x();
	ypl = data.get_planetary_system().get_planet(0).get_y();

	mdcplocal = 0.0;
	mdcptotal = 0.0;

	for (unsigned int n_radial = Zero_or_active; n_radial < Max_or_active; ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal < data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
			unsigned int cell = n_radial*data[t_data::DENSITY].get_size_azimuthal() + n_azimuthal;
			dist = sqrt ( (abs[cell]-xpl)*(abs[cell]-xpl) +	(ord[cell]-ypl)*(ord[cell]-ypl) );
			if ( dist < HillRadius ) {
				mdcplocal += Surf[n_radial] * data[t_data::DENSITY](n_radial, n_azimuthal);
			}
		}
	}

	MPI_Allreduce(&mdcplocal, &mdcptotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return mdcptotal;
}
