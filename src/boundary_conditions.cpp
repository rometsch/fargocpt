/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/

#include "boundary_conditions.h"
#include "parameters.h"
#include "logging.h"
#include "global.h"
#include <math.h>
#include "util.h"
#include "Theo.h"
#include <cstring>

// temporary
#include "SideEuler.h"
#include "LowTasks.h"
extern boolean OuterSourceMass;
extern bool Adiabatic;

namespace boundary_conditions {

void apply_boundary_condition(t_data &data, double dt, bool final)
{
	// if this is the final boundary condition and damping is enable, do it
	if (final && parameters::damping_enabled) {
		damping(data, dt);
	}

	// inner boundary
	switch (parameters::boundary_inner) {
		case parameters::boundary_condition_open:
			open_boundary_inner(data);
			break;
		case parameters::boundary_condition_reflecting:
			reflecting_boundary_inner(data);
			break;
		case  parameters::boundary_condition_boundary_layer:
			boundary_layer_inner_boundary(data);
			break;
		case parameters::boundary_condition_nonreflecting:
			NonReflectingBoundary_inner(data, &data[t_data::V_RADIAL], &data[t_data::DENSITY], &data[t_data::ENERGY]);
			break;
		case parameters::boundary_condition_viscous_outflow:
			viscous_outflow_boundary_inner(data);
			break;
		case parameters::boundary_condition_evanescent: // evanescent works only for inner and outer together until now
			if (parameters::boundary_outer == parameters::boundary_condition_evanescent) {
				EvanescentBoundary(data, dt);
			} else {
				logging::print_master(LOG_ERROR "Different EvanescentBoundary Parameters. Old .par File?\n");
				die("inner/outer evanescent boundary not implemented yet");
			}
			break;
		case parameters::boundary_condition_keplerian:
			keplerian2d_boundary_inner(data);
			break;
	}

	// outer boundary
	switch (parameters::boundary_outer) {
		case parameters::boundary_condition_open:
			open_boundary_outer(data);
			break;
		case parameters::boundary_condition_reflecting:
			reflecting_boundary_outer(data);
			break;
		case  parameters::boundary_condition_boundary_layer:
			boundary_layer_outer_boundary(data);
			break;
		case parameters::boundary_condition_nonreflecting:
			NonReflectingBoundary_outer(data, &data[t_data::V_RADIAL], &data[t_data::DENSITY], &data[t_data::ENERGY]);
            break;
		case parameters::boundary_condition_viscous_outflow:
			die("outer viscous outflow boundary not implemented");
			break;
		case parameters::boundary_condition_evanescent:
			//EvanescentBoundary_outer(VRadial, VAzimuthal, Density, Energy, dt, sys);
			break;
		case parameters::boundary_condition_keplerian:
			keplerian2d_boundary_outer(data);
			break;

	}

	if (CPU_Rank == CPU_Highest) {
		if ((parameters::domegadr_zero) && (parameters::boundary_outer != parameters::boundary_condition_boundary_layer)) {
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::V_AZIMUTHAL].get_max_azimuthal(); ++n_azimuthal) {
				// this is a work around as long as V_AZIMUTHAL is defined as a vector
				if (data[t_data::V_AZIMUTHAL].is_vector()) {
					data[t_data::V_AZIMUTHAL](data[t_data::V_AZIMUTHAL].get_max_radial()-1, n_azimuthal) = Rmed[data[t_data::V_AZIMUTHAL].get_max_radial()-1]/Rmed[data[t_data::V_AZIMUTHAL].get_max_radial()-2]*data[t_data::V_AZIMUTHAL](data[t_data::V_AZIMUTHAL].get_max_radial()-2, n_azimuthal);
				} else {
					data[t_data::V_AZIMUTHAL](data[t_data::V_AZIMUTHAL].get_max_radial(), n_azimuthal) = Rmed[data[t_data::V_AZIMUTHAL].get_max_radial()]/Rmed[data[t_data::V_AZIMUTHAL].get_max_radial()-1]*data[t_data::V_AZIMUTHAL](data[t_data::V_AZIMUTHAL].get_max_radial()-1, n_azimuthal);
				}
			}
		} else {
			/* seems to be done by ApplySubKeplerianBoundary
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::V_AZIMUTHAL].get_max_azimuthal(); ++n_azimuthal) {
					data[t_data::V_AZIMUTHAL](data[t_data::V_AZIMUTHAL].get_max_radial(), n_azimuthal) = data[t_data::V_AZIMUTHAL0](data[t_data::V_AZIMUTHAL0].get_max_radial(), n_azimuthal);
			}
			*/
		}
	}



	if (OuterSourceMass)
		ApplyOuterSourceMass(&data[t_data::DENSITY], &data[t_data::V_RADIAL]);

	if (parameters::massoverflow)
		mass_overflow(data);
}

/**
	inner open boundary condition
*/
void open_boundary_inner(t_data &data)
{
	if (CPU_Rank != 0)
		return;

	for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
		// copy first ring into ghost ring
		data[t_data::DENSITY](0, n_azimuthal) = data[t_data::DENSITY](1, n_azimuthal);
		data[t_data::ENERGY](0, n_azimuthal) = data[t_data::ENERGY](1, n_azimuthal);

		// set velocity to min(v[+1],0) (allow only outflow)
		if (data[t_data::V_RADIAL](2, n_azimuthal) > 0.0) {
			data[t_data::V_RADIAL](1, n_azimuthal) = 0.0;
			data[t_data::V_RADIAL](0, n_azimuthal) = 0.0;
		} else {
			data[t_data::V_RADIAL](1, n_azimuthal) = data[t_data::V_RADIAL](2, n_azimuthal);
			data[t_data::V_RADIAL](0, n_azimuthal) = data[t_data::V_RADIAL](2, n_azimuthal);
		}
	}
}

/**
	outer open boundary condition
 */
void open_boundary_outer(t_data &data)
{
	if (CPU_Rank != CPU_Highest)
		return;

	for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
		// copy last ring into ghost ring
		data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial(), n_azimuthal) = data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial()-1, n_azimuthal);
		data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial(), n_azimuthal) = data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial()-1, n_azimuthal);

		// set velocity to min(v[+1],0) (allow only outflow)
		if (data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial()-2, n_azimuthal) < 0.0) {
			data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial()-1, n_azimuthal) = 0.0;
			data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial(), n_azimuthal) = 0.0;
		} else {
			data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial()-1, n_azimuthal) = data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial()-2, n_azimuthal);
			data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial(), n_azimuthal) = data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial()-2, n_azimuthal);
		}
	}
}

/**
	inner reflecting boundary condition
*/
void reflecting_boundary_inner(t_data &data)
{
	if (CPU_Rank == 0) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
			// copy first ring into ghost ring
			data[t_data::DENSITY](0, n_azimuthal) = data[t_data::DENSITY](1, n_azimuthal);
			data[t_data::ENERGY](0, n_azimuthal) = data[t_data::ENERGY](1, n_azimuthal);

			data[t_data::V_RADIAL](1, n_azimuthal) = 0.0;
			data[t_data::V_RADIAL](0, n_azimuthal) = -data[t_data::V_RADIAL](2, n_azimuthal);
		}
	}
}

/**
	outer reflecting boundary condition
*/
void reflecting_boundary_outer(t_data &data)
{
	if (CPU_Rank == CPU_Highest) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
			// copy last ring into ghost ring
			data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial(), n_azimuthal) = data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial()-1, n_azimuthal);
			data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial(), n_azimuthal) = data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial()-1, n_azimuthal);

			data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial()-1, n_azimuthal) = 0.0;
			data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial(), n_azimuthal) = -data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial()-2, n_azimuthal);
		}
	}
}

void viscous_outflow_boundary_inner(t_data &data)
{
	if (CPU_Rank == 0) {
		double nu_med = 0.0;

		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::VISCOSITY].get_max_azimuthal(); ++n_azimuthal) {
			nu_med += data[t_data::VISCOSITY](1, n_azimuthal);
			data[t_data::DENSITY](0, n_azimuthal) = data[t_data::DENSITY](1, n_azimuthal);
			data[t_data::ENERGY](0, n_azimuthal) = data[t_data::ENERGY](1, n_azimuthal);
		}

		nu_med /= data[t_data::DENSITY].get_size_azimuthal();

		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::V_RADIAL].get_max_azimuthal(); ++n_azimuthal) {
			// V_rad =  - 1.5 / r * Nu (Kley, Papaloizou and Ogilvie, 2008)
			data[t_data::V_RADIAL](1, n_azimuthal) = -1.5 / Rinf[1] * nu_med;
			data[t_data::V_RADIAL](0, n_azimuthal) = -1.5 / Rinf[0] * nu_med;
		}
	}
}

static void damping_single(t_polargrid &quantity, t_polargrid &quantity0, double dt) {
	// use the correct radius array corresponding to quantity
	t_radialarray &radius = quantity.is_scalar() ? Rb : Ra;
	double delta;
	double X, X0, Xnew;
	bool is_density = strcmp( quantity.get_name(), "dens");

	// is this CPU in the inner damping domain?
	if ((parameters::damping_inner_limit > 1.0) && (radius[0] < RMIN*parameters::damping_inner_limit)) {
		// find range
		unsigned int limit = 0;
		while ((radius[limit] < RMIN*parameters::damping_inner_limit) && (limit <= quantity.get_max_radial())) {
			limit++;
		}
		limit--;

		double tau = parameters::damping_time_factor*2.0*PI/omega_kepler(RMIN);

		for (unsigned int n_radial = 0; n_radial <= limit; ++n_radial) {
			double factor = pow2((radius[n_radial]-RMIN*parameters::damping_inner_limit)/(RMIN-RMIN*parameters::damping_inner_limit));
			double exp_factor = exp(-dt*factor/tau);

			for (unsigned int n_azimuthal = 0; n_azimuthal <= quantity.get_max_azimuthal(); ++n_azimuthal) {
				X = quantity(n_radial,n_azimuthal);
				X0 = quantity0(n_radial,n_azimuthal);
				Xnew = (X-X0) * exp_factor + X0;
				quantity(n_radial, n_azimuthal) = Xnew;
				delta = Xnew - X;
				if ( is_density ) {
					if (delta > 0) {
						MassDelta.WaveDampingPositive += delta*Surf[n_radial];
					} else {
						MassDelta.WaveDampingNegative += delta*Surf[n_radial];
					}
				}
			}
		}
	}

	// is this CPU in the outer damping domain?
	if ((parameters::damping_outer_limit < 1.0) && (radius[quantity.get_max_radial()] > RMAX*parameters::damping_outer_limit)) {
		// find range
		unsigned int limit = quantity.get_max_radial();
		while ((radius[limit] > RMAX*parameters::damping_outer_limit) && (limit > 0)) {
			limit--;
		}
		limit++;

		double tau = parameters::damping_time_factor*2.0*PI/omega_kepler(RMAX);

		for (unsigned int n_radial = limit; n_radial <= quantity.get_max_radial(); ++n_radial) {
			double factor = pow2((radius[n_radial]-RMAX*parameters::damping_outer_limit)/(RMAX-RMAX*parameters::damping_outer_limit));
			double exp_factor = exp(-dt*factor/tau);

			for (unsigned int n_azimuthal = 0; n_azimuthal <= quantity.get_max_azimuthal(); ++n_azimuthal) {
				X = quantity(n_radial,n_azimuthal);
				X0 = quantity0(n_radial,n_azimuthal);
				Xnew = (X-X0) * exp_factor + X0;
				quantity(n_radial, n_azimuthal) = Xnew;
				delta = Xnew - X;
				if ( is_density ) {
					if (delta > 0) {
						MassDelta.WaveDampingPositive += delta*Surf[n_radial];
					} else {
						MassDelta.WaveDampingNegative += delta*Surf[n_radial];
					}
				}
			}
		}
	}
}

/**
	damping of velocities
*/
void damping(t_data &data, double dt)
{
	if (parameters::damping_v_radial) {
		damping_single(data[t_data::V_RADIAL], data[t_data::V_RADIAL0], dt);
	}

	if (parameters::damping_v_azimuthal) {
		damping_single(data[t_data::V_AZIMUTHAL], data[t_data::V_AZIMUTHAL0], dt);
	}

	if (parameters::damping_surface_density) {
		damping_single(data[t_data::DENSITY], data[t_data::DENSITY0], dt);
	}

	if ((Adiabatic) && (parameters::damping_energy)) {
		damping_single(data[t_data::ENERGY], data[t_data::ENERGY0], dt);
	}
}

void mass_overflow(t_data &data)
{
	if (CPU_Rank != CPU_Highest)
		return;
	// get location of binary star
	if (parameters::mof_planet+1 > data.get_planetary_system().get_number_of_planets()) {
		logging::print_master(LOG_ERROR "Wrong Planet/Star for Mass Overflow specified! Old .par File?\n");
		die("Wrong Planet/Star for Mass Overflow specified! Old .par File?");
	}
	double xplanet = data.get_planetary_system().get_planet(parameters::mof_planet).get_x();
	double yplanet = data.get_planetary_system().get_planet(parameters::mof_planet).get_y();
	// get grid cell where binary star is nearest
	// atan2(y,x) is from -PI to PI
	double angle = atan2(yplanet,xplanet)/2/PI;
	if (angle < 0) {
		angle +=1.0;
	}
	int maxcells = data[t_data::DENSITY].get_size_azimuthal();
	int nearest_grid_cell = maxcells * angle; // nearest gridcell to binary star
	double sigma = parameters::mof_sigma;
	int number_of_cells = maxcells * 3 * sigma; // walk through 3 sigma, so we get 99.7% accuracy
	if (number_of_cells<1) number_of_cells = 1;
	double sigmabar = maxcells * sigma;

	double weight_factor = 0.0;
	double check = 0.0;
	int gridcell = 0;
	int maxradial = data[t_data::DENSITY].get_max_radial();
	double dist = data.get_planetary_system().get_planet(parameters::mof_planet).get_distance();
	double mass_stream = 0.0;
	for (int i = -number_of_cells; i <= number_of_cells; i++) {
		gridcell = (nearest_grid_cell + i + maxcells)%maxcells;

		// adapt gauss profile
		weight_factor = 1.0 / (sigmabar*sqrt(2.0*PI)) * exp(-1.0/2.0 * pow2(i/sigmabar));
		check += weight_factor;
		mass_stream = weight_factor * parameters::mof_value / 2.0 / PI * DT;

		// overflow is given in flow / year
		// scalar quantities maxradial -1, vectorial quantities maxrad-2
		data[t_data::DENSITY](maxradial-1,gridcell) += mass_stream;
		// calculate gas velocities
		data[t_data::V_RADIAL](maxradial-2,gridcell) = -1.0 * ( data[t_data::DENSITY](maxradial-1,gridcell) * data[t_data::V_RADIAL](maxradial-2,gridcell) + mass_stream * 0.1 * omega_kepler(dist) * Rmed[maxradial-1]) / (data[t_data::DENSITY](maxradial-1,gridcell) + mass_stream); // speed inwards

		data[t_data::V_AZIMUTHAL](maxradial-2,gridcell) = (omega_kepler(dist) - OmegaFrame) * Rmed[maxradial-1];

		#ifndef NDEBUG
		logging::print(LOG_VERBOSE "dens %lE, WF %lE , angle %lf, nearest_grid_cell %i, mass_stream %lE, gridcell %i, maxcells %i , noc %i , i %i \n",data[t_data::DENSITY](maxradial-1,gridcell), weight_factor, angle, nearest_grid_cell,  mass_stream, gridcell, maxcells, number_of_cells*2+1, i);
		#endif
	}

	// check if mass overflow is constant
		if (!(check > 0.99) || !(check < 1.01)) {
			logging::print_master(LOG_ERROR "weight_factor %lf should be 0.997 \n",weight_factor);

		}


}

void apply_boundary_condition_temperature(t_data &data) {
	if (CPU_Rank == 0) {
		for (unsigned int n_radial = 0; n_radial <= GHOSTCELLS_B; ++n_radial) {
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal(); ++n_azimuthal) {
				data[t_data::ENERGY](n_radial,n_azimuthal) = data[t_data::TEMPERATURE](n_radial,n_azimuthal)*data[t_data::DENSITY](n_radial,n_azimuthal)/(ADIABATICINDEX-1.0)/parameters::MU*constants::R;
			}
		}

		switch (parameters::boundary_inner) {
			case parameters::boundary_condition_open:
				open_boundary_inner(data);
				break;
			case parameters::boundary_condition_reflecting:
				reflecting_boundary_inner(data);
				break;
			default:
				logging::print_master(LOG_ERROR "Boundary condition at inner boundary not supported for temperature!\n");
				break;
		}

		for (unsigned int n_radial = 0; n_radial <= GHOSTCELLS_B; ++n_radial) {
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal(); ++n_azimuthal) {
				data[t_data::TEMPERATURE](n_radial,n_azimuthal) = data[t_data::ENERGY](n_radial,n_azimuthal)/data[t_data::DENSITY](n_radial,n_azimuthal)*(ADIABATICINDEX-1.0)*parameters::MU/constants::R;
			}
		}
	}

	if (CPU_Rank == CPU_Highest) {
		for (unsigned int n_radial = data[t_data::ENERGY].get_max_radial()-GHOSTCELLS_B; n_radial <= data[t_data::ENERGY].get_max_radial(); ++n_radial) {
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal(); ++n_azimuthal) {
				data[t_data::ENERGY](n_radial,n_azimuthal) = data[t_data::TEMPERATURE](n_radial,n_azimuthal)*data[t_data::DENSITY](n_radial,n_azimuthal)/(ADIABATICINDEX-1.0)/parameters::MU*constants::R;
			}
		}

		// outer boundary
		switch (parameters::boundary_outer) {
			case parameters::boundary_condition_open:
				open_boundary_outer(data);
				break;
			case parameters::boundary_condition_reflecting:
				reflecting_boundary_outer(data);
				break;
			default:
				logging::print_master(LOG_ERROR "Boundary condition at outer boundary not supported for temperature!\n");
				break;
		}

		for (unsigned int n_radial = data[t_data::ENERGY].get_max_radial()-GHOSTCELLS_B; n_radial <= data[t_data::ENERGY].get_max_radial(); ++n_radial) {
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal(); ++n_azimuthal) {
				data[t_data::TEMPERATURE](n_radial,n_azimuthal) = data[t_data::ENERGY](n_radial,n_azimuthal)/data[t_data::DENSITY](n_radial,n_azimuthal)*(ADIABATICINDEX-1.0)*parameters::MU/constants::R;
			}
		}
	}
}

/**
	Boundary conditions for calculation of the boundary layer (BL) starting here:
	TODO: Stellar radiative flux into disk via implicit routine!
*/

/**
	Inner boundary: zero gradient & fixed velocities
*/

void boundary_layer_inner_boundary(t_data &data) {

	if (CPU_Rank != 0)
		return;

	for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
		// zero gradient
		data[t_data::DENSITY](0, n_azimuthal) = data[t_data::DENSITY](1, n_azimuthal);
		data[t_data::ENERGY](0, n_azimuthal) = data[t_data::ENERGY](1, n_azimuthal);

		// set vrad to fraction of Keplerian velocity
		data[t_data::V_RADIAL](1, n_azimuthal) = -1. * parameters::vrad_fraction_of_kepler * sqrt(1./Ra[1]);
		data[t_data::V_RADIAL](0, n_azimuthal) = data[t_data::V_RADIAL](1, n_azimuthal);

		// set vphi to stellar rotation rate
		data[t_data::V_AZIMUTHAL](0, n_azimuthal) = parameters::stellar_rotation_rate * sqrt(1./Rb[0]);
	}
}

/**
	Outer boundary: floating boundary conditions & pressure correction for Omega
*/

void boundary_layer_outer_boundary(t_data &data) {

	if (CPU_Rank != CPU_Highest)
		return;

	for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
		// floating BCs
		data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial(), n_azimuthal) = data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial()-1, n_azimuthal) * sqrt(Ra[data[t_data::DENSITY].get_max_radial()-1]/Ra[data[t_data::DENSITY].get_max_radial()]);
		data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial(), n_azimuthal) = data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial()-1, n_azimuthal) * pow(Ra[data[t_data::ENERGY].get_max_radial()-1]/Ra[data[t_data::ENERGY].get_max_radial()], 1.25);

		data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial()-1, n_azimuthal) = -1.*fabs(data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial()-2, n_azimuthal))*sqrt(Ra[data[t_data::V_RADIAL].get_max_radial()-2]/Ra[data[t_data::V_RADIAL].get_max_radial()-1]);
		data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial(), n_azimuthal) = -1.*fabs(data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial()-2, n_azimuthal))*sqrt(Ra[data[t_data::V_RADIAL].get_max_radial()-2]/Ra[data[t_data::V_RADIAL].get_max_radial()]);

		// Omega at outer boundary equals omega_kepler (plus leading order pressure correction)
		data[t_data::V_AZIMUTHAL](data[t_data::V_AZIMUTHAL].get_max_radial(), n_azimuthal) = 1./sqrt(Rb[data[t_data::DENSITY].get_max_radial()]);
		// TODO: Include pressure correction, like in uphi[*jN] = 1./sqrt(Rb[*jN]) + 0.5/Sigma[*jN]*sqrt(pow3(Rb[*jN])*pow2(Rb[*jN]))*.5/Rb[*jN]*(P[*jN+1]-P[*jN-1])/DeltaRa[*jN+1];
	}
}

/**
	End BL Boundary Conditions!
*/


void keplerian2d_boundary_inner(t_data &data) {
	for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal(); ++n_azimuthal) {
					data[t_data::DENSITY](1, n_azimuthal) = parameters::sigma0*pow(Rmed[1],-SIGMASLOPE);
					data[t_data::DENSITY](0, n_azimuthal) = parameters::sigma0*pow(Rmed[0],-SIGMASLOPE);
					data[t_data::ENERGY](1, n_azimuthal) = 1.0/(ADIABATICINDEX-1.0)*parameters::sigma0*pow2(ASPECTRATIO)*pow(Rmed[1],-SIGMASLOPE-1.0+2.0*FLARINGINDEX);
					data[t_data::ENERGY](0, n_azimuthal) = 1.0/(ADIABATICINDEX-1.0)*parameters::sigma0*pow2(ASPECTRATIO)*pow(Rmed[0],-SIGMASLOPE-1.0+2.0*FLARINGINDEX);
					data[t_data::TEMPERATURE](1, n_azimuthal) = data[t_data::ENERGY](1, n_azimuthal)/data[t_data::DENSITY](1, n_azimuthal)*(ADIABATICINDEX-1.0)*parameters::MU*constants::R;
					data[t_data::TEMPERATURE](0, n_azimuthal) = data[t_data::ENERGY](0, n_azimuthal)/data[t_data::DENSITY](0, n_azimuthal)*(ADIABATICINDEX-1.0)*parameters::MU*constants::R;
					data[t_data::V_RADIAL](1 ,n_azimuthal) = 0.0;
					data[t_data::V_RADIAL](0, n_azimuthal) = -data[t_data::V_RADIAL](2, n_azimuthal);
					data[t_data::V_AZIMUTHAL](1, n_azimuthal) = sqrt(constants::G*M/Rmed[1]);
					data[t_data::V_AZIMUTHAL](0, n_azimuthal) = sqrt(constants::G*M/Rmed[0]);
	}
}

void keplerian2d_boundary_outer(t_data &data) {
	for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal(); ++n_azimuthal) {
					data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial(), n_azimuthal) = parameters::sigma0*pow(Rmed[data[t_data::DENSITY].get_max_radial()],-SIGMASLOPE);
					data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial()-1, n_azimuthal) = parameters::sigma0*pow(Rmed[data[t_data::DENSITY].get_max_radial()-1],-SIGMASLOPE);
					data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial(), n_azimuthal) = 1.0/(ADIABATICINDEX-1.0)*parameters::sigma0*pow2(ASPECTRATIO)*pow(Rmed[data[t_data::DENSITY].get_max_radial()],-SIGMASLOPE-1.0+2.0*FLARINGINDEX);
					data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial()-1, n_azimuthal) = 1.0/(ADIABATICINDEX-1.0)*parameters::sigma0*pow2(ASPECTRATIO)*pow(Rmed[data[t_data::DENSITY].get_max_radial()-1],-SIGMASLOPE-1.0+2.0*FLARINGINDEX);
					data[t_data::TEMPERATURE](data[t_data::TEMPERATURE].get_max_radial(), n_azimuthal) = data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial(), n_azimuthal)/data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial(), n_azimuthal)*(ADIABATICINDEX-1.0)*parameters::MU*constants::R;
					data[t_data::TEMPERATURE](data[t_data::TEMPERATURE].get_max_radial()-1, n_azimuthal) = data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial()-1, n_azimuthal)/data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial()-1, n_azimuthal)*(ADIABATICINDEX-1.0)*parameters::MU*constants::R;
					data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial(), n_azimuthal) = -data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial()-2, n_azimuthal);
					data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial()-1,n_azimuthal) = 0.0;
					data[t_data::V_AZIMUTHAL](data[t_data::V_AZIMUTHAL].get_max_radial(), n_azimuthal) = sqrt(constants::G*M/Rmed[data[t_data::DENSITY].get_max_radial()]);
					data[t_data::V_AZIMUTHAL](data[t_data::V_AZIMUTHAL].get_max_radial()-1, n_azimuthal) = sqrt(constants::G*M/Rmed[data[t_data::DENSITY].get_max_radial()-1]);
	}
}

}
