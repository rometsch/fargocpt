#include <stdio.h>
#include <float.h>
#include <math.h>
#include <assert.h>

#include "parameters.h"
#include "selfgravity.h"
#include "viscosity.h"
#include "SideEuler.h"
#include "Pframeforce.h"
#include "global.h"
#include "init.h"
#include "nongnu.h"
#include "logging.h"
#include "SourceEuler.h"
#include "Theo.h"
#include "LowTasks.h"
#include "output.h"
#include "axilib.h"
#include "util.h"
#include "quantities.h"
#include "TransportEuler.h"
#include "boundary_conditions.h"

#include "options.h"
#include "open-simplex-noise.h"

extern int Restart;
extern boolean Corotating;

/**
	resize all (global) radialarrays
*/
void resize_radialarrays(unsigned int size)
{
	EnergyMed.resize(size);
	SigmaMed.resize(size);
	Rmed.resize(size);
	InvRmed.resize(size);

	// these are size+1 so that we always can use Rinf[]
	// e.g.: Rinf[max_radial] instead of Rsup[max_radial-1]
	Rinf.resize(size+1);
	InvRinf.resize(size+1);
	Rsup.resize(size);
	Surf.resize(size);
	InvSurf.resize(size);
	InvDiffRmed.resize(size);
	InvDiffRsup.resize(size);
	Radii.resize(size);
	GlobalRmed.resize(size);
	SigmaInf.resize(size);
	GLOBAL_bufarray.resize(size);
}

/**
	Fills 1D arrays like Rmed, Rinf, Rsup, Surf, ...
*/
void init_radialarrays()
{
	FILE* fd_input;
	char* fd_input_filename;
	unsigned int nRadial;

	if (asprintf(&fd_input_filename, "%s%s", OUTPUTDIR, "radii.dat") == -1) {
		logging::print_master(LOG_ERROR "Not enough memory for string buffer.\n");
		PersonalExit(1);
	}

	fd_input = fopen(fd_input_filename, "r");

	// free up fdInputFilename
	free(fd_input_filename);

	if (fd_input == NULL) {
		logging::print_master(LOG_INFO "Warning : no `radii.dat' file found. Using default.\n");
		double interval, slope;
		switch (parameters::radial_grid_type) {
			case parameters::logarithmic_spacing:
				for (nRadial = 0; nRadial <= GlobalNRadial; ++nRadial) {
					Radii[nRadial] = RMIN*exp((double)(nRadial-1.0)/(double)(GlobalNRadial-2.0)*log(RMAX/RMIN));
				}
				break;
			case parameters::arithmetic_spacing:
				interval = (RMAX-RMIN)/(double)(GlobalNRadial-2.0);
				for (nRadial = 0; nRadial <= GlobalNRadial; ++nRadial) {
					Radii[nRadial] = RMIN+interval*(double)(nRadial-1.0);
				}
				break;
			case parameters::exponential_spacing:
				slope = ( .011*5. - .01 )*100./((double)(GlobalNRadial-2));
				for (nRadial = 0; nRadial <= GlobalNRadial; ++nRadial) {
					Radii[nRadial] = (RMIN-RMAX)/(exp(slope)-exp(slope*((double)GlobalNRadial-1)))*exp(slope*(double)(nRadial)) + RMIN - (RMIN - RMAX)/(1.-exp(slope*((double)GlobalNRadial-2)));
				}
				break;
			default:
				die("Invalid setting for RadialSpacing");
		}
	} else {
		logging::print_master(LOG_INFO "Reading 'radii.dat' file.\n");
		for (nRadial = 0; nRadial <= GlobalNRadial; ++nRadial) {
			double temp;

			if (fscanf(fd_input, "%lf", &temp)==1) {
				Radii[nRadial] = temp;
			} else {
				logging::print_master(LOG_ERROR "Reading 'radii.dat' file: No data left to read :(\n");
				PersonalExit(1);
			}
		}
	}

	/* if input file is open, close it */
	if (fd_input != NULL)
		fclose(fd_input);

	for (nRadial = 0; nRadial < GlobalNRadial; ++nRadial) {
		// Rmed is in the center of the cell where the center of mass is
		// Rmed = 1/2 * [ (4/3 Pi r_sup^3) - (4/3 Pi r_inf^3) ] / [ (Pi r_sup^2) - (Pi r_inf^2) ]
		GlobalRmed[nRadial] = 2.0/3.0*(Radii[nRadial+1]*Radii[nRadial+1]*Radii[nRadial+1]-Radii[nRadial]*Radii[nRadial]*Radii[nRadial]);
		GlobalRmed[nRadial] = GlobalRmed[nRadial] / (Radii[nRadial+1]*Radii[nRadial+1]-Radii[nRadial]*Radii[nRadial]);
	}

	logging::print_master(LOG_VERBOSE "Active %s grid is ranging from %g to %g. Total grid is range from %g to %g.\n", parameters::radial_grid_names[parameters::radial_grid_type] ,Radii[1],Radii[GlobalNRadial-1],Radii[0],Radii[GlobalNRadial]);

	for (nRadial = 0; nRadial < NRadial; ++nRadial) {
		Rinf[nRadial] = Radii[nRadial+IMIN];
		Rsup[nRadial] = Radii[nRadial+IMIN+1];

		Rmed[nRadial] = 2.0/3.0*(Rsup[nRadial]*Rsup[nRadial]*Rsup[nRadial]-Rinf[nRadial]*Rinf[nRadial]*Rinf[nRadial]);
		Rmed[nRadial] = Rmed[nRadial] / (Rsup[nRadial]*Rsup[nRadial]-Rinf[nRadial]*Rinf[nRadial]);

		// TODO: Is already calculated a few lines above. assert should check this
		assert( (Rmed[nRadial] - GlobalRmed[nRadial+IMIN]) < DBL_EPSILON );

		Surf[nRadial] = PI*(pow2(Rsup[nRadial])-pow2(Rinf[nRadial]))/(double)NAzimuthal;

		InvRmed[nRadial] = 1.0/Rmed[nRadial];
		InvSurf[nRadial] = 1.0/Surf[nRadial];
		InvDiffRsup[nRadial] = 1.0/(Rsup[nRadial]-Rinf[nRadial]);
		InvRinf[nRadial] = 1.0/Rinf[nRadial];
	}

	Rinf[NRadial]=Radii[NRadial+IMIN];
	InvRinf[NRadial]= 1.0/Rinf[NRadial];

	for (nRadial = 1; nRadial < NRadial; ++nRadial) {
		InvDiffRmed[nRadial] = 1.0/(Rmed[nRadial]-Rmed[nRadial-1]);
	}

	/* output radii to used_rad.dat (on master only) */
	if (CPU_Master) {
		FILE* fd_output;
		char* fd_output_filename;

		if (asprintf(&fd_output_filename, "%s%s", OUTPUTDIR, "used_rad.dat") == -1) {
			logging::print_master(LOG_ERROR "Not enough memory for string buffer.\n");
			PersonalExit(1);
		}

		fd_output = fopen(fd_output_filename, "w");

		if (fd_output == NULL) {
			logging::print_master(LOG_ERROR "Can't write %s.\nProgram stopped.\n", fd_output_filename);
			PersonalExit(1);
		}
		// free up fdOutputFilename
		free(fd_output_filename);

		for (nRadial = 0; nRadial <= GlobalNRadial; ++nRadial) {
			fprintf(fd_output, "%.18g\n", Radii[nRadial]);
		}

		/* if output file is open, close it */
		if (fd_output != NULL)
			fclose(fd_output);
	}
    MPI_Barrier(MPI_COMM_WORLD);

}

/**
	Initialises all physics
*/

void init_physics(t_data &data) {
	double foostep = 0.;

	if( (parameters::sigma_initialize_condition == parameters::initialize_condition_shakura_sunyaev) && (parameters::energy_initialize_condition == parameters::initialize_condition_shakura_sunyaev) ) {
		init_shakura_sunyaev(data);
		return;
	} else if( (parameters::sigma_initialize_condition == parameters::initialize_condition_shakura_sunyaev) || (parameters::energy_initialize_condition == parameters::initialize_condition_shakura_sunyaev) ) {
		die("Both Sigma and Energy have to be initialised by Shakura & Sunyaev Standard-Solution. Other initialisation not yet implemented!");
	}

	// gas density initialization
	init_gas_density(data);

	// if energy equation is taken into account, we initialize the gas thermal energy
	if (parameters::Adiabatic) {
		init_gas_energy(data);
	}

	if (parameters::self_gravity) {
		// if SelfGravity = YES or Z, planets are initialized feeling disk potential. Only the surface density is required to calculate the radial self-gravity acceleration. The disk radial and azimutal velocities are not updated
    selfgravity::init(data);
		selfgravity::init_planetary_system(data);
		logging::print_master(LOG_INFO "sg initialised\n");
	}

	// ListPlanets(sys);
	data.get_planetary_system().list_planets();

	OmegaFrame = OMEGAFRAME;

	if (Corotating) {
		OmegaFrame = data.get_planetary_system().get_planet(0).get_omega();
	}

	FrameAngle = 0;

	// only gas velocities remain to be initialized
	init_euler(data);
	init_gas_velocities(data);
	boundary_conditions::apply_boundary_condition(data, 0.0, false);
}

/**
	Wrapper for initialisation of physics according to Shakura & Sunyaev 1974
*/

void init_shakura_sunyaev(t_data &data) {
	double factor;

	if( !parameters::Adiabatic )
		die("Isothermal equation of state and Shakura & Sunyaev starting conditions has not yet been implemented!");

	for( unsigned int n_radial = 0; n_radial < data[t_data::TEMPERATURE].Nrad; ++n_radial ) {
		for( unsigned int n_azimuthal = 0; n_azimuthal < data[t_data::TEMPERATURE].Nsec; ++n_azimuthal ) {

			factor = pow(1. - sqrt(parameters::star_radius/(Rb[n_radial] + 2.*(RMAX - RMIN)/GlobalNRadial )), 0.25);

			data[t_data::DENSITY](n_radial,n_azimuthal) = (5.2*pow(ALPHAVISCOSITY,-4./5.)*pow(parameters::mass_accretion_rate*units::mass_accretion_rate.get_cgs_factor()/1.e16, 7./10.)*pow(parameters::M0, 0.25)*pow(Rb[n_radial]*units::length.get_cgs_factor()/1.e10,-0.75)*pow(factor,14./5.))/units::surface_density.get_cgs_factor();
			data[t_data::ASPECTRATIO](n_radial,n_azimuthal) = (1.7e8*pow(ALPHAVISCOSITY,-1./10.)*pow(parameters::mass_accretion_rate*units::mass_accretion_rate.get_cgs_factor()/1.e16,3./20.)*pow(parameters::M0,-3./8.)*pow(Rb[n_radial]*units::length.get_cgs_factor()/1.e10,9./8.)*pow(factor,3./5.))/(Rb[n_radial]*units::length.get_cgs_factor());
			data[t_data::TEMPERATURE](n_radial,n_azimuthal) = (1.4e4*pow(ALPHAVISCOSITY,-1./5.)*pow(parameters::mass_accretion_rate*units::mass_accretion_rate.get_cgs_factor()/1.e16, 3./10.)*pow(parameters::M0, 0.25)*pow(Rb[n_radial]*units::length.get_cgs_factor()/1.e10,-0.75)*pow(factor,6./5.))/units::temperature.get_cgs_factor();
			data[t_data::V_RADIAL](n_radial,n_azimuthal) = -(2.7e4*pow(ALPHAVISCOSITY,4./5.)*pow(parameters::mass_accretion_rate*units::mass_accretion_rate.get_cgs_factor()/1.e16, 3./10.)*pow(parameters::M0, -0.25)*pow(Rb[n_radial]*units::length.get_cgs_factor()/1.e10,-0.25)*pow(factor,-14./5.))/units::velocity.get_cgs_factor();

			data[t_data::SOUNDSPEED](n_radial,n_azimuthal) = sqrt(constants::R/parameters::MU * ADIABATICINDEX * data[t_data::TEMPERATURE](n_radial, n_azimuthal));
			data[t_data::ENERGY](n_radial, n_azimuthal) = constants::R/parameters::MU * 1./(ADIABATICINDEX-1.)*data[t_data::DENSITY](n_radial, n_azimuthal)*data[t_data::TEMPERATURE](n_radial, n_azimuthal);
			data[t_data::PRESSURE](n_radial, n_azimuthal) = (ADIABATICINDEX-1.)*data[t_data::ENERGY](n_radial, n_azimuthal);
			data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) = omega_kepler(Rb[n_radial])*Rb[n_radial];

			// if( Rb[n_radial] < parameters::star_radius ){
			// 	data[t_data::TEMPERATURE](n_radial, n_azimuthal) = pow( 3.*constants::G*M*parameters::mass_accretion_rate / (8.*PI*pow(Rb[n_radial],3)*constants::sigma.get_code_value() ) * (1. - pow(parameters::star_radius/Rb[n_radial+1],.5)), 0.25);
			// 	data[t_data::SOUNDSPEED](n_radial, n_azimuthal) = sqrt( constants::R/parameters::MU * ADIABATICINDEX * data[t_data::TEMPERATURE](n_radial, n_azimuthal) );
			// 	data[t_data::ASPECTRATIO](n_radial, n_azimuthal) = sqrt( constants::R/parameters::MU * data[t_data::TEMPERATURE](n_radial, n_azimuthal) ) / (omega_kepler(Rb[n_radial]) * Rb[n_radial]);
			// 	data[t_data::DENSITY](n_radial, n_azimuthal) = parameters::mass_accretion_rate/(3.*PI*ALPHAVISCOSITY*pow(data[t_data::SOUNDSPEED](n_radial, n_azimuthal),2)/omega_kepler(Rb[n_radial]))*(1. - pow(parameters::star_radius/Rb[n_radial+1],.5));
			// 	data[t_data::ENERGY](n_radial, n_azimuthal) = constants::R/parameters::MU * 1./(ADIABATICINDEX-1.)*data[t_data::DENSITY](n_radial, n_azimuthal)*data[t_data::TEMPERATURE](n_radial, n_azimuthal);
			// 	data[t_data::PRESSURE](n_radial, n_azimuthal) = (ADIABATICINDEX-1.)*data[t_data::ENERGY](n_radial, n_azimuthal);
			// 	data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) = omega_kepler(Rb[n_radial])*Rb[n_radial];
			// 	data[t_data::V_RADIAL](n_radial, n_azimuthal) = - parameters::mass_accretion_rate/(2.*PI*data[t_data::DENSITY](n_radial, n_azimuthal)*Ra[n_radial]);
			// } else {
			// 	data[t_data::TEMPERATURE](n_radial, n_azimuthal) = pow( 3.*constants::G*M*parameters::mass_accretion_rate / (8.*PI*pow(Rb[n_radial],3)*constants::sigma.get_code_value() ) * (1. - pow(parameters::star_radius/Rb[n_radial],.5)), 0.25);
			// 	data[t_data::SOUNDSPEED](n_radial, n_azimuthal) = sqrt( constants::R/parameters::MU * ADIABATICINDEX * data[t_data::TEMPERATURE](n_radial, n_azimuthal) );
			// 	data[t_data::ASPECTRATIO](n_radial, n_azimuthal) = sqrt( constants::R/parameters::MU * data[t_data::TEMPERATURE](n_radial, n_azimuthal) ) / (omega_kepler(Rb[n_radial]) * Rb[n_radial]);
			// 	data[t_data::DENSITY](n_radial, n_azimuthal) = parameters::mass_accretion_rate/(3.*PI*ALPHAVISCOSITY*pow(data[t_data::SOUNDSPEED](n_radial, n_azimuthal),2)/omega_kepler(Rb[n_radial]))*(1. - pow(parameters::star_radius/Rb[n_radial],.5));
			// 	//data[t_data::DENSITY](n_radial, n_azimuthal) = parameters::mass_accretion_rate/(3.*PI*ALPHAVISCOSITY*data[t_data::SOUNDSPEED](n_radial, n_azimuthal)*data[t_data::ASPECTRATIO](n_radial, n_azimuthal)*Rb[n_radial])*(1. - pow(parameters::star_radius/Rb[n_radial],.5));
			// 	data[t_data::ENERGY](n_radial, n_azimuthal) = constants::R/parameters::MU * 1./(ADIABATICINDEX-1.)*data[t_data::DENSITY](n_radial, n_azimuthal)*data[t_data::TEMPERATURE](n_radial, n_azimuthal);
			// 	data[t_data::PRESSURE](n_radial, n_azimuthal) = (ADIABATICINDEX-1.)*data[t_data::ENERGY](n_radial, n_azimuthal);
			// 	data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) = omega_kepler(Rb[n_radial])*Rb[n_radial];
			// 	data[t_data::V_RADIAL](n_radial, n_azimuthal) = - parameters::mass_accretion_rate/(2.*PI*data[t_data::DENSITY](n_radial, n_azimuthal)*Ra[n_radial]);
			// }
		}
	}

	RefillSigma(&data[t_data::DENSITY]);
	RefillEnergy(&data[t_data::ENERGY]);

	if( parameters::self_gravity )
		die("Self-gravity and Shakura-Sunyaev starting values has not yet been implemented!");

	/** init_euler w/o updates that have already been done above **/
	InitTransport();

	InitComputeAccel();

	viscosity::update_viscosity(data);
	/** end init_euler **/

	if( CentrifugalBalance )
		die("CentrifugalBalance and Shakura-Sunyaev starting values has not yet been implemented!");
}


/**
	Initializes gas density.
*/
void init_gas_density(t_data &data)
{
	switch (parameters::sigma_initialize_condition) {
		case parameters::initialize_condition_profile:
			logging::print_master(LOG_INFO "Initializing Sigma(r) = %g = %g %s * [r/(%g AU)]^(%g)\n", parameters::sigma0, parameters::sigma0*units::surface_density.get_cgs_factor(), units::surface_density.get_cgs_symbol(), units::length.get_cgs_factor()/units::cgs_AU, -SIGMASLOPE);

			for (unsigned int n_radial = 0; n_radial < data[t_data::DENSITY].Nrad; ++n_radial) {
				for (unsigned int n_azimuthal = 0; n_azimuthal < data[t_data::DENSITY].Nsec; ++n_azimuthal) {
					data[t_data::DENSITY](n_radial, n_azimuthal) = parameters::sigma0*pow(Rmed[n_radial],-SIGMASLOPE);
				}
			}
			break;

		case parameters::initialize_condition_read1D:
			logging::print_master(LOG_INFO "Loading Sigma from '%s' (1D).\n", parameters::sigma_filename);
			data[t_data::DENSITY].read1D(parameters::sigma_filename,true);
			break;

		case parameters::initialize_condition_read2D:
			logging::print_master(LOG_INFO "Loading Sigma from '%s' (2D).\n", parameters::sigma_filename);
			data[t_data::DENSITY].read2D(parameters::sigma_filename);
			break;

		case parameters::initialize_condition_shakura_sunyaev:
			die("Bad choice!");	//TODO: better explanation!
			break;
// 		case parameters::initialize_condition_shakura_sunyaev:
// 			logging::print_master(LOG_INFO "Initializing Sigma from Shakura and Sunyaev 1973 standard solution (cf. A&A, 24, 337)");
//
// 			for (unsigned int n_radial = 0; n_radial < data[t_data::DENSITY].Nrad; ++n_radial) {
// 				for (unsigned int n_azimuthal = 0; n_azimuthal < data[t_data::DENSITY].Nsec; ++n_azimuthal) {
// 					data[t_data::DENSITY](n_radial, n_azimuthal) =
// 					data[t_data::DENSITY](n_radial, n_azimuthal) = parameters::sigma0*pow(Rmed[n_radial],-SIGMASLOPE);
// 				}
// 			}
// 			break;
	}

	if (parameters::sigma_randomize) {

        struct osn_context *osn;
		if(open_simplex_noise(parameters::random_seed + 325582 , &osn) != 0) // add random number to seed to keep seed unique, because it is reused several times in the code.
        {
            die("Bad open simplex noise!");
        }


		logging::print_master(LOG_INFO "Randomizing Sigma by %.2f %%.\n",parameters::sigma_random_factor*100);
		for (unsigned int n_radial = 0; n_radial <= data[t_data::DENSITY].get_max_radial(); ++n_radial) {
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {

                double r = Rmed[n_radial];
                double angle = (double)n_azimuthal/(double)data[t_data::V_RADIAL].get_size_azimuthal()*2.0*PI;
                double x = r*cos(angle);
                double y = r*sin(angle);


                double f = parameters::sigma_feature_size;

                int simplex_order = 11;
                double noise = 0.0;
                for(int i = 0; i < simplex_order; ++i)
                {
                    double feature_factor= double(1<<i);
                    double weight_factor = double(1<<(simplex_order-i-1));

                    noise += open_simplex_noise2(osn, feature_factor*x/f, feature_factor*y/f) * weight_factor;
                }

                noise /= double((1<<simplex_order)-1);

                data[t_data::DENSITY](n_radial, n_azimuthal) *= (1 + parameters::sigma_random_factor*noise);
			}
		}

        open_simplex_noise_free(osn);
	}

	// profile damping?
	if (parameters::profile_damping) {
		logging::print_master(LOG_INFO "Damping Sigma for r > %g %s over a range from %g %s\n",parameters::profile_damping_point*units::length.get_cgs_factor(),units::length.get_cgs_symbol(),parameters::profile_damping_width*units::length.get_cgs_factor(),units::length.get_cgs_symbol());
		for (unsigned int n_radial = 0; n_radial <= data[t_data::DENSITY].get_max_radial(); ++n_radial) {
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
				// damp density to 0 for r > profile_damping_point
				data[t_data::DENSITY](n_radial, n_azimuthal) *= cutoff(parameters::profile_damping_point, parameters::profile_damping_width,Rmed[n_radial]);
				// set density to SIGMA0*SIGMA_FLOOR for r > r > profile_damping_point
				//data[t_data::DENSITY](n_radial, n_azimuthal) += (1-cutoff(parameters::profile_damping_point+parameters::profile_damping_width, parameters::profile_damping_width,Rmed[n_radial]))* parameters::sigma0 * parameters::sigma_floor;
			}
		}
	}

	// renormalize sigma0?
	if (parameters::sigma_adjust) {
		double total_mass = quantities::gas_total_mass(data);
		parameters::sigma0 *= parameters::sigma_discmass/total_mass;
		logging::print_master(LOG_INFO "Setting Sigma0=%g %s to set disc mass of %g to %g.\n",parameters::sigma0*units::surface_density.get_cgs_factor(),units::surface_density.get_cgs_symbol(),total_mass,parameters::sigma_discmass);

		// update density grid
		for (unsigned int n_radial = 0; n_radial <= data[t_data::DENSITY].get_max_radial(); ++n_radial) {
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
				data[t_data::DENSITY](n_radial, n_azimuthal) *= parameters::sigma_discmass/total_mass;
			}
		}
	} else {
		double total_mass = quantities::gas_total_mass(data);
		logging::print_master(LOG_INFO "Total disk is mass is %g = %g %s.\n", total_mass, total_mass*units::mass.get_cgs_factor(),units::mass.get_cgs_symbol());
	}

	// set SigmaMed/SigmaInf
	RefillSigma(&data[t_data::DENSITY]);
}

/**
	Initializes gas energy.
*/
void init_gas_energy(t_data &data)
{
	if (ADIABATICINDEX == 1.0) {
		logging::print_master(LOG_ERROR "The adiabatic index must differ from unity to initialize the gas internal energy. I must exit.\n");
		PersonalExit(1);
	}

	switch (parameters::energy_initialize_condition) {
		case parameters::initialize_condition_profile:
			logging::print_master(LOG_INFO "Initializing Energy=%g %s * [r/(%.1f AU)]^(%g). Flaring index is %g. T=%g %s * [r/(%.1f AU)]^(%g).\n", 1.0/((ADIABATICINDEX-1.0))*parameters::sigma0*pow2(ASPECTRATIO_REF)*units::energy.get_cgs_factor(),units::energy.get_cgs_symbol(), units::length.get_cgs_factor()/units::cgs_AU, -SIGMASLOPE-1.0+2.0*FLARINGINDEX, FLARINGINDEX, parameters::MU/constants::R*  pow2(ASPECTRATIO_REF)*units::temperature.get_cgs_factor(),units::temperature.get_cgs_symbol(), units::length.get_cgs_factor()/units::cgs_AU, -1.0+2.0*FLARINGINDEX);

			for (unsigned int n_radial = 0; n_radial <= data[t_data::ENERGY].get_max_radial(); ++n_radial) {
				for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal(); ++n_azimuthal) {
					data[t_data::ENERGY](n_radial,n_azimuthal) = 1.0/(ADIABATICINDEX-1.0)*parameters::sigma0*pow2(ASPECTRATIO_REF)*pow(Rmed[n_radial],-SIGMASLOPE-1.0+2.0*FLARINGINDEX);
				}
			}
			break;

		case parameters::initialize_condition_read1D:
			logging::print_master(LOG_INFO "Loading Energy from '%s' (1D).\n", parameters::energy_filename);
			data[t_data::ENERGY].read1D(parameters::energy_filename,true);
			break;

		case parameters::initialize_condition_read2D:
			logging::print_master(LOG_INFO "Loading Energy from '%s' (2D).\n", parameters::energy_filename);
			data[t_data::ENERGY].read2D(parameters::energy_filename);
			break;

		case parameters::initialize_condition_shakura_sunyaev:
			die("Bad choice!"); //TODO: better explanation!
			break;
	}

	// profile damping?
	if (parameters::profile_damping) {
		logging::print_master(LOG_INFO "Damping Energy for r > %g %s over a range of %g %s\n",parameters::profile_damping_point*units::length.get_cgs_factor(),units::length.get_cgs_symbol(),parameters::profile_damping_width*units::length.get_cgs_factor(),units::length.get_cgs_symbol());
		for (unsigned int n_radial = 0; n_radial <= data[t_data::ENERGY].get_max_radial(); ++n_radial) {
			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal(); ++n_azimuthal) {
				// damp energy to 0 for r > profile_damping_point
				data[t_data::ENERGY](n_radial, n_azimuthal) *= cutoff(parameters::profile_damping_point, parameters::profile_damping_width,Rmed[n_radial]);
			}
		}
	}

	// set EnergyMed
	RefillEnergy(&data[t_data::ENERGY]);
}

/**
	Initialize gas velocities. This is described in 3.1.1 of Clement Baruteau's PhD Thesis
*/
void init_gas_velocities(t_data &data)
{
	double r, ri;
	double t1, t2, r1, r2;
	double vt_cent[MAX1D];

	// Check if pure keplerian initialization is set
	if (parameters::initialize_pure_keplerian) {
		for (unsigned int n_radial = 0; n_radial <= data[t_data::V_AZIMUTHAL].get_max_radial(); ++n_radial) {
			// TODO: This should be Rinf but seems to produce incorrect data
			if (n_radial == data[t_data::V_AZIMUTHAL].Nrad) {
				r = Rmed[data[t_data::V_AZIMUTHAL].Nrad-1];
			} else {
				r = Rmed[n_radial];
			}

			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::V_AZIMUTHAL].get_max_azimuthal(); ++n_azimuthal) {
				data[t_data::V_RADIAL](n_radial,n_azimuthal) = 0.0;
				data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) = sqrt(constants::G*M/r);
			}
		}

		return;
	}

	/* Pressure is already initialized: cf initeuler in SourceEuler.c ... */
	/* --------- */

	// Initialization of azimutal velocity with exact centrifugal balance
	/* --------- */
	if (CentrifugalBalance) {
		double vt_int[MAX1D];

		/* vt_int \equiv rOmega� = grad(P)/sigma +  \partial(phi)/\partial(r)  -  acc_sg_radial */
		mpi_make1Dprofile(data[t_data::PRESSURE].Field, GLOBAL_bufarray);

		// get global SigmaMed
		t_radialarray SigmaMedGlobal(GlobalNRadial);
		MPI_Status status;
		MPI_Request request;

		for (int i = 0; i < CPU_Number; ++i) {
			if (i == CPU_Rank) {
				// receive data from prev CPU
				if (CPU_Rank > 0) {
					MPI_Irecv(&SigmaMedGlobal[0], IMIN+CPUOVERLAP-1, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &request);
					MPI_Wait(&request, &status);
				}
				// insert own data
				for (unsigned int j = IMIN+(CPU_Rank == 0 ? 0 : CPUOVERLAP); j <= IMAX-(CPU_Rank == CPU_Highest ? 0 : CPUOVERLAP); ++j) {
					SigmaMedGlobal[j] = SigmaMed[Zero_or_active+j-IMIN-(CPU_Rank == 0 ? 0 : CPUOVERLAP)];
				}
				// send data to next CPU
				if (CPU_Rank < CPU_Highest) {
					MPI_Isend(&SigmaMedGlobal[0], IMAX-CPUOVERLAP, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &request);
					MPI_Wait(&request, &status);
				}
			}
		}

		// now CPU_Highst has all data. broadcast to everybody
		MPI_Bcast(&SigmaMedGlobal[0], GlobalNRadial, MPI_DOUBLE, CPU_Highest, MPI_COMM_WORLD);


		/* global axisymmetric pressure field, known by all cpus */
		for (unsigned int i = 1; i < GlobalNRadial; i++) {
			vt_int[i] = ( GLOBAL_bufarray[i] - GLOBAL_bufarray[i-1] ) / (.5*(SigmaMedGlobal[i]+SigmaMedGlobal[i-1]))/(GlobalRmed[i]-GlobalRmed[i-1]) + constants::G*(1.0/GlobalRmed[i-1]-1.0/GlobalRmed[i])/(GlobalRmed[i]-GlobalRmed[i-1]);
		}

		/* Case of a disk with self-gravity */
		if (parameters::self_gravity) { // Better test with CL rigid!
			double *GLOBAL_AxiSGAccr = (double*)malloc(sizeof(double) * GlobalNRadial);
			mpi_make1Dprofile(selfgravity::g_radial, GLOBAL_AxiSGAccr);

			for (unsigned int i = 1; i < GlobalNRadial; i++) {
				vt_int[i] -= ( (Radii[i] - GlobalRmed[i-1])*GLOBAL_AxiSGAccr[i] + (GlobalRmed[i] - Radii[i])*GLOBAL_AxiSGAccr[i-1] ) / (GlobalRmed[i]-GlobalRmed[i-1]);
			}
			free(GLOBAL_AxiSGAccr);
		}

		for (unsigned int i = 1; i < GlobalNRadial; i++)
			vt_int[i] = sqrt(vt_int[i]*Radii[i])-Radii[i]*OmegaFrame;

		t1 = vt_cent[0] = vt_int[1]+.75*(vt_int[1]-vt_int[2]);
		r1 = ConstructSequence (vt_cent, vt_int, GlobalNRadial);
		vt_cent[0] += .25*(vt_int[1]-vt_int[2]);
		t2 = vt_cent[0];
		r2 = ConstructSequence (vt_cent, vt_int, GlobalNRadial);
		t1 = t1-r1/(r2-r1)*(t2-t1);
		vt_cent[0] = t1;
		ConstructSequence(vt_cent, vt_int, GlobalNRadial);
		vt_cent[GlobalNRadial] = vt_cent[GlobalNRadial-1];
	}

	// Initialization with self-gravity, without exact centrifugal balance
	if (parameters::self_gravity && !CentrifugalBalance)
		selfgravity::init_azimuthal_velocity(data[t_data::V_AZIMUTHAL]);

	for (unsigned int n_radial = 0; n_radial <= data[t_data::V_AZIMUTHAL].get_max_radial(); ++n_radial) {
		if (n_radial == data[t_data::V_AZIMUTHAL].Nrad) {
			r = Rmed[data[t_data::V_AZIMUTHAL].Nrad-1];
			ri= Rinf[data[t_data::V_AZIMUTHAL].Nrad-1];
		} else {
			r = Rmed[n_radial];
			ri= Rinf[n_radial];
		}

		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::V_AZIMUTHAL].get_max_azimuthal(); ++n_azimuthal) {
			if (!parameters::self_gravity) {
				// v_azimuthal = Omega_K * r * (...)
				data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) = r * omega_kepler(r) * sqrt(1.0 - pow2(ASPECTRATIO_REF) * pow(r,2.0*FLARINGINDEX) * (1.+SIGMASLOPE-2.0*FLARINGINDEX) );
			}

			data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) -= OmegaFrame*r;

			if (CentrifugalBalance)
				data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) = vt_cent[n_radial+IMIN];

			if (n_radial == data[t_data::V_RADIAL].Nrad) {
				data[t_data::V_RADIAL](n_radial, n_azimuthal) = 0.0;
			} else {
				data[t_data::V_RADIAL](n_radial, n_azimuthal) = IMPOSEDDISKDRIFT*parameters::sigma0/SigmaInf[n_radial]/ri;

				if (!parameters::initialize_vradial_zero) {
					if (ViscosityAlpha) {
                        data[t_data::V_RADIAL](n_radial, n_azimuthal) -= 3.0*data[t_data::VISCOSITY](n_radial, n_azimuthal)/r*(-SIGMASLOPE + 2.0*FLARINGINDEX + 1.0);
					} else {
						data[t_data::V_RADIAL](n_radial, n_azimuthal) -= 3.0*data[t_data::VISCOSITY](n_radial, n_azimuthal)/r*(-SIGMASLOPE+.5);
					}
				}
			}
		}
	}

	/* set VRadial for innermost and outermost ring to 0 */
	for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::V_RADIAL].get_max_azimuthal(); ++n_azimuthal) {
		data[t_data::V_RADIAL](0,n_azimuthal) = 0.0;
		data[t_data::V_RADIAL](data[t_data::V_RADIAL].Nrad,n_azimuthal) = 0.0;
	}

}
