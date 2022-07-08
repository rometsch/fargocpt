#include <assert.h>
#include <cstring>
#include <float.h>
#include <math.h>
#include <stdio.h>

#include "LowTasks.h"
#include "Pframeforce.h"
#include "SideEuler.h"
#include "SourceEuler.h"
#include "Theo.h"
#include "TransportEuler.h"
#include "axilib.h"
#include "boundary_conditions.h"
#include "find_cell_id.h"
#include "global.h"
#include "init.h"
#include "logging.h"
#include "nongnu.h"
#include "output.h"
#include "parameters.h"
#include "pvte_law.h"
#include "quantities.h"
#include "selfgravity.h"
#include "util.h"
#include "viscosity.h"
#include <gsl/gsl_sf_bessel.h>

#include "open-simplex-noise.h"
#include "options.h"

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
    Rinf.resize(size + 1);
    InvRinf.resize(size + 1);
    Rsup.resize(size);
    Surf.resize(size);
    InvSurf.resize(size);
    InvDiffRmed.resize(size);
    InvDiffRsup.resize(size);
    InvDiffRsupRb.resize(size);
    TwoDiffRaSq.resize(size);
    TwoDiffRbSq.resize(size);
    FourThirdInvRbInvdphiSq.resize(size);
    Radii.resize(size);
    GlobalRmed.resize(size);
    SigmaInf.resize(size);
    GLOBAL_bufarray.resize(size);
    GLOBAL_AxiSGAccr.resize(size);
}

/**
	Fills 1D arrays like Rmed, Rinf, Rsup, Surf, ...
*/
void init_radialarrays()
{
    FILE *fd_input;
    char *fd_input_filename;
    unsigned int nRadial;

    if (asprintf(&fd_input_filename, "%s%s", OUTPUTDIR, "radii.dat") == -1) {
	logging::print_master(LOG_ERROR
			      "Not enough memory for string buffer.\n");
	PersonalExit(1);
    }

    fd_input = fopen(fd_input_filename, "r");

    // free up fdInputFilename
    free(fd_input_filename);

    double first_cell_size = 0.0;
    double cell_growth_factor = 0.0;
    if (fd_input == NULL) {
	logging::print_master(
	    LOG_INFO "Warning : no `radii.dat' file found. Using default.\n");
	switch (parameters::radial_grid_type) {
	case parameters::logarithmic_spacing: {
	    cell_growth_factor =
		std::pow((RMAX / RMIN), 1.0 / ((double)GlobalNRadial - 2.0));
	    for (nRadial = 0; nRadial <= GlobalNRadial + 1; ++nRadial) {

		Radii[nRadial] =
		    RMIN * std::pow(cell_growth_factor, (double)nRadial - 1.0);
	    }
	    break;
	}
	case parameters::arithmetic_spacing: {
	    cell_growth_factor = ((double)GlobalNRadial - 2.0) / (RMAX - RMIN);
	    const double interval =
		(RMAX - RMIN) / (double)(GlobalNRadial - 2.0);
	    for (nRadial = 0; nRadial <= GlobalNRadial + 1; ++nRadial) {
		Radii[nRadial] = RMIN + interval * (double)(nRadial - 1.0);
	    }
	    break;
	}
	case parameters::exponential_spacing: {
	    // Don't touch the rest
	    cell_growth_factor =
		std::pow((RMAX / RMIN), 1.0 / ((double)GlobalNRadial - 2.0));

	    first_cell_size = RMIN * (cell_growth_factor - 1.0) *
			      parameters::exponential_cell_size_factor;
	    const double f = (RMAX - RMIN) / first_cell_size;
	    double exp_growth_factor = 1.02;
	    const double Nr = (double)GlobalNRadial - 2.0;
	    for (int i = 0; i < 500000; ++i) {
		exp_growth_factor =
		    exp_growth_factor -
		    ((std::pow(exp_growth_factor, Nr) - exp_growth_factor * f +
		      f - 1)) /
			(Nr * std::pow(exp_growth_factor, Nr - 1.0) - f);
	    }
	    cell_growth_factor = exp_growth_factor;
	    for (nRadial = 0; nRadial <= GlobalNRadial + 1; ++nRadial) {
		Radii[nRadial] = RMIN + first_cell_size *
					    (std::pow(exp_growth_factor,
						      (double)nRadial - 1.0) -
					     1.0) /
					    (exp_growth_factor - 1.0);
	    }
	    break;
	}
	case parameters::custom_spacing:
	    break;
	default:
	    die("Invalid setting for RadialSpacing");
	}
    } else {
	logging::print_master(LOG_INFO "Reading 'radii.dat' file.\n");
	parameters::radial_grid_type = parameters::custom_spacing;
	for (nRadial = 0; nRadial <= GlobalNRadial; ++nRadial) {
	    double temp;

	    if (fscanf(fd_input, "%lf", &temp) == 1) {
		Radii[nRadial] = temp;
	    } else {
		logging::print_master(
		    LOG_ERROR
		    "Reading 'radii.dat' file: No data left to read :(\n");
		PersonalExit(1);
	    }
	}
    }

    init_cell_finder(cell_growth_factor, first_cell_size);

    /* if input file is open, close it */
    if (fd_input != NULL) {
	fclose(fd_input);
    }

    for (nRadial = 0; nRadial < GlobalNRadial + 1; ++nRadial) {
	// Rmed is in the center of the cell where the center of mass is
	// Rmed = 1/2 * [ (4/3 Pi r_sup^3) - (4/3 Pi r_inf^3) ] / [ (Pi r_sup^2)
	// - (Pi r_inf^2) ]
	GlobalRmed[nRadial] =
	    2.0 / 3.0 *
	    (Radii[nRadial + 1] * Radii[nRadial + 1] * Radii[nRadial + 1] -
	     Radii[nRadial] * Radii[nRadial] * Radii[nRadial]);
	GlobalRmed[nRadial] =
	    GlobalRmed[nRadial] / (Radii[nRadial + 1] * Radii[nRadial + 1] -
				   Radii[nRadial] * Radii[nRadial]);
    }

    logging::print_master(
	LOG_VERBOSE
	"Active %s grid is ranging from %g to %g. Total grid is range from %g to %g.\n",
	parameters::radial_grid_names[parameters::radial_grid_type], Radii[1],
	Radii[GlobalNRadial - 1], Radii[0], Radii[GlobalNRadial]);

    for (nRadial = 0; nRadial < NRadial + 1; ++nRadial) {
	Rinf[nRadial] = Radii[nRadial + IMIN];
	Rsup[nRadial] = Radii[nRadial + IMIN + 1];

	Rmed[nRadial] = 2.0 / 3.0 *
			(Rsup[nRadial] * Rsup[nRadial] * Rsup[nRadial] -
			 Rinf[nRadial] * Rinf[nRadial] * Rinf[nRadial]);
	Rmed[nRadial] = Rmed[nRadial] / (Rsup[nRadial] * Rsup[nRadial] -
					 Rinf[nRadial] * Rinf[nRadial]);

	// TODO: Is already calculated a few lines above. assert should check
	// this
	assert((Rmed[nRadial] - GlobalRmed[nRadial + IMIN]) < DBL_EPSILON);

	Surf[nRadial] =
	    M_PI * (std::pow(Rsup[nRadial], 2) - std::pow(Rinf[nRadial], 2)) /
	    (double)NAzimuthal;

	InvRmed[nRadial] = 1.0 / Rmed[nRadial];
	InvSurf[nRadial] = 1.0 / Surf[nRadial];
	InvDiffRsup[nRadial] = 1.0 / (Rsup[nRadial] - Rinf[nRadial]);
	InvDiffRsupRb[nRadial] =
	    1.0 / ((Rsup[nRadial] - Rinf[nRadial]) * Rmed[nRadial]);
	TwoDiffRaSq[nRadial] = 2.0 / (Rsup[nRadial] * Rsup[nRadial] -
				      Rinf[nRadial] * Rinf[nRadial]);
	FourThirdInvRbInvdphiSq[nRadial] =
	    4.0 / 3.0 / Rmed[nRadial] * invdphi * invdphi;
	InvRinf[nRadial] = 1.0 / Rinf[nRadial];
    }

    Rinf[NRadial] = Radii[NRadial + IMIN];
    InvRinf[NRadial] = 1.0 / Rinf[NRadial];

    for (nRadial = 1; nRadial < NRadial + 1; ++nRadial) {
	InvDiffRmed[nRadial] = 1.0 / (Rmed[nRadial] - Rmed[nRadial - 1]);
	TwoDiffRbSq[nRadial] = 2.0 / (Rmed[nRadial] * Rmed[nRadial] -
				      Rmed[nRadial - 1] * Rmed[nRadial - 1]);
    }

    /* output radii to used_rad.dat (on master only) */
    if (CPU_Master) {
	FILE *fd_output;
	char *fd_output_filename;

	if (asprintf(&fd_output_filename, "%s%s", OUTPUTDIR, "used_rad.dat") ==
	    -1) {
	    logging::print_master(LOG_ERROR
				  "Not enough memory for string buffer.\n");
	    PersonalExit(1);
	}

	fd_output = fopen(fd_output_filename, "w");

	if (fd_output == NULL) {
	    logging::print_master(LOG_ERROR
				  "Can't write %s.\nProgram stopped.\n",
				  fd_output_filename);
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

void init_physics(t_data &data)
{
    if ((parameters::sigma_initialize_condition ==
	 parameters::initialize_condition_shakura_sunyaev) &&
	(parameters::energy_initialize_condition ==
	 parameters::initialize_condition_shakura_sunyaev)) {
	init_shakura_sunyaev(data);
	return;
    } else if ((parameters::sigma_initialize_condition ==
		parameters::initialize_condition_shakura_sunyaev) ||
	       (parameters::energy_initialize_condition ==
		parameters::initialize_condition_shakura_sunyaev)) {
	die("Both Sigma and Energy have to be initialised by Shakura & Sunyaev Standard-Solution. Other initialisation not yet implemented!");
    }

    if (parameters::ShockTube == 1) {
	init_shock_tube_test(data);
    } else if (parameters::ShockTube == 2) {
	init_PVTE_shock_tube_test(data);
	init_eos_arrays(data);
    } else {
	// gas density initialization
	init_gas_density(data);

	// if energy equation is taken into account, we initialize the gas
	// thermal energy
	if (parameters::Adiabatic) {
	    init_gas_energy(data);
	}

	if (parameters::do_init_secondary_disk) {
	    init_secondary_disk_densities(data);
	    if (parameters::Adiabatic) {
		init_secondary_disk_energies(data);
	    }
	}

	if (parameters::variableGamma) {
	    init_eos_arrays(data);
	}

	renormalize_sigma_and_report(data);
    }

    if (parameters::self_gravity) {
	// if SelfGravity = YES or Z, planets are initialized feeling disk
	// potential. Only the surface density is required to calculate the
	// radial self-gravity acceleration. The disk radial and azimutal
	// velocities are not updated
	selfgravity::init(data);
	logging::print_master(LOG_INFO "sg initialised\n");
    }

    // ListPlanets(sys);
    data.get_planetary_system().list_planets();
    // data.get_planetary_system().correct_planet_accretion();

	if (data.get_planetary_system().get_number_of_planets() < 2) {
	/*if (parameters::boundary_outer ==
	    parameters::boundary_condition_center_of_mass_initial) {
	    die("Do not use 'Nbody center of mass outer boundary' with only one body!\n");
	}*/

	if (ASPECTRATIO_MODE > 0) {
	    die("Do not use Nbody aspectratio mode with only 1 body!\n");
	}
    }

    OmegaFrame = OMEGAFRAME;

    if (Corotating) {
	OmegaFrame = data.get_planetary_system()
			 .get_planet(parameters::corotation_reference_body)
			 .get_omega();
    }

    FrameAngle = 0;

    // only gas velocities remain to be initialized
    init_euler(data);
    init_gas_velocities(data);
    if (parameters::do_init_secondary_disk) {
	init_secondary_disk_velocities(data);
    }

    // boundary_conditions::apply_boundary_condition(data, 0.0, false);
}

/**
	Wrapper for initialisation of physics according to Shakura & Sunyaev
   1974
*/

void init_shakura_sunyaev(t_data &data)
{
    double factor;

    if (!parameters::Adiabatic) {
	die("Isothermal equation of state and Shakura & Sunyaev starting conditions has not yet been implemented!");
    }

    if (ASPECTRATIO_MODE > 0) {
	die("ASPECTRATIO_NBODY and Shakura & Sunyaev starting conditions has not yet been implemented!");
    }

    for (unsigned int n_radial = 0; n_radial < data[t_data::TEMPERATURE].Nrad;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::TEMPERATURE].Nsec; ++n_azimuthal) {

	    factor =
		std::pow(1. - std::sqrt(parameters::star_radius /
					(Rb[n_radial] +
					 2. * (RMAX - RMIN) / GlobalNRadial)),
			 0.25);

	    data[t_data::DENSITY](n_radial, n_azimuthal) =
		(5.2 * std::pow(ALPHAVISCOSITY, -4. / 5.) *
		 std::pow(parameters::mass_accretion_rate *
			      units::mass_accretion_rate.get_cgs_factor() /
			      1.e16,
			  7. / 10.) *
		 std::pow(parameters::M0, 0.25) *
		 std::pow(Rb[n_radial] * units::length.get_cgs_factor() / 1.e10,
			  -0.75) *
		 std::pow(factor, 14. / 5.)) /
		units::surface_density.get_cgs_factor();
	    data[t_data::SCALE_HEIGHT](n_radial, n_azimuthal) =
		(1.7e8 * std::pow(ALPHAVISCOSITY, -1. / 10.) *
		 std::pow(parameters::mass_accretion_rate *
			      units::mass_accretion_rate.get_cgs_factor() /
			      1.e16,
			  3. / 20.) *
		 std::pow(parameters::M0, -3. / 8.) *
		 std::pow(Rb[n_radial] * units::length.get_cgs_factor() / 1.e10,
			  9. / 8.) *
		 std::pow(factor, 3. / 5.)) /
		(Rb[n_radial] * units::length.get_cgs_factor()) * Rb[n_radial];
	    data[t_data::TEMPERATURE](n_radial, n_azimuthal) =
		(1.4e4 * std::pow(ALPHAVISCOSITY, -1. / 5.) *
		 std::pow(parameters::mass_accretion_rate *
			      units::mass_accretion_rate.get_cgs_factor() /
			      1.e16,
			  3. / 10.) *
		 std::pow(parameters::M0, 0.25) *
		 std::pow(Rb[n_radial] * units::length.get_cgs_factor() / 1.e10,
			  -0.75) *
		 std::pow(factor, 6. / 5.)) /
		units::temperature.get_cgs_factor();
	    data[t_data::V_RADIAL](n_radial, n_azimuthal) =
		-(2.7e4 * std::pow(ALPHAVISCOSITY, 4. / 5.) *
		  std::pow(parameters::mass_accretion_rate *
			       units::mass_accretion_rate.get_cgs_factor() /
			       1.e16,
			   3. / 10.) *
		  std::pow(parameters::M0, -0.25) *
		  std::pow(Rb[n_radial] * units::length.get_cgs_factor() /
			       1.e10,
			   -0.25) *
		  std::pow(factor, -14. / 5.)) /
		units::velocity.get_cgs_factor();

	    data[t_data::SOUNDSPEED](n_radial, n_azimuthal) =
		std::sqrt(constants::R / parameters::MU * ADIABATICINDEX *
			  data[t_data::TEMPERATURE](n_radial, n_azimuthal));
	    data[t_data::ENERGY](n_radial, n_azimuthal) =
		constants::R / parameters::MU * 1. / (ADIABATICINDEX - 1.) *
		data[t_data::DENSITY](n_radial, n_azimuthal) *
		data[t_data::TEMPERATURE](n_radial, n_azimuthal);
	    data[t_data::PRESSURE](n_radial, n_azimuthal) =
		(ADIABATICINDEX - 1.) *
		data[t_data::ENERGY](n_radial, n_azimuthal);
	    data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) =
		calculate_omega_kepler(Rb[n_radial]) * Rb[n_radial];
	}
    }

    RefillSigma(&data[t_data::DENSITY]);
    RefillEnergy(&data[t_data::ENERGY]);

    if (parameters::self_gravity) {
	die("Self-gravity and Shakura-Sunyaev starting values has not yet been implemented!");
    }

    InitCellCenterCoordinates();

    /** init_euler w/o updates that have already been done above **/
    InitTransport();

    viscosity::update_viscosity(data);
    /** end init_euler **/

    if (CentrifugalBalance)
	die("CentrifugalBalance and Shakura-Sunyaev starting values has not yet been implemented!");
}

/**
	Initializes density and energy for the spreading ring test.
	Intented to be used with spreading_ring.par
	See R. Speith and W. Kley 2003: Stability of the viscously spreading
   ring
*/
void init_spreading_ring_test(t_data &data)
{

    const double R0 = 1.0;
    double R0_ = R0;

    unsigned int R0_id = 0;
    for (unsigned int n_radial = 0; n_radial < data[t_data::DENSITY].Nrad;
	 ++n_radial) {
	if (Rsup[n_radial] > R0 && R0 > Rinf[n_radial]) {
	    R0_ = Rmed[n_radial];
	    R0_id = n_radial;
	}
    }

    if (R0_id != 0)
	logging::print(
	    LOG_INFO "Initializing Spreading Ring at radius = %.5e\n", R0_);

    const double Disk_Mass = parameters::sigma_discmass;
    const double tau0 = 0.016;

    const double x = Rmed[R0_id] / R0;
    const double I = gsl_sf_bessel_Inu(0.25, 2.0 * x / tau0);
    const double Sigma0 = Disk_Mass / (M_PI * R0 * R0) * 1.0 /
			  (tau0 * std::pow(x, 0.25)) * I *
			  std::exp(-(1.0 + x * x) / tau0);

    for (unsigned int n_radial = 0; n_radial < data[t_data::DENSITY].Nrad;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::DENSITY].Nsec; ++n_azimuthal) {
	    const double density_floor = Sigma0 * parameters::sigma_floor;
	    const double energy = 0.0;

	    const double x = Rmed[n_radial] / R0;
	    const double I = gsl_sf_bessel_Inu(0.25, 2.0 * x / tau0);
	    double density = Disk_Mass / (M_PI * R0 * R0) * 1.0 /
			     (tau0 * std::pow(x, 0.25)) * I *
			     std::exp(-(1.0 + x * x) / tau0);

	    density = std::max(density, density_floor);

	    data[t_data::DENSITY](n_radial, n_azimuthal) = density;
	    data[t_data::ENERGY](n_radial, n_azimuthal) = energy;
	}
    }

    // set SigmaMed/SigmaInf
    RefillSigma(&data[t_data::DENSITY]);
    RefillEnergy(&data[t_data::ENERGY]);
}

/**
	Initializes density and energy for the spreading ring test.
	Intented to be used with spreading_ring.par
	See R. Speith and W. Kley 2003: Stability of the viscously spreading
   ring
   Init function copied to be equal to Jibin's pluto simulations
*/
void init_spreading_ring_test_jibin(t_data &data)
{

    const double Disk_Mass = parameters::sigma_discmass;
    const double Ring_Mass = 1.0e-4;
    const double h = ASPECTRATIO_REF;
    const double p = SIGMASLOPE;
    const double q = 2.0 * FLARINGINDEX - 1.0;
    const double tau0 = 0.018;
    const double Rmin = RMIN;
    const double Rmax = RMAX;

    const double sig0 = Disk_Mass / (2.0 * M_PI) * (p + 2.0) /
			(std::pow(Rmax, p + 2.0) - std::pow(Rmin, p + 2.0));
    parameters::sigma0 = sig0;

    logging::print_master(LOG_INFO "Initializing viscous spreading ring\n");

    logging::print_master(
	LOG_INFO "spreading ring sig0code = %.5e	sig0cgs = %.5e\n", sig0,
	sig0 * units::surface_density.get_cgs_factor());

    for (unsigned int n_radial = 0; n_radial < data[t_data::DENSITY].Nrad;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::DENSITY].Nsec; ++n_azimuthal) {

	    const double R = Rmed[n_radial];
	    const double OmegaK = 1.0 / (R * std::sqrt(R));
	    const double arg = 2.0 * R / tau0;
	    const double bessel = gsl_sf_bessel_Inu(0.25, arg);

	    const double sig_ring = Ring_Mass * bessel *
				    std::exp((-1.0 - R * R) / tau0) /
				    (M_PI * tau0 * std::pow(R, 0.25));
	    const double sig_disk = sig0 * std::pow(R, p);
	    const double energy = 0.0;

	    const double sig_noise =
		0.05 * sig_disk * (1.0 - 2.0 * (rand() / (double)RAND_MAX));
	    const double vr = 0.0;
	    const double corr = std::sqrt(1.0 + (p + q) * h * h);
	    const double vaz = R * OmegaK * corr - R * OmegaFrame;

	    const double sig = sig_ring + sig_disk + sig_noise;

	    data[t_data::DENSITY](n_radial, n_azimuthal) = sig;
	    data[t_data::ENERGY](n_radial, n_azimuthal) = energy;
	    data[t_data::V_RADIAL](n_radial, n_azimuthal) = vr;
	    data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) = vaz;
	}
    }

    for (unsigned int n_azimuthal = 0;
	 n_azimuthal < data[t_data::V_RADIAL].Nsec; ++n_azimuthal) {
	data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial(),
			       n_azimuthal) = 0.0;
    }

    // set SigmaMed/SigmaInf
    RefillSigma(&data[t_data::DENSITY]);
    RefillEnergy(&data[t_data::ENERGY]);
}

/**
	Initializes density and energy for Shock Tube test.
	Intented to be used with shock_tube.par
*/
void init_shock_tube_test(t_data &data)
{
    logging::print_master(LOG_INFO "Initializing ShockTube\n");

    for (unsigned int n_radial = 0; n_radial < data[t_data::DENSITY].Nrad;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::DENSITY].Nsec; ++n_azimuthal) {
	    double density = 1.0;
	    double energy = 2.5;

	    if (Rmed[n_radial] - GlobalRmed[0] > 0.5) {
		density = 0.125;
		energy = 2.0 * 0.125;
	    }

	    data[t_data::DENSITY](n_radial, n_azimuthal) = density;
	    data[t_data::ENERGY](n_radial, n_azimuthal) = energy;
	}
    }

    units::length.set_cgs_factor(1.0);
    units::length.set_cgs_symbol("1");

    units::mass.set_cgs_factor(1.0);
    units::mass.set_cgs_symbol("1");

    units::time.set_cgs_factor(1.0);
    units::time.set_cgs_symbol("1");

    units::energy.set_cgs_factor(1.0);
    units::energy.set_cgs_symbol("1");

    units::energy_density.set_cgs_factor(1.0);
    units::energy_density.set_cgs_symbol("1");

    units::temperature.set_cgs_factor(1.0);
    units::temperature.set_cgs_symbol("1");

    units::density.set_cgs_factor(1.0);
    units::density.set_cgs_symbol("1");

    units::surface_density.set_cgs_factor(1.0);
    units::surface_density.set_cgs_symbol("1");

    units::opacity.set_cgs_factor(1.0);
    units::opacity.set_cgs_symbol("1");

    units::energy_flux.set_cgs_factor(1.0);
    units::energy_flux.set_cgs_symbol("1");

    units::velocity.set_cgs_factor(1.0);
    units::velocity.set_cgs_symbol("1");

    units::acceleration.set_cgs_factor(1.0);
    units::acceleration.set_cgs_symbol("1");

    units::angular_momentum.set_cgs_factor(1.0);
    units::angular_momentum.set_cgs_symbol("1");

    units::kinematic_viscosity.set_cgs_factor(1.0);
    units::kinematic_viscosity.set_cgs_symbol("1");

    units::dynamic_viscosity.set_cgs_factor(1.0);
    units::dynamic_viscosity.set_cgs_symbol("1");

    units::stress.set_cgs_factor(1.0);
    units::stress.set_cgs_symbol("1");

    units::pressure.set_cgs_factor(1.0);
    units::pressure.set_cgs_symbol("1");

    units::power.set_cgs_factor(1.0);
    units::power.set_cgs_symbol("1");

    units::potential.set_cgs_factor(1.0);
    units::potential.set_cgs_symbol("1");

    units::torque.set_cgs_factor(1.0);
    units::torque.set_cgs_symbol("1");

    units::force.set_cgs_factor(1.0);
    units::force.set_cgs_symbol("1");

    units::mass_accretion_rate.set_cgs_factor(1.0);
    units::mass_accretion_rate.set_cgs_symbol("1");

    // after all units have calculated, calculate constants in code units
    constants::calculate_constants_in_code_units();

    // set SigmaMed/SigmaInf
    RefillSigma(&data[t_data::DENSITY]);
    RefillEnergy(&data[t_data::ENERGY]);
}

void init_PVTE_shock_tube_test(t_data &data)
{
    logging::print_master(LOG_INFO "Initializing PVTE-ShockTube\n");

    for (unsigned int n_radial = 0; n_radial < data[t_data::DENSITY].Nrad;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::DENSITY].Nsec; ++n_azimuthal) {
	    double density = 1.0;
	    double energy = 10.361627466581034;

	    if (Rmed[n_radial] - GlobalRmed[0] > 0.5) {
		density = 0.125;
		energy = 0.9110851732216827;
	    }

	    data[t_data::DENSITY](n_radial, n_azimuthal) = density;
	    data[t_data::ENERGY](n_radial, n_azimuthal) = energy;
	}
    }

    units::length.set_cgs_factor(1.0);
    units::length.set_cgs_symbol("1");

    units::mass.set_cgs_factor(1.0);
    units::mass.set_cgs_symbol("1");

    units::time.set_cgs_factor(1.0);
    units::time.set_cgs_symbol("1");

    units::energy.set_cgs_factor(1.0);
    units::energy.set_cgs_symbol("1");

    units::energy_density.set_cgs_factor(0.00000004576860232875);
    units::energy_density.set_cgs_symbol("1");

    units::temperature.set_cgs_factor(3341.5268389972975);
    units::temperature.set_cgs_symbol("1");

    units::density.set_cgs_factor(1.66053886e-19);
    units::density.set_cgs_symbol("1");

    units::surface_density.set_cgs_factor(1.0);
    units::surface_density.set_cgs_symbol("1");

    units::opacity.set_cgs_factor(1.0);
    units::opacity.set_cgs_symbol("1");

    units::energy_flux.set_cgs_factor(1.0);
    units::energy_flux.set_cgs_symbol("1");

    units::velocity.set_cgs_factor(5.25e5);
    units::velocity.set_cgs_symbol("1");

    units::acceleration.set_cgs_factor(1.0);
    units::acceleration.set_cgs_symbol("1");

    units::angular_momentum.set_cgs_factor(1.0);
    units::angular_momentum.set_cgs_symbol("1");

    units::kinematic_viscosity.set_cgs_factor(1.0);
    units::kinematic_viscosity.set_cgs_symbol("1");

    units::dynamic_viscosity.set_cgs_factor(1.0);
    units::dynamic_viscosity.set_cgs_symbol("1");

    units::stress.set_cgs_factor(1.0);
    units::stress.set_cgs_symbol("1");

    units::pressure.set_cgs_factor(1.0);
    units::pressure.set_cgs_symbol("1");

    units::power.set_cgs_factor(1.0);
    units::power.set_cgs_symbol("1");

    units::torque.set_cgs_factor(1.0);
    units::torque.set_cgs_symbol("1");

    units::force.set_cgs_factor(1.0);
    units::force.set_cgs_symbol("1");

    units::mass_accretion_rate.set_cgs_factor(1.0);
    units::mass_accretion_rate.set_cgs_symbol("1");

    // after all units have calculated, calculate constants in code units
    constants::calculate_constants_in_code_units();

    compute_pressure(data, false);

    // set SigmaMed/SigmaInf
    RefillSigma(&data[t_data::DENSITY]);
    RefillEnergy(&data[t_data::ENERGY]);
}

void init_secondary_disk_densities(t_data &data)
{

    logging::print_master(LOG_INFO "Initializing Secondary disk densities\n");

    if (data.get_planetary_system().get_number_of_planets() < 2) {
	die("Error: cannot initialize secondary disk with only %d nbody objects!\n",
	    data.get_planetary_system().get_number_of_planets());
    }

    const auto &planet = data.get_planetary_system().get_planet(1);
    const double disk_size = parameters::profile_cutoff_point_outer *
			     planet.get_dimensionless_roche_radius() /
			     (1.0 - planet.get_dimensionless_roche_radius());
    const double cutoff_width = parameters::profile_cutoff_width_outer *
				planet.get_dimensionless_roche_radius() /
				(1.0 - planet.get_dimensionless_roche_radius());
    const double mass_q = planet.get_mass() /
			  data.get_planetary_system().get_planet(0).get_mass();
    const double compute_radius =
	eggleton_1983(mass_q, planet.get_distance_to_primary());
    const double scaling_factor = std::sqrt(planet.get_mass());

    const double min_dist = RMIN / 3.0;

    for (unsigned int n_radial = 0; n_radial < data[t_data::DENSITY].Nrad;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::DENSITY].Nsec; ++n_azimuthal) {

	    const double phi = (double)n_azimuthal * dphi;
	    const double rmed = Rmed[n_radial];
	    const double x = rmed * std::cos(phi) - planet.get_x();
	    const double y = rmed * std::sin(phi) - planet.get_y();
	    const double r = std::max(std::sqrt(x * x + y * y), min_dist);

	    if (r < compute_radius) {
		const double density = parameters::sigma0 * scaling_factor *
				       std::pow(r, -SIGMASLOPE) *
				       cutoff_outer(disk_size, cutoff_width, r);

		const double density_old =
		    std::max(data[t_data::DENSITY](n_radial, n_azimuthal),
			     parameters::sigma_floor * parameters::sigma0);

		data[t_data::DENSITY](n_radial, n_azimuthal) =
		    std::max(density, density_old);
	    }
	}
    }
}

void init_secondary_disk_energies(t_data &data)
{

    if (!parameters::Adiabatic) {
	return;
    }
    logging::print_master(LOG_INFO "Initializing Secondary disk energies\n");

    if (data.get_planetary_system().get_number_of_planets() < 2) {
	die("Error: cannot initialize secondary disk with only %d nbody objects!\n",
	    data.get_planetary_system().get_number_of_planets());
    }

    const auto &planet = data.get_planetary_system().get_planet(1);
    const double disk_size = parameters::profile_cutoff_point_outer *
			     planet.get_dimensionless_roche_radius() /
			     (1.0 - planet.get_dimensionless_roche_radius());
    const double cutoff_width = parameters::profile_cutoff_width_outer *
				planet.get_dimensionless_roche_radius() /
				(1.0 - planet.get_dimensionless_roche_radius());
    const double mass_q = planet.get_mass() /
			  data.get_planetary_system().get_planet(0).get_mass();
    const double compute_radius =
	eggleton_1983(mass_q, planet.get_distance_to_primary());
    const double scaling_factor = std::sqrt(planet.get_mass());

    const double min_dist = RMIN / 3.0;

    for (unsigned int n_radial = 0;
	 n_radial < data[t_data::ENERGY].get_size_radial(); ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::ENERGY].get_size_azimuthal();
	     ++n_azimuthal) {

	    const double phi = (double)n_azimuthal * dphi;
	    const double rmed = Rmed[n_radial];
	    const double x = rmed * std::cos(phi) - planet.get_x();
	    const double y = rmed * std::sin(phi) - planet.get_y();
	    const double r = std::max(std::sqrt(x * x + y * y), min_dist);

	    if (r < compute_radius) {
		const double energy =
		    1.0 / (ADIABATICINDEX - 1.0) * parameters::sigma0 *
		    scaling_factor * std::pow(ASPECTRATIO_REF, 2) *
		    std::pow(r, -SIGMASLOPE - 1.0 + 2.0 * FLARINGINDEX) *
		    constants::G * planet.get_mass() *
		    cutoff_outer(disk_size, cutoff_width, r);

		const double temperature_floor =
		    parameters::minimum_temperature *
		    units::temperature.get_inverse_cgs_factor();
		const double energy_floor =
		    temperature_floor *
		    data[t_data::DENSITY](n_radial, n_azimuthal) /
		    parameters::MU * constants::R / (ADIABATICINDEX - 1.0);

		const double temperature_ceil =
		    parameters::maximum_temperature *
		    units::temperature.get_inverse_cgs_factor();
		const double energy_ceil =
		    temperature_ceil *
		    data[t_data::DENSITY](n_radial, n_azimuthal) /
		    parameters::MU * constants::R / (ADIABATICINDEX - 1.0);

		const double energy_old = std::max(
		    data[t_data::ENERGY](n_radial, n_azimuthal), energy_floor);
		data[t_data::ENERGY](n_radial, n_azimuthal) =
		    std::min(std::max(energy, energy_old), energy_ceil);
	    }
	}
    }
}

void init_secondary_disk_velocities(t_data &data)
{

    logging::print_master(LOG_INFO "Initializing Secondary disk velocities\n");

    if (CentrifugalBalance) {
	logging::print_master(
	    LOG_INFO "CentrifugalBalance not tested with secondary disk!");
    }

    if (parameters::self_gravity) {
	logging::print_master(LOG_INFO
			      "Self gravity not tested with secondary disk!");
    }

    if (data.get_planetary_system().get_number_of_planets() < 2) {
	die("Error: cannot initialize secondary disk with only %d nbody objects!\n",
	    data.get_planetary_system().get_number_of_planets());
    }

    const auto &planet = data.get_planetary_system().get_planet(1);
    const double mass_q = planet.get_mass() /
			  data.get_planetary_system().get_planet(0).get_mass();
    const double compute_radius =
	eggleton_1983(mass_q, planet.get_distance_to_primary());

    const double min_dist = RMIN / 3.0;

    for (unsigned int n_radial = 0;
	 n_radial <= data[t_data::V_RADIAL].get_max_radial(); ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::V_RADIAL].get_max_azimuthal();
	     ++n_azimuthal) {

	    const double phi = (double)n_azimuthal * dphi;
	    double rinf;
	    if (n_radial == data[t_data::V_AZIMUTHAL].Nrad) {
		rinf = Rinf[data[t_data::V_AZIMUTHAL].Nrad - 1];
	    } else {
		rinf = Rinf[n_radial];
	    }

	    const double cell_x = rinf * std::cos(phi);
	    const double cell_y = rinf * std::sin(phi);

	    // Position in secondary frame
	    const double x_sec = cell_x - planet.get_x();
	    const double y_sec = cell_y - planet.get_y();
	    const double r_sec =
		std::max(std::sqrt(x_sec * x_sec + y_sec * y_sec), min_dist);

	    if (r_sec < compute_radius) {

		// pressure support correction
		double corr;
		if (parameters::initialize_pure_keplerian) {
		    corr = 1.0;
		} else {
		    corr = std::sqrt(
			1.0 - std::pow(ASPECTRATIO_REF, 2) *
				  std::pow(r_sec, 2.0 * FLARINGINDEX) *
				  (1. + SIGMASLOPE - 2.0 * FLARINGINDEX));
		}

		// Velocities in center of mass frame
		const double vr_sec = 0.0;
		const double vaz_sec =
		    std::sqrt(constants::G * planet.get_mass() / r_sec) * corr;

		const double vx_sec =
		    (vr_sec * x_sec - vaz_sec * y_sec) / r_sec;
		const double vy_sec =
		    (vr_sec * y_sec + vaz_sec * x_sec) / r_sec;

		// shift velocity from center of mass frame to primary frame
		const double vx = vx_sec - planet.get_vx();
		const double vy = vy_sec - planet.get_vy();

		const double vr = vx * std::cos(phi) + vy * std::sin(phi);
		data[t_data::V_RADIAL](n_radial, n_azimuthal) = vr;
	    }
	}
    }

    for (unsigned int n_radial = 0;
	 n_radial < data[t_data::V_AZIMUTHAL].get_size_radial(); ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::V_AZIMUTHAL].get_size_azimuthal();
	     ++n_azimuthal) {

	    const double phi = ((double)n_azimuthal - 0.5) * dphi;
	    double rmed;
	    if (n_radial == data[t_data::V_AZIMUTHAL].Nrad) {
		rmed = Rmed[data[t_data::V_AZIMUTHAL].Nrad - 1];
	    } else {
		rmed = Rmed[n_radial];
	    }

	    const double cell_x = rmed * std::cos(phi);
	    const double cell_y = rmed * std::sin(phi);

	    // Position in secondary frame
	    const double x_sec = cell_x - planet.get_x();
	    const double y_sec = cell_y - planet.get_y();
	    const double r_sec =
		std::max(std::sqrt(x_sec * x_sec + y_sec * y_sec), min_dist);

	    if (r_sec < compute_radius) {

		// pressure support correction
		double corr;
		if (parameters::initialize_pure_keplerian) {
		    corr = 1.0;
		} else {
		    corr = std::sqrt(
			1.0 - std::pow(ASPECTRATIO_REF, 2) *
				  std::pow(r_sec, 2.0 * FLARINGINDEX) *
				  (1. + SIGMASLOPE - 2.0 * FLARINGINDEX));
		}

		// Velocities in center of mass frame
		const double vr_sec = 0.0;
		const double vaz_sec =
		    std::sqrt(constants::G * planet.get_mass() / r_sec) * corr;

		const double vx_sec =
		    (vr_sec * x_sec - vaz_sec * y_sec) / r_sec;
		const double vy_sec =
		    (vr_sec * y_sec + vaz_sec * x_sec) / r_sec;

		// shift velocities from center of mass frame to primary frame
		const double vx = vx_sec - planet.get_vx();
		const double vy = vy_sec - planet.get_vy();

		const double vaz = vy * std::cos(phi) - vx * std::sin(phi);
		data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) = vaz;
	    }
	}
    }
}

void init_gas_density(t_data &data)
{
    switch (parameters::sigma_initialize_condition) {
    case parameters::initialize_condition_profile:
	logging::print_master(
	    LOG_INFO "Initializing Sigma(r) = %g = %g %s * [r/(%g AU)]^(%g)\n",
	    parameters::sigma0,
	    parameters::sigma0 * units::surface_density.get_cgs_factor(),
	    units::surface_density.get_cgs_symbol(),
	    units::length.get_cgs_factor() / units::cgs_AU, -SIGMASLOPE);

	for (unsigned int n_radial = 0; n_radial < data[t_data::DENSITY].Nrad;
	     ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::DENSITY].Nsec; ++n_azimuthal) {
		const double density =
		    parameters::sigma0 * std::pow(Rmed[n_radial], -SIGMASLOPE);
		const double density_floor =
		    parameters::sigma_floor * parameters::sigma0;
		data[t_data::DENSITY](n_radial, n_azimuthal) =
		    std::max(density, density_floor);
	    }
	}
	break;

    case parameters::initialize_condition_profile_Nbody_centered: {
	logging::print_master(
	    LOG_INFO
	    "Initializing from CMS Sigma(r) = %g = %g %s * [r/(%g AU)]^(%g)\n",
	    parameters::sigma0,
	    parameters::sigma0 * units::surface_density.get_cgs_factor(),
	    units::surface_density.get_cgs_symbol(),
	    units::length.get_cgs_factor() / units::cgs_AU, -SIGMASLOPE);

	Pair cms = data.get_planetary_system().get_center_of_mass();
	const double cms_x = cms.x;
	const double cms_y = cms.y;

	for (unsigned int n_radial = 0;
	     n_radial < data[t_data::DENSITY].get_size_radial(); ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::DENSITY].get_size_azimuthal();
		 ++n_azimuthal) {

		const double phi = (double)n_azimuthal * dphi;
		const double rmed = Rmed[n_radial];
		const double x = rmed * std::cos(phi) - cms_x;
		const double y = rmed * std::sin(phi) - cms_y;
		const double r = std::sqrt(x * x + y * y);

		const double density =
		    parameters::sigma0 * std::pow(r, -SIGMASLOPE);
		const double density_floor =
		    parameters::sigma_floor * parameters::sigma0;
		data[t_data::DENSITY](n_radial, n_azimuthal) =
		    std::max(density, density_floor);
	    }
	}
    } break;

    case parameters::initialize_condition_read1D:
	logging::print_master(LOG_INFO "Loading Sigma from '%s' (1D).\n",
			      parameters::sigma_filename);
	data[t_data::DENSITY].read1D(parameters::sigma_filename, true);
	break;

    case parameters::initialize_condition_read2D:
	logging::print_master(LOG_INFO "Loading Sigma from '%s' (2D).\n",
			      parameters::sigma_filename);
	data[t_data::DENSITY].read2D(parameters::sigma_filename);
	break;

    case parameters::initialize_condition_shakura_sunyaev:
	die("Bad choice!"); // TODO: better explanation!
	break;
	// 		case parameters::initialize_condition_shakura_sunyaev:
	// 			logging::print_master(LOG_INFO "Initializing
	// Sigma from Shakura and Sunyaev 1973 standard solution (cf. A&A, 24,
	// 337)");
	//
	// 			for (unsigned int n_radial = 0; n_radial <
	// data[t_data::DENSITY].Nrad; ++n_radial) { for (unsigned int
	// n_azimuthal = 0; n_azimuthal < data[t_data::DENSITY].Nsec;
	// ++n_azimuthal) {
	// data[t_data::DENSITY](n_radial, n_azimuthal) =
	// 					data[t_data::DENSITY](n_radial,
	// n_azimuthal) = parameters::sigma0*pow(Rmed[n_radial],-SIGMASLOPE);
	// 				}
	// 			}
	// 			break;
    }

    if (parameters::SpreadingRing) {
	init_spreading_ring_test_jibin(data);
    }

    if (parameters::sigma_randomize) {

	struct osn_context *osn;
	if (open_simplex_noise(parameters::random_seed + 325582, &osn) !=
	    0) // add random number to seed to keep seed unique, because it is
	       // reused several times in the code.
	{
	    die("Bad open simplex noise!");
	}

	logging::print_master(LOG_INFO "Randomizing Sigma by %.2f %%.\n",
			      parameters::sigma_random_factor * 100);
	for (unsigned int n_radial = 0;
	     n_radial <= data[t_data::DENSITY].get_max_radial(); ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
		 ++n_azimuthal) {

		double r = Rmed[n_radial];
		double angle =
		    (double)n_azimuthal /
		    (double)data[t_data::V_RADIAL].get_size_azimuthal() * 2.0 *
		    M_PI;
		double x = r * std::cos(angle);
		double y = r * std::sin(angle);

		double f = parameters::sigma_feature_size;

		int simplex_order = 11;
		double noise = 0.0;
		for (int i = 0; i < simplex_order; ++i) {
		    double feature_factor = double(1 << i);
		    double weight_factor = double(1 << (simplex_order - i - 1));

		    noise += open_simplex_noise2(osn, feature_factor * x / f,
						 feature_factor * y / f) *
			     weight_factor;
		}

		noise /= double((1 << simplex_order) - 1);

		data[t_data::DENSITY](n_radial, n_azimuthal) *=
		    (1 + parameters::sigma_random_factor * noise);
	    }
	}

	open_simplex_noise_free(osn);
    }

    // profile cutoff at outer boundary?
    if (parameters::profile_cutoff_outer) {
	logging::print_master(
	    LOG_INFO "Cutoff Sigma for r > %g %s over a range from %g %s\n",
	    parameters::profile_cutoff_point_outer *
		units::length.get_cgs_factor(),
	    units::length.get_cgs_symbol(),
	    parameters::profile_cutoff_width_outer *
		units::length.get_cgs_factor(),
	    units::length.get_cgs_symbol());
	for (unsigned int n_radial = 0;
	     n_radial <= data[t_data::DENSITY].get_max_radial(); ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
		 ++n_azimuthal) {
		// Cutoff density to 0 for r > profile_cutoff_point_outer
		double r;
		if (parameters::sigma_initialize_condition ==
		    parameters::initialize_condition_profile_Nbody_centered) {
		    Pair cms = data.get_planetary_system().get_center_of_mass();
		    const double cms_x = cms.x;
		    const double cms_y = cms.y;

		    const double phi = (double)n_azimuthal * dphi;
		    const double rmed = Rmed[n_radial];
		    const double x = rmed * std::cos(phi) - cms_x;
		    const double y = rmed * std::sin(phi) - cms_y;
		    r = std::sqrt(x * x + y * y);
		} else {
		    r = Rmed[n_radial];
		}
		const double density_damped =
		    data[t_data::DENSITY](n_radial, n_azimuthal) *
		    cutoff_outer(parameters::profile_cutoff_point_outer,
				 parameters::profile_cutoff_width_outer, r);
		const double density_floor =
		    parameters::sigma_floor * parameters::sigma0;
		data[t_data::DENSITY](n_radial, n_azimuthal) =
		    std::max(density_damped, density_floor);
	    }
	}
    }

    // profile cutoff at inner boundary?
    if (parameters::profile_cutoff_inner) {
	logging::print_master(
	    LOG_INFO "Cutoff Sigma for r < %g %s over a range from %g %s\n",
	    parameters::profile_cutoff_point_inner *
		units::length.get_cgs_factor(),
	    units::length.get_cgs_symbol(),
	    parameters::profile_cutoff_width_inner *
		units::length.get_cgs_factor(),
	    units::length.get_cgs_symbol());
	for (unsigned int n_radial = 0;
	     n_radial <= data[t_data::DENSITY].get_max_radial(); ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
		 ++n_azimuthal) {
		// Cutoff density to 0 for r < profile_cutoff_point_inner
		double r;
		if (parameters::sigma_initialize_condition ==
		    parameters::initialize_condition_profile_Nbody_centered) {
		    Pair cms = data.get_planetary_system().get_center_of_mass();
		    const double cms_x = cms.x;
		    const double cms_y = cms.y;

		    const double phi = (double)n_azimuthal * dphi;
		    const double rmed = Rmed[n_radial];
		    const double x = rmed * std::cos(phi) - cms_x;
		    const double y = rmed * std::sin(phi) - cms_y;
		    r = std::sqrt(x * x + y * y);
		} else {
		    r = Rmed[n_radial];
		}
		const double density_damped =
		    data[t_data::DENSITY](n_radial, n_azimuthal) *
		    cutoff_inner(parameters::profile_cutoff_point_inner,
				 parameters::profile_cutoff_width_inner, r);
		const double density_floor =
		    parameters::sigma_floor * parameters::sigma0;
		data[t_data::DENSITY](n_radial, n_azimuthal) =
		    std::max(density_damped, density_floor);
	    }
	}
    }
}

void renormalize_sigma_and_report(t_data &data)
{
    // renormalize sigma0?
    if (parameters::sigma_adjust) {
	double total_mass = quantities::gas_total_mass(data, 2.0 * RMAX);
	parameters::sigma0 *= parameters::sigma_discmass / total_mass;
	logging::print_master(
	    LOG_INFO "Setting Sigma0=%g %s to set disc mass of %g to %g.\n",
	    parameters::sigma0 * units::surface_density.get_cgs_factor(),
	    units::surface_density.get_cgs_symbol(), total_mass,
	    parameters::sigma_discmass);

	// update density grid
	for (unsigned int n_radial = 0;
	     n_radial <= data[t_data::DENSITY].get_max_radial(); ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
		 ++n_azimuthal) {
		data[t_data::DENSITY](n_radial, n_azimuthal) *=
		    parameters::sigma_discmass / total_mass;
	    }
	}
    } else {
	double total_mass = quantities::gas_total_mass(data, 2.0 * RMAX);
	logging::print_master(LOG_INFO "Total disk is mass is %g = %g %s.\n",
			      total_mass,
			      total_mass * units::mass.get_cgs_factor(),
			      units::mass.get_cgs_symbol());
    }

    // set SigmaMed/SigmaInf
    RefillSigma(&data[t_data::DENSITY]);
}

void init_eos_arrays(t_data &data)
{

    logging::print_master(LOG_INFO "Generating lookup tables \n");
    pvte::initializeLookupTables();
    logging::print_master(LOG_INFO "Lookup tables generated \n");

    for (unsigned int n_rad = 0;
	 n_rad < data[t_data::GAMMAEFF].get_size_radial(); ++n_rad) {
	for (unsigned int n_az = 0;
	     n_az < data[t_data::GAMMAEFF].get_size_azimuthal(); ++n_az) {
	    data[t_data::GAMMAEFF](n_rad, n_az) = ADIABATICINDEX;
	    data[t_data::GAMMA1](n_rad, n_az) = ADIABATICINDEX;
	    data[t_data::MU](n_rad, n_az) = parameters::MU;
	}
    }
}

void init_gas_energy(t_data &data)
{
    if (ADIABATICINDEX == 1.0) {
	logging::print_master(
	    LOG_ERROR
	    "The adiabatic index must differ from unity to initialize the gas internal energy. I must exit.\n");
	PersonalExit(1);
    }

    switch (parameters::energy_initialize_condition) {
    case parameters::initialize_condition_profile:
	logging::print_master(
	    LOG_INFO
	    "Initializing Energy=%g %s * [r/(%.1f AU)]^(%g). Flaring index is %g. T=%g %s * [r/(%.1f AU)]^(%g).\n",
	    1.0 / ((ADIABATICINDEX - 1.0)) * parameters::sigma0 *
		std::pow(ASPECTRATIO_REF, 2) * units::energy.get_cgs_factor(),
	    units::energy.get_cgs_symbol(),
	    units::length.get_cgs_factor() / units::cgs_AU,
	    -SIGMASLOPE - 1.0 + 2.0 * FLARINGINDEX, FLARINGINDEX,
	    parameters::MU / constants::R * std::pow(ASPECTRATIO_REF, 2) *
		constants::G * hydro_center_mass *
		units::temperature.get_cgs_factor(),
	    units::temperature.get_cgs_symbol(),
	    units::length.get_cgs_factor() / units::cgs_AU,
	    -1.0 + 2.0 * FLARINGINDEX);

	for (unsigned int n_radial = 0;
	     n_radial <= data[t_data::ENERGY].get_max_radial(); ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal();
		 ++n_azimuthal) {
		const double energy =
		    1.0 / (ADIABATICINDEX - 1.0) * parameters::sigma0 *
		    std::pow(ASPECTRATIO_REF, 2) *
		    std::pow(Rmed[n_radial],
			     -SIGMASLOPE - 1.0 + 2.0 * FLARINGINDEX) *
		    constants::G * hydro_center_mass;
		const double temperature_floor =
		    parameters::minimum_temperature *
		    units::temperature.get_inverse_cgs_factor();
		const double energy_floor =
		    temperature_floor *
		    data[t_data::DENSITY](n_radial, n_azimuthal) /
		    parameters::MU * constants::R / (ADIABATICINDEX - 1.0);

		data[t_data::ENERGY](n_radial, n_azimuthal) =
		    std::max(energy, energy_floor);
	    }
	}
	break;

    case parameters::initialize_condition_profile_Nbody_centered: {
	const double mass = data.get_planetary_system().get_mass();
	logging::print_master(
	    LOG_INFO
	    "Initializing CMS Energy=%g %s * [r/(%.1f AU)]^(%g). Flaring index is %g. T=%g %s * [r/(%.1f AU)]^(%g).\n",
	    1.0 / ((ADIABATICINDEX - 1.0)) * parameters::sigma0 *
		std::pow(ASPECTRATIO_REF, 2) * units::energy.get_cgs_factor(),
	    units::energy.get_cgs_symbol(),
	    units::length.get_cgs_factor() / units::cgs_AU,
	    -SIGMASLOPE - 1.0 + 2.0 * FLARINGINDEX, FLARINGINDEX,
	    parameters::MU / constants::R * std::pow(ASPECTRATIO_REF, 2) *
		constants::G * mass * units::temperature.get_cgs_factor(),
	    units::temperature.get_cgs_symbol(),
	    units::length.get_cgs_factor() / units::cgs_AU,
	    -1.0 + 2.0 * FLARINGINDEX);

	for (unsigned int n_radial = 0;
	     n_radial <= data[t_data::ENERGY].get_max_radial(); ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal();
		 ++n_azimuthal) {

		const Pair cms =
		    data.get_planetary_system().get_center_of_mass();
		const double cms_x = cms.x;
		const double cms_y = cms.y;

		const double phi = (double)n_azimuthal * dphi;
		const double rmed = Rmed[n_radial];
		const double x = rmed * std::cos(phi) - cms_x;
		const double y = rmed * std::sin(phi) - cms_y;
		const double r = std::sqrt(x * x + y * y);

		const double energy =
		    1.0 / (ADIABATICINDEX - 1.0) * parameters::sigma0 *
		    std::pow(ASPECTRATIO_REF, 2) *
		    std::pow(r, -SIGMASLOPE - 1.0 + 2.0 * FLARINGINDEX) *
		    constants::G * mass;
		const double temperature_floor =
		    parameters::minimum_temperature *
		    units::temperature.get_inverse_cgs_factor();
		const double energy_floor =
		    temperature_floor *
		    data[t_data::DENSITY](n_radial, n_azimuthal) /
		    parameters::MU * constants::R / (ADIABATICINDEX - 1.0);

		data[t_data::ENERGY](n_radial, n_azimuthal) =
		    std::max(energy, energy_floor);
	    }
	}
    } break;

    case parameters::initialize_condition_read1D:
	logging::print_master(LOG_INFO "Loading Energy from '%s' (1D).\n",
			      parameters::energy_filename);
	data[t_data::ENERGY].read1D(parameters::energy_filename, true);
	break;

    case parameters::initialize_condition_read2D:
	logging::print_master(LOG_INFO "Loading Energy from '%s' (2D).\n",
			      parameters::energy_filename);
	data[t_data::ENERGY].read2D(parameters::energy_filename);
	break;

    case parameters::initialize_condition_shakura_sunyaev:
	die("Bad choice!"); // TODO: better explanation!
	break;
    }

    // profile damping outer?
    if (parameters::profile_cutoff_outer) {
	logging::print_master(
	    LOG_INFO "Damping Energy for r > %g %s over a range of %g %s\n",
	    parameters::profile_cutoff_point_outer *
		units::length.get_cgs_factor(),
	    units::length.get_cgs_symbol(),
	    parameters::profile_cutoff_width_outer *
		units::length.get_cgs_factor(),
	    units::length.get_cgs_symbol());
	for (unsigned int n_radial = 0;
	     n_radial <= data[t_data::ENERGY].get_max_radial(); ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal();
		 ++n_azimuthal) {
		// damp energy to 0 for r > profile_cutoff_point_outer
		double r;
		if (parameters::energy_initialize_condition ==
		    parameters::initialize_condition_profile_Nbody_centered) {
		    Pair cms = data.get_planetary_system().get_center_of_mass();
		    const double cms_x = cms.x;
		    const double cms_y = cms.y;

		    const double phi = (double)n_azimuthal * dphi;
		    const double rmed = Rmed[n_radial];
		    const double x = rmed * std::cos(phi) - cms_x;
		    const double y = rmed * std::sin(phi) - cms_y;
		    r = std::sqrt(x * x + y * y);
		} else {
		    r = Rmed[n_radial];
		}
		const double energy_damped =
		    data[t_data::ENERGY](n_radial, n_azimuthal) *
		    cutoff_outer(parameters::profile_cutoff_point_outer,
				 parameters::profile_cutoff_width_outer, r);
		const double temperature_floor =
		    parameters::minimum_temperature *
		    units::temperature.get_inverse_cgs_factor();
		const double energy_floor =
		    temperature_floor *
		    data[t_data::DENSITY](n_radial, n_azimuthal) /
		    parameters::MU * constants::R / (ADIABATICINDEX - 1.0);

		data[t_data::ENERGY](n_radial, n_azimuthal) =
		    std::max(energy_damped, energy_floor);
	    }
	}
    }

    // profile damping inner?
    if (parameters::profile_cutoff_inner) {
	logging::print_master(
	    LOG_INFO "Damping Energy for r < %g %s over a range of %g %s\n",
	    parameters::profile_cutoff_point_inner *
		units::length.get_cgs_factor(),
	    units::length.get_cgs_symbol(),
	    parameters::profile_cutoff_width_inner *
		units::length.get_cgs_factor(),
	    units::length.get_cgs_symbol());
	for (unsigned int n_radial = 0;
	     n_radial < data[t_data::ENERGY].get_size_radial(); ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::ENERGY].get_size_azimuthal();
		 ++n_azimuthal) {
		double r;
		if (parameters::energy_initialize_condition ==
		    parameters::initialize_condition_profile_Nbody_centered) {
		    Pair cms = data.get_planetary_system().get_center_of_mass();
		    const double cms_x = cms.x;
		    const double cms_y = cms.y;

		    const double phi = (double)n_azimuthal * dphi;
		    const double rmed = Rmed[n_radial];
		    const double x = rmed * std::cos(phi) - cms_x;
		    const double y = rmed * std::sin(phi) - cms_y;
		    r = std::sqrt(x * x + y * y);
		} else {
		    r = Rmed[n_radial];
		}
		// damp energy to 0 for r < profile_cutoff_point_inner
		const double energy_damped =
		    data[t_data::ENERGY](n_radial, n_azimuthal) *
		    cutoff_inner(parameters::profile_cutoff_point_inner,
				 parameters::profile_cutoff_width_inner, r);
		const double temperature_floor =
		    parameters::minimum_temperature *
		    units::temperature.get_inverse_cgs_factor();
		const double energy_floor =
		    temperature_floor *
		    data[t_data::DENSITY](n_radial, n_azimuthal) /
		    parameters::MU * constants::R / (ADIABATICINDEX - 1.0);

		data[t_data::ENERGY](n_radial, n_azimuthal) =
		    std::max(energy_damped, energy_floor);
	    }
	}
    }

    // set EnergyMed
    RefillEnergy(&data[t_data::ENERGY]);
}

/**
	Initialize gas velocities. This is described in 3.1.1 of Clement
   Baruteau's PhD Thesis
*/
void init_gas_velocities(t_data &data)
{
    double r, ri;
    double t1, t2, r1, r2;
    double vt_cent[MAX1D];

    if (parameters::sigma_initialize_condition ==
	parameters::initialize_condition_profile_Nbody_centered) {
	const double mass = data.get_planetary_system().get_mass();
	for (unsigned int n_radial = 0;
	     n_radial < data[t_data::V_RADIAL].get_size_radial(); ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::V_RADIAL].get_size_azimuthal();
		 ++n_azimuthal) {

		Pair cms = data.get_planetary_system().get_center_of_mass();
		const double cms_x = cms.x;
		const double cms_y = cms.y;

		const double phi = (double)n_azimuthal * dphi;
		double r;
		if (n_radial == data[t_data::V_AZIMUTHAL].Nrad) {
		    r = Rinf[data[t_data::V_AZIMUTHAL].Nrad - 1];
		} else {
		    r = Rinf[n_radial];
		}

		const double cell_x = r * std::cos(phi);
		const double cell_y = r * std::sin(phi);

		// Position in center of mass frame
		const double x_com = cell_x - cms_x;
		const double y_com = cell_y - cms_y;
		const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

		// pressure support correction
		double corr;
		double vr0;
		if (parameters::initialize_pure_keplerian) {
		    corr = 1.0;
		    vr0 = 0.0;
		} else {
		    corr = std::sqrt(
			1.0 - std::pow(ASPECTRATIO_REF, 2) *
				  std::pow(r_com, 2.0 * FLARINGINDEX) *
				  (1. + SIGMASLOPE - 2.0 * FLARINGINDEX));

		    /// Viscous speed
			const double h = ASPECTRATIO_REF * std::pow(r_com, FLARINGINDEX);
		    const double cs_iso =
			h *std::sqrt(constants::G * hydro_center_mass / r_com);
			const double H = h * r_com;
		    const double nu = ALPHAVISCOSITY * cs_iso * H;
		    vr0 = -3.0 * nu / r_com *
			  (-SIGMASLOPE + 2.0 * FLARINGINDEX + 1.0);
		}

		// Velocities in center of mass frame
		Pair v_cms =
		    data.get_planetary_system().get_center_of_mass_velocity();
		const double vr_com = vr0;
		const double vaz_com =
		    std::sqrt(constants::G * mass / r_com) * corr;

		const double vx_com =
		    (vr_com * x_com - vaz_com * y_com) / r_com;
		const double vy_com =
		    (vr_com * y_com + vaz_com * x_com) / r_com;

		// shift velocity from center of mass frame to primary frame
		const double vx = vx_com + v_cms.x;
		const double vy = vy_com + v_cms.y;

		const double vr = vx * std::cos(phi) + vy * std::sin(phi);
		data[t_data::V_RADIAL](n_radial, n_azimuthal) = vr;
	    }
	}

	for (unsigned int n_radial = 0;
	     n_radial <= data[t_data::V_AZIMUTHAL].get_max_radial();
	     ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::V_AZIMUTHAL].get_max_azimuthal();
		 ++n_azimuthal) {

		Pair cms = data.get_planetary_system().get_center_of_mass();
		const double cms_x = cms.x;
		const double cms_y = cms.y;

		const double phi = ((double)n_azimuthal - 0.5) * dphi;
		double r;
		if (n_radial == data[t_data::V_AZIMUTHAL].Nrad) {
		    r = Rmed[data[t_data::V_AZIMUTHAL].Nrad - 1];
		} else {
		    r = Rmed[n_radial];
		}

		const double cell_x = r * std::cos(phi);
		const double cell_y = r * std::sin(phi);

		// Position in center of mass frame
		const double x_com = cell_x - cms_x;
		const double y_com = cell_y - cms_y;
		const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

		// pressure support correction
		double corr;
		double vr0;

		if (parameters::initialize_pure_keplerian) {
		    corr = 1.0;
		    vr0 = 0.0;
		} else {
		    corr = std::sqrt(
			1.0 - std::pow(ASPECTRATIO_REF, 2) *
				  std::pow(r_com, 2.0 * FLARINGINDEX) *
				  (1. + SIGMASLOPE - 2.0 * FLARINGINDEX));

		    /// Viscous speed
			const double h = ASPECTRATIO_REF * std::pow(r_com, FLARINGINDEX);
			const double cs_iso =
			h *std::sqrt(constants::G * hydro_center_mass / r_com);
			const double H = h * r_com;
			const double nu = ALPHAVISCOSITY * cs_iso * H;
			vr0 = -3.0 * nu / r_com *
			  (-SIGMASLOPE + 2.0 * FLARINGINDEX + 1.0);
		}

		// Velocities in center of mass frame
		Pair v_cms =
		    data.get_planetary_system().get_center_of_mass_velocity();
		const double vr_com = vr0;
		const double vaz_com =
		    std::sqrt(constants::G * mass / r_com) * corr;

		const double vx_com =
		    (vr_com * x_com - vaz_com * y_com) / r_com;
		const double vy_com =
		    (vr_com * y_com + vaz_com * x_com) / r_com;

		// shift velocities from center of mass frame to primary frame
		const double vx = vx_com + v_cms.x;
		const double vy = vy_com + v_cms.y;

		const double vaz = vy * std::cos(phi) - vx * std::sin(phi);
		data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) = vaz;
	    }
	}

	return;
    }

    // Check if pure keplerian initialization is set
    if (parameters::initialize_pure_keplerian) {
	for (unsigned int n_radial = 0;
	     n_radial <= data[t_data::V_AZIMUTHAL].get_max_radial();
	     ++n_radial) {
	    // TODO: This should be Rinf but seems to produce incorrect data
	    if (n_radial == data[t_data::V_AZIMUTHAL].Nrad) {
		r = Rmed[data[t_data::V_AZIMUTHAL].Nrad - 1];
	    } else {
		r = Rmed[n_radial];
	    }

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::V_AZIMUTHAL].get_max_azimuthal();
		 ++n_azimuthal) {
		data[t_data::V_RADIAL](n_radial, n_azimuthal) = 0.0;
		data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) =
		    std::sqrt(constants::G * hydro_center_mass / r);
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
	std::memset(vt_int, 0, MAX1D * sizeof(*vt_int));

	/* vt_int \equiv rOmega� = grad(P)/sigma +  \partial(phi)/\partial(r)  -
	 * acc_sg_radial */
	mpi_make1Dprofile(data[t_data::PRESSURE].Field, GLOBAL_bufarray);

	// get global SigmaMed
	t_radialarray SigmaMedGlobal(GlobalNRadial);
	MPI_Status status;
	MPI_Request request;

	for (int i = 0; i < CPU_Number; ++i) {
	    if (i == CPU_Rank) {
		// receive data from prev CPU
		if (CPU_Rank > 0) {
		    MPI_Irecv(&SigmaMedGlobal[0], IMIN + CPUOVERLAP - 1,
			      MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD,
			      &request);
		    MPI_Wait(&request, &status);
		}
		// insert own data
		for (unsigned int j = IMIN + (CPU_Rank == 0 ? 0 : CPUOVERLAP);
		     j <= IMAX - (CPU_Rank == CPU_Highest ? 0 : CPUOVERLAP);
		     ++j) {
		    SigmaMedGlobal[j] =
			SigmaMed[Zero_or_active + j - IMIN -
				 (CPU_Rank == 0 ? 0 : CPUOVERLAP)];
		}
		// send data to next CPU
		if (CPU_Rank < CPU_Highest) {
		    MPI_Isend(&SigmaMedGlobal[0], IMAX - CPUOVERLAP, MPI_DOUBLE,
			      CPU_Next, 0, MPI_COMM_WORLD, &request);
		    MPI_Wait(&request, &status);
		}
	    }
	}

	// now CPU_Highst has all data. broadcast to everybody
	MPI_Bcast(&SigmaMedGlobal[0], GlobalNRadial, MPI_DOUBLE, CPU_Highest,
		  MPI_COMM_WORLD);

	/* global axisymmetric pressure field, known by all cpus */
	for (unsigned int i = 1; i < GlobalNRadial; i++) {
	    vt_int[i] = (GLOBAL_bufarray[i] - GLOBAL_bufarray[i - 1]) /
			    (.5 * (SigmaMedGlobal[i] + SigmaMedGlobal[i - 1])) /
			    (GlobalRmed[i] - GlobalRmed[i - 1]) +
			constants::G * hydro_center_mass *
			    (1.0 / GlobalRmed[i - 1] - 1.0 / GlobalRmed[i]) /
			    (GlobalRmed[i] - GlobalRmed[i - 1]);
	}

	/* Case of a disk with self-gravity */
	if (parameters::self_gravity) { // Better test with CL rigid!
	    double *GLOBAL_AxiSGAccr =
		(double *)malloc(sizeof(double) * GlobalNRadial);
	    mpi_make1Dprofile(selfgravity::g_radial, GLOBAL_AxiSGAccr);

	    for (unsigned int i = 1; i < GlobalNRadial; i++) {
		vt_int[i] -=
		    ((Radii[i] - GlobalRmed[i - 1]) * GLOBAL_AxiSGAccr[i] +
		     (GlobalRmed[i] - Radii[i]) * GLOBAL_AxiSGAccr[i - 1]) /
		    (GlobalRmed[i] - GlobalRmed[i - 1]);
	    }
	    free(GLOBAL_AxiSGAccr);
	}

	for (unsigned int i = 1; i < GlobalNRadial; i++)
	    vt_int[i] = std::sqrt(vt_int[i] * Radii[i]) - Radii[i] * OmegaFrame;

	t1 = vt_cent[0] = vt_int[1] + .75 * (vt_int[1] - vt_int[2]);
	r1 = ConstructSequence(vt_cent, vt_int, GlobalNRadial);
	vt_cent[0] += .25 * (vt_int[1] - vt_int[2]);
	t2 = vt_cent[0];
	r2 = ConstructSequence(vt_cent, vt_int, GlobalNRadial);
	t1 = t1 - r1 / (r2 - r1) * (t2 - t1);
	vt_cent[0] = t1;
	ConstructSequence(vt_cent, vt_int, GlobalNRadial);
	vt_cent[GlobalNRadial] = vt_cent[GlobalNRadial - 1];
    }

    // Initialization with self-gravity, without exact centrifugal balance
	if (parameters::self_gravity && !CentrifugalBalance)
	selfgravity::init_azimuthal_velocity(data[t_data::V_AZIMUTHAL]);

	for (unsigned int n_radial = 1;
	 n_radial <= data[t_data::V_AZIMUTHAL].get_max_radial(); ++n_radial) {
	if (n_radial == data[t_data::V_AZIMUTHAL].Nrad) {
	    r = Rmed[data[t_data::V_AZIMUTHAL].Nrad - 1];
	    ri = Rinf[data[t_data::V_AZIMUTHAL].Nrad - 1];
	} else {
	    r = Rmed[n_radial];
	    ri = Rinf[n_radial];
	}

	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::V_AZIMUTHAL].get_max_azimuthal();
	     ++n_azimuthal) {
	    if (!parameters::self_gravity) {
		// v_azimuthal = Omega_K * r * (...)
		data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) =
		    r * calculate_omega_kepler(r) *
		    std::sqrt(1.0 - std::pow(ASPECTRATIO_REF, 2) *
					std::pow(r, 2.0 * FLARINGINDEX) *
					(1. + SIGMASLOPE - 2.0 * FLARINGINDEX));
	    }

	    data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) -= OmegaFrame * r;

	    if (CentrifugalBalance)
		data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) =
		    vt_cent[n_radial + IMIN];

	    if (n_radial == data[t_data::V_RADIAL].Nrad) {
		data[t_data::V_RADIAL](n_radial, n_azimuthal) = 0.0;
	    } else {
		data[t_data::V_RADIAL](n_radial, n_azimuthal) =
		    IMPOSEDDISKDRIFT * parameters::sigma0 / SigmaInf[n_radial] /
		    ri;

		if (!parameters::initialize_vradial_zero) {
		    if (ViscosityAlpha) {
			const double nu = 0.5 * (data[t_data::VISCOSITY](n_radial, n_azimuthal) + data[t_data::VISCOSITY](n_radial-1, n_azimuthal));
			data[t_data::V_RADIAL](n_radial, n_azimuthal) -=
				3.0 * nu / Rinf[n_radial] *
			    (-SIGMASLOPE + 2.0 * FLARINGINDEX + 1.0);

		    } else {
			data[t_data::V_RADIAL](n_radial, n_azimuthal) -=
			    3.0 *
			    data[t_data::VISCOSITY](n_radial, n_azimuthal) / r *
			    (-SIGMASLOPE + .5);
		    }
		}
	    }
	}
    }

    /// set VRadial for innermost and outermost ring to 0
    for (unsigned int n_azimuthal = 0;
	 n_azimuthal <= data[t_data::V_RADIAL].get_max_azimuthal();
	 ++n_azimuthal) {
	data[t_data::V_RADIAL](0, n_azimuthal) = 0.0;
	data[t_data::V_RADIAL](data[t_data::V_RADIAL].Nrad, n_azimuthal) = 0.0;
    }
}
