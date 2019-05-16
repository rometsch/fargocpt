/**
	\file parameters.cpp
	\author Tobias Mueller <Tobias_Mueller@twam.info>
*/
#include <ctype.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <math.h>

#include "parameters.h"
#include "output.h"
#include "config.h"
#include "logging.h"
#include "LowTasks.h"
#include "viscosity.h"
#include "util.h"
#include "string.h"
#include "global.h"
#include "constants.h"
#include "boundary_conditions.h"

namespace parameters {


t_radial_grid radial_grid_type;
const char* radial_grid_names[] = {"arithmetic", "logarithmic", "exponential"};

t_boundary_condition boundary_inner;
t_boundary_condition boundary_outer;
bool domegadr_zero;

bool damping;
double damping_inner_limit;
double damping_outer_limit;
double damping_time_factor;


int damping_energy_id;
std::vector<t_DampingType> damping_vector;

double minimum_temperature;
double maximum_temperature;
double MU;
bool heating_viscous_enabled;
double heating_viscous_factor;
bool heating_star_enabled;
double heating_star_factor;
double heating_star_ramping_time;
bool heating_star_simple;

double cooling_radiative_factor;
bool cooling_radiative_enabled;
bool cooling_beta_enabled;
double cooling_beta;

bool radiative_diffusion_enabled;
double radiative_diffusion_omega;
bool radiative_diffusion_omega_auto_enabled;
unsigned int radiative_diffusion_max_iterations;

t_initialize_condition sigma_initialize_condition;
char *sigma_filename = NULL;
int random_seed;
bool sigma_randomize;
double sigma_random_factor;
double sigma_floor;
double sigma_feature_size;
bool sigma_adjust;
double sigma_discmass;
double sigma0;
t_initialize_condition energy_initialize_condition;
char *energy_filename = NULL;

t_artificial_viscosity artificial_viscosity;
double artificial_viscosity_factor;
bool artificial_viscosity_dissipation;

bool calculate_disk;

bool massoverflow;
unsigned int mof_planet;
double mof_sigma;
double mof_value;

bool profile_damping;
double profile_damping_point;
double profile_damping_width;

bool disk_feedback;

bool integrate_planets;

double density_factor;
double tau_factor;
double kappa_factor;

bool self_gravity;

bool write_torques;

bool write_disk_quantities;
bool write_lightcurves;
bool write_at_every_timestep;
std::vector<double> lightcurves_radii;
bool write_massflow;

unsigned int log_after_steps;
double log_after_real_seconds;	

t_opacity opacity;

double thickness_smoothing;
double thickness_smoothing_sg;

bool initialize_pure_keplerian;
bool initialize_vradial_zero;

double star_radius;
double star_temperature;

double radial_viscosity_factor;
double vrad_fraction_of_kepler;
double stellar_rotation_rate;
double mass_accretion_rate;

unsigned int zbuffer_size;
double zbuffer_maxangle;

double CFL;

double L0;
double M0;

unsigned int number_of_particles;
bool integrate_particles;
double particle_radius;
double particle_density;
double particle_slope;
double particle_minimum_radius;
double particle_maximum_radius;
double particle_escape_radius;
bool particle_gas_drag_enabled;
bool particle_disk_gravity_enabled;

// for constant opacity
double kappa_const = 1.0;

	void exitOnDeprecatedSetting(std::string setting_name, std::string reason, std::string instruction) {
	if (config::key_exists(setting_name.c_str())) {
		std::string warn = "Deprecated setting '" + setting_name + "' found in .par file. " + reason + " " + instruction+ "\n";
		logging::print_master((LOG_ERROR + warn).c_str());
		die("Config Error");
	}
}

t_DampingType write_damping_type(t_damping_type type_inner, t_damping_type type_outer, t_data::t_polargrid_type quantity, t_data::t_polargrid_type quantity0, std::string description)
{
	t_DampingType damping_type;
	damping_type.array_to_damp = quantity;
	damping_type.array_with_damping_values = quantity0;
	std::string description_inner;
	std::string description_outer;

	switch (type_inner)
	{
		case damping_none:
			damping_type.inner_damping_function = nullptr;
			description_inner = "Damping " + description + " is disabled at inner boundary.";
			break;
		case damping_mean:
			damping_type.inner_damping_function = &boundary_conditions::damping_single_inner_mean;
			description_inner = "Damping " + description + " to mean value at inner boundary.";
			break;
		case damping_initial:
			damping_type.inner_damping_function = &boundary_conditions::damping_single_inner;
			description_inner = "Damping " + description + " to initial value at inner boundary.";
			break;
		case damping_zero:
			damping_type.inner_damping_function = &boundary_conditions::damping_single_inner_zero;
			description_inner = "Damping " + description + " to zero at inner boundary.";
			break;
	}

	damping_type.description_inner = description_inner;

	switch (type_outer)
	{
		case damping_none:
			damping_type.outer_damping_function = nullptr;
			description_outer = "Damping " + description + " is disabled at outer boundary.";
			break;
		case damping_mean:
			damping_type.outer_damping_function = &boundary_conditions::damping_single_outer_mean;
			description_outer = "Damping " + description + " to mean value at outer boundary.";
			break;
		case damping_initial:
			damping_type.outer_damping_function = &boundary_conditions::damping_single_outer;
			description_outer = "Damping " + description + " to initial value at outer boundary.";
			break;
		case damping_zero:
			damping_type.outer_damping_function = &boundary_conditions::damping_single_outer_zero;
			description_outer = "Damping " + description + " to zero at outer boundary.";
			break;
	}

	damping_type.description_outer = description_outer;

	return damping_type;
}


void read(char* filename, t_data &data)
{
	// read config from
	if (config::read_config_from_file(filename) == -1) {
		logging::print_master(LOG_ERROR "Cannot read config file '%s'!\n", filename);
		PersonalExit(EXIT_FAILURE);
	}

	/* grid */
	NRadial = config::value_as_unsigned_int_default("NRAD",64);
	NAzimuthal = config::value_as_unsigned_int_default("NSEC",64);
	RMIN = config::value_as_double_default("RMIN",1.0);
	RMAX = config::value_as_double_default("RMAX",1.0);

	switch (tolower(*config::value_as_string_default("RadialSpacing","ARITHMETIC"))) {
		case 'a': // arithmetic
			radial_grid_type = arithmetic_spacing;
			break;
		case 'l': // logarithmic;
			radial_grid_type = logarithmic_spacing;
			break;
		case 'e': // exponential;
			radial_grid_type = exponential_spacing;
			break;
		default:
			die("Invalid setting for RadialSpacing");
	}

	// units
	L0 = config::value_as_double_default("l0",1.0);
	M0 = config::value_as_double_default("m0",1.0);

	// output settings
	data[t_data::DENSITY].set_write(config::value_as_bool_default("WriteDensity", true));
	data[t_data::V_RADIAL].set_write(config::value_as_bool_default("WriteVelocity", true));
	data[t_data::V_AZIMUTHAL].set_write(config::value_as_bool_default("WriteVelocity", true));
	data[t_data::ENERGY].set_write(config::value_as_bool_default("WriteEnergy", true));
	data[t_data::TEMPERATURE].set_write(config::value_as_bool_default("WriteTemperature", false));
	data[t_data::SOUNDSPEED].set_write(config::value_as_bool_default("WriteSoundSpeed", false));
	data[t_data::PRESSURE].set_write(config::value_as_bool_default("WritePressure", false));
	data[t_data::TOOMRE].set_write(config::value_as_bool_default("WriteToomre", false));
	data[t_data::TOOMRE_1D].set_write(config::value_as_bool_default("WriteRadialToomre", false));
	data[t_data::QPLUS].set_write(config::value_as_bool_default("WriteQPlus", false));
	data[t_data::QMINUS].set_write(config::value_as_bool_default("WriteQMinus", false));
	data[t_data::KAPPA].set_write(config::value_as_bool_default("WriteKappa", false));
	data[t_data::TAU_COOL].set_write(config::value_as_bool_default("WriteTauCool", false));
	data[t_data::ALPHA_GRAV].set_write(config::value_as_bool_default("WriteAlphaGrav", false));
	data[t_data::ALPHA_GRAV_1D].set_write(config::value_as_bool_default("WriteRadialAlphaGrav", false));
	data[t_data::ALPHA_GRAV_MEAN].set_write(config::value_as_bool_default("WriteAlphaGravMean", false));
	data[t_data::ALPHA_GRAV_MEAN_1D].set_write(config::value_as_bool_default("WriteRadialAlphaGravMean", false));
	data[t_data::ALPHA_REYNOLDS].set_write(config::value_as_bool_default("WriteAlphaReynolds", false));
	data[t_data::ALPHA_REYNOLDS_1D].set_write(config::value_as_bool_default("WriteRadialAlphaReynolds", false));
	data[t_data::ALPHA_REYNOLDS_MEAN].set_write(config::value_as_bool_default("WriteAlphaReynoldsMean", false));
	data[t_data::ALPHA_REYNOLDS_MEAN_1D].set_write(config::value_as_bool_default("WriteRadialAlphaReynoldsMean", false));
	data[t_data::VISCOSITY].set_write(config::value_as_bool_default("WriteViscosity", false));
	data[t_data::DIV_V].set_write(config::value_as_bool_default("WriteDivV", false));
	data[t_data::ECCENTRICITY].set_write(config::value_as_bool_default("WriteEccentricity", false));
	data[t_data::T_REYNOLDS].set_write(config::value_as_bool_default("WriteTReynolds", false));
	data[t_data::T_GRAVITATIONAL].set_write(config::value_as_bool_default("WriteTGravitational", false));
	data[t_data::P_DIVV].set_write(config::value_as_bool_default("WritepDV", false));
	data[t_data::TAU].set_write(config::value_as_bool_default("WriteTau", false));
	data[t_data::ASPECTRATIO].set_write(config::value_as_bool_default("WriteAspectRatio", false));
	data[t_data::VISIBILITY].set_write(config::value_as_bool_default("WriteVisibility", false));
	data[t_data::LUMINOSITY_1D].set_write(config::value_as_bool_default("WriteRadialLuminosity", false));
	data[t_data::DISSIPATION_1D].set_write(config::value_as_bool_default("WriteRadialDissipation", false));
	data[t_data::TAU_EFF].set_write(config::value_as_bool_default("WriteVerticalOpticalDepth", false));

	write_torques = config::value_as_bool_default("WriteTorques", false);

	write_disk_quantities = config::value_as_bool_default("WriteDiskQuantities", false);
	write_at_every_timestep = config::value_as_bool_default("WriteAtEveryTimestep", false);
	write_lightcurves = config::value_as_bool_default("WriteLightCurves", false);

	write_massflow = config::value_as_bool_default("WriteMassFlow", false);

	log_after_steps = config::value_as_unsigned_int_default("LogAfterSteps",0);
	log_after_real_seconds = config::value_as_double_default("LogAfterRealSeconds",600.0);

	// parse light curve radii
	if (config::key_exists("WriteLightCurvesRadii")) {
		// get light curves radii string
		char* lightcurves_radii_string = new char[strlen(config::value_as_string("WriteLightCurvesRadii"))];
		strcpy(lightcurves_radii_string, config::value_as_string("WriteLightCurvesRadii"));

		char *ptr = strtok(lightcurves_radii_string, " ,");
		while (ptr != NULL) {
			double value;
			value = strtod(ptr, NULL);
			if ((value > RMIN) && (value < RMAX)) {
				lightcurves_radii.push_back(value);
			}
			ptr = strtok(NULL, " ,");
		}

		lightcurves_radii.push_back(RMIN);
		lightcurves_radii.push_back(RMAX);

		delete [] lightcurves_radii_string;

		// sort vector
		sort(lightcurves_radii.begin(), lightcurves_radii.end());
	}

	// boundary conditions
	switch (tolower(*config::value_as_string_default("InnerBoundary","Open"))) {
		case 'o':
			boundary_inner = boundary_condition_open;
			break;
		case 'n':
			boundary_inner = boundary_condition_nonreflecting;
			break;
		case 'e':
			boundary_inner = boundary_condition_evanescent;
			break;
		case 'r':
			boundary_inner = boundary_condition_reflecting;
			break;
		case 'v':
			boundary_inner = boundary_condition_viscous_outflow;
			break;
		case 'b':
			boundary_inner = boundary_condition_boundary_layer;
			break;
		case 'k':
			boundary_inner = boundary_condition_keplerian;
			break;
		default:
			die("Invalid setting for InnerBoundary: %s",config::value_as_string_default("InnerBoundary","Open"));
	}

	switch (tolower(*config::value_as_string_default("OuterBoundary","Open"))) {
		case 'o':
			boundary_outer = boundary_condition_open;
			break;
		case 'n':
			boundary_outer = boundary_condition_nonreflecting;
			break;
		case 'e':
			boundary_outer = boundary_condition_evanescent;
			break;
		case 'r':
			boundary_outer = boundary_condition_reflecting;
			break;
		case 'v':
			boundary_outer = boundary_condition_viscous_outflow;
			break;
		case 'b':
			boundary_outer = boundary_condition_boundary_layer;
			break;
		case 'k':
			boundary_outer = boundary_condition_keplerian;
			break;
		default:
			die("Invalid setting for OuterBoundary: %s",config::value_as_string_default("OuterBoundary","Open"));
	}

	domegadr_zero = config::value_as_bool_default("DomegaDrZero", false);
	damping = config::value_as_bool_default("Damping", false);

	damping_inner_limit = config::value_as_double_default("DampingInnerLimit",1.05);
	if (damping_inner_limit < 1) {
		die("DampingInnerLimit must not be <1\n");
	}
	damping_outer_limit = config::value_as_double_default("DampingOuterLimit",0.95);
	if (damping_outer_limit > 1) {
		die("DampingOuterLimit must not be >1\n");
	}
	damping_time_factor = config::value_as_double_default("DampingTimeFactor",1.0);


	t_damping_type tmp_damping_inner;
	t_damping_type tmp_damping_outer;

	if(config::key_exists("DampingVRadial"))
		die("DampingVRadial flag is decrepated used DampingVRadialInner and DampingVRadialOuter instead!");

	tmp_damping_inner = config::value_as_boudary_damping_default("DampingVRadialInner", "None");
	tmp_damping_outer = config::value_as_boudary_damping_default("DampingVRadialOuter", "None");

	damping_vector.push_back(write_damping_type(tmp_damping_inner, tmp_damping_outer, t_data::V_RADIAL, t_data::V_RADIAL0, "VRadial"));


	if(config::key_exists("DampingVAzimuthal"))
		die("DampingVRadial flag is decrepated used DampingVAzimuthalInner and DampingVAzimuthalOuter instead!");

	tmp_damping_inner = config::value_as_boudary_damping_default("DampingVAzimuthalInner", "None");
	tmp_damping_outer= config::value_as_boudary_damping_default("DampingVAzimuthalOuter", "None");

	damping_vector.push_back(write_damping_type(tmp_damping_inner, tmp_damping_outer, t_data::V_AZIMUTHAL, t_data::V_AZIMUTHAL0, "VAzimuthal"));


	if(config::key_exists("DampingSurfaceDensity"))
		die("DampingSurfaceDensity flag is decrepated used DampingSurfaceDensityInner and DampingSurfaceDensityOuter instead!");

	tmp_damping_inner = config::value_as_boudary_damping_default("DampingSurfaceDensityInner", "None");
	tmp_damping_outer = config::value_as_boudary_damping_default("DampingSurfaceDensityOuter", "None");


	damping_vector.push_back(write_damping_type(tmp_damping_inner, tmp_damping_outer, t_data::DENSITY, t_data::DENSITY0, "SurfaceDensity"));


	tmp_damping_inner = config::value_as_boudary_damping_default("DampingEnergy", "None");
	if(config::key_exists("DampingEnergy"))
		die("DampingEnergy flag is decrepated used DampingEnergyInner and DampingEnergyOuter instead!");

	tmp_damping_inner = config::value_as_boudary_damping_default("DampingEnergyInner", "None");
	tmp_damping_outer = config::value_as_boudary_damping_default("DampingEnergyOuter", "None");

	damping_vector.push_back(write_damping_type(tmp_damping_inner, tmp_damping_outer, t_data::ENERGY, t_data::ENERGY0, "Energy"));
	damping_energy_id = (int)damping_vector.size() - 1;



	calculate_disk = config::value_as_bool_default("DISK", 1);

	MU = config::value_as_double_default("mu",1.0);
	minimum_temperature = config::value_as_double_default("MinimumTemperature", 0);
    maximum_temperature = config::value_as_double_default("MaximumTemperature", 1.0e300);

	// TODO: remove temporary warning
	if (config::key_exists("HeatingViscous") == false) {
		die("please specify HeatingViscous in config file");
	}
	heating_viscous_enabled = config::value_as_bool_default("HeatingViscous", false);
	heating_viscous_factor = config::value_as_double_default("HeatingViscousFactor", 1.0);
	heating_star_enabled = config::value_as_bool_default("HeatingStar", false);
	heating_star_factor = config::value_as_double_default("HeatingStarFactor", 1.0);
	heating_star_ramping_time = config::value_as_double_default("HeatingStarRampingTime", 0.0);
	heating_star_simple = config::value_as_bool_default("HeatingStarSimple", false);

	radiative_diffusion_enabled = config::value_as_bool_default("RadiativeDiffusion", false);
	radiative_diffusion_omega = config::value_as_double_default("RadiativeDiffusionOmega", 1.5);
	radiative_diffusion_omega_auto_enabled = config::value_as_bool_default("RadiativeDiffusionAutoOmega", false);
	radiative_diffusion_max_iterations = config::value_as_unsigned_int_default("RadiativeDiffusionMaxIterations", 50000);

	zbuffer_size = config::value_as_unsigned_int_default("zbufferSize", 100);
	zbuffer_maxangle = config::value_as_double_default("zbufferMaxAngle", 10.0/180.0*PI);

	cooling_radiative_factor = config::value_as_double_default("CoolingRadiativeFactor", 1.0);
	cooling_radiative_enabled = config::value_as_bool_default("CoolingRadiativeLocal", false);
	cooling_beta_enabled = config::value_as_bool_default("CoolingBetaLocal", false);
	cooling_beta = config::value_as_double_default("CoolingBeta", 1.0);

	// initialisation
	initialize_pure_keplerian = config::value_as_bool_default("InitializePureKeplerian", false);
	initialize_vradial_zero = config::value_as_bool_default("InitializeVradialZero", false);

	switch (tolower(*config::value_as_string_default("SigmaCondition","Profile"))) {
		case 'p': // Profile
			sigma_initialize_condition = initialize_condition_profile;
			break;
		case '1': // 1D
			sigma_initialize_condition = initialize_condition_read1D;
			break;
		case '2': // 2D
			sigma_initialize_condition = initialize_condition_read2D;
			break;
		case 's': // Boundary Layer: initialize w/ SS-73-Standard-Solution
			sigma_initialize_condition = initialize_condition_shakura_sunyaev;
			break;
		default:
			die("Invalid setting for SigmaCondition: %s", config::value_as_string_default("SigmaCondition","Profile"));
	}

	if (config::key_exists("SigmaFilename")) {
		if (strlen(config::value_as_string_default("SigmaFilename","")) > 0) {
			sigma_filename = new char[strlen(config::value_as_string_default("SigmaFilename",""))];
			strcpy(sigma_filename, config::value_as_string_default("SigmaFilename",""));
		}
	}

	switch (tolower(*config::value_as_string_default("EnergyCondition","Profile"))) {
		case 'p': // Profile
			energy_initialize_condition = initialize_condition_profile;
			break;
		case '1': // 1D
			energy_initialize_condition = initialize_condition_read1D;
			break;
		case '2': // 2D
			energy_initialize_condition = initialize_condition_read2D;
			break;
		case 's': // Boundary Layer: initialize w/ SS-73-Standard-Solution
			energy_initialize_condition = initialize_condition_shakura_sunyaev;
			break;
		default:
			die("Invalid setting for EnergyCondition: %s", config::value_as_string_default("EnergyCondition","Profile"));
	}

	if (config::key_exists("EnergyFilename")) {
		if (strlen(config::value_as_string_default("EnergyFilename","")) > 0) {
			energy_filename = new char[strlen(config::value_as_string_default("EnergyFilename",""))];
			strcpy(energy_filename, config::value_as_string_default("EnergyFilename",""));
		}
	}

    random_seed = config::value_as_int_default("RandomSeed", 0);
	sigma_randomize = config::value_as_bool_default("RandomSigma", 0);
	sigma_random_factor = config::value_as_double_default("RandomFactor", 0.1);
    sigma_feature_size = config::value_as_double_default("FeatureSize", (RMAX - RMIN)/150);
	sigma_floor = config::value_as_double_default("SigmaFloor", 1e-9);
	sigma0 = config::value_as_double_default("SIGMA0", 173.);
	sigma_adjust = config::value_as_bool_default("SetSigma0", false);
	sigma_discmass = config::value_as_double_default("discmass", 0.01);
	density_factor = config::value_as_double_default("DensityFactor", 2.0);

	tau_factor = config::value_as_double_default("TauFactor", 1.0);
	kappa_factor = config::value_as_double_default("KappaFactor", 1.0);

	// artificial visocisty
	switch (tolower(*config::value_as_string_default("ArtificialViscosity","SN"))) {
		case 'n': // none
			artificial_viscosity = artificial_viscosity_none;
			break;
		case 't': // TW
			artificial_viscosity = artificial_viscosity_TW;
			break;
		case 's': // SN
			artificial_viscosity = artificial_viscosity_SN;
			break;
		default:
			die("Invalid setting for ArtificialViscosity: %s",config::value_as_string_default("ArtificialViscosity","SN"));
	}
	artificial_viscosity_dissipation = config::value_as_bool_default("ArtificialViscosityDissipation", true);
	artificial_viscosity_factor = config::value_as_double_default("ArtificialViscosityFactor", 1.41);
	// warning
	if (config::key_exists("CVNR")) {
			die("Parameter CVNR has been renamed to ArtificialViscosityFactor");
	}

	//
	thickness_smoothing = config::value_as_double_default("ThicknessSmoothing", 0.0);
	thickness_smoothing_sg = config::value_as_double_default("ThicknessSmoothingSG", thickness_smoothing);
	integrate_planets = config::value_as_bool_default("IntegratePlanets", true);

	// mass overflow
	massoverflow = config::value_as_bool_default("massoverflow", false);
	mof_planet = config::value_as_int_default("mofplanet", 0);
	mof_sigma = config::value_as_double_default("mofsigma", 1.0);
	mof_value = config::value_as_double_default("mofvalue", 10E-9);

	// profile damping
	profile_damping = config::value_as_bool_default("ProfileDamping", false);
	profile_damping_point = config::value_as_double_default("ProfileDampingPoint", 0.0);
	profile_damping_width = config::value_as_double_default("ProfileDampingWidth", 1.0);

    exitOnDeprecatedSetting("FeelsDisk", "Replaced by parameter DiskFeedback for clarification since it affects more than just the star.", "Please replace 'FeelsDisk' by 'DiskFeedback'.");
	disk_feedback = config::value_as_bool_default("DiskFeedback", true);

	// self gravity
	self_gravity = config::value_as_bool_default("SelfGravity",0);

	// opacity
	switch (tolower(*config::value_as_string_default("Opacity","Lin"))) {
		case 'l': // Lin
			opacity = opacity_lin;
			break;
		case 'b': // Bell
			opacity = opacity_bell;
			break;
		case 'z': // Zhu
			opacity = opacity_zhu;
			break;
		case 'k': // Kramers
			opacity = opacity_kramers;
			break;
		case 'c': // Constant
			opacity = opacity_const_op;
			kappa_const = config::value_as_double_default("KappaConst", 1.0);
			break;
		default:
			die("Invalid setting for Opacity: %s",config::value_as_string_default("Opacity","Lin"));
	}

	// star parameters
	star_temperature = config::value_as_double_default("StarTemperature", 5778);
	star_radius = config::value_as_double_default("StarRadius", 0.009304813);

	// boundary layer parameters
	radial_viscosity_factor = config::value_as_double_default("RadialViscosityFactor", 1.);
	vrad_fraction_of_kepler = config::value_as_double_default("VRadIn", 1.6e-3);
	stellar_rotation_rate = config::value_as_double_default("StellarRotation", 0.1);
	mass_accretion_rate = config::value_as_double_default("MassAccretionRate", 1.e-9);

	CFL = config::value_as_double_default("CFL", 0.5);

	// particles
	integrate_particles = config::value_as_bool_default("IntegrateParticles", false);
	number_of_particles = config::value_as_unsigned_int_default("NumberOfParticles", 0);
	particle_radius = config::value_as_double_default("ParticleRadius", 100.0);
	particle_density = config::value_as_double_default("ParticleDensity", 2.65);
	particle_slope = config::value_as_double_default("ParticleSlope", 0.0);
	particle_minimum_radius = config::value_as_double_default("ParticleMinimumRadius", RMIN);
	particle_maximum_radius = config::value_as_double_default("ParticleMaximumRadius", RMAX);
	particle_escape_radius = config::value_as_double_default("ParticleEscapeRadius", particle_maximum_radius);
	particle_gas_drag_enabled = config::value_as_bool_default("ParticleGasDragEnabled", true);
	particle_disk_gravity_enabled = config::value_as_bool_default("ParticleDiskGravityEnabled", false);
}

void apply_units() {
	star_temperature /= units::temperature.get_cgs_factor();
	mass_accretion_rate = mass_accretion_rate*(units::cgs_Msol/units::cgs_Year)*1./units::mass_accretion_rate.get_cgs_factor();
	sigma0 /= units::surface_density.get_cgs_factor();
	particle_radius /= units::length.get_cgs_factor();
	particle_density /= units::density.get_cgs_factor();
}

void summarize_parameters()
{
	// artifical viscosity
	switch (artificial_viscosity) {
		case artificial_viscosity_none:
			logging::print_master(LOG_INFO "Using no artificial viscosity.\n");
			break;
		case artificial_viscosity_TW:
			logging::print_master(LOG_INFO "Using Tscharnuter-Winkler (1979) artificial viscosity with C = %lf.\n", parameters::artificial_viscosity_factor);
			logging::print_master(LOG_INFO "Artificial viscosity is %s for dissipation.\n", parameters::artificial_viscosity_dissipation ? "used" : "not used");
			break;
		case artificial_viscosity_SN:
			logging::print_master(LOG_INFO "Using Stone-Norman (1991, ZEUS-2D) artificial viscosity with C = %lf.\n", parameters::artificial_viscosity_factor);
			logging::print_master(LOG_INFO "Artificial viscosity is %s for dissipation.\n", parameters::artificial_viscosity_dissipation ? "used" : "not used");
			break;
	}

	// boundary conditions
	switch (boundary_inner) {
		case boundary_condition_open:
			logging::print_master(LOG_INFO "Using 'open boundary condition' at inner boundary.\n");
			break;
		case boundary_condition_reflecting:
			logging::print_master(LOG_INFO "Using 'reflecting boundary condition' at inner boundary.\n");
			break;
		case boundary_condition_nonreflecting:
			logging::print_master(LOG_INFO "Using 'nonreflecting boundary condition' at inner boundary.\n");
			break;
		case boundary_condition_evanescent:
			logging::print_master(LOG_INFO "Using 'evanescent boundary condition' at inner boundary.\n");
			break;
		case boundary_condition_viscous_outflow:
			logging::print_master(LOG_INFO "Using 'viscous outflow boundary condition' at inner boundary.\n");
			break;
		case boundary_condition_boundary_layer:
			logging::print_master(LOG_INFO "Using 'boundary layer boundary conditions' at inner boundary.\n");
			break;
        case boundary_condition_keplerian:
            logging::print_master(LOG_INFO "Using 'keplarian boundary conditions' at inner boundary.\n");
            break;
	}

	switch (boundary_outer) {
		case boundary_condition_open:
			logging::print_master(LOG_INFO "Using 'open boundary condition' at outer boundary.\n");
			break;
		case boundary_condition_reflecting:
			logging::print_master(LOG_INFO "Using 'reflecting boundary condition' at outer boundary.\n");
			break;
		case boundary_condition_nonreflecting:
			logging::print_master(LOG_INFO "Using 'nonreflecting boundary condition' at outer boundary.\n");
			break;
		case boundary_condition_evanescent:
			logging::print_master(LOG_INFO "Using 'evanescent boundary condition' at outer boundary.\n");
			break;
		case boundary_condition_viscous_outflow:
			logging::print_master(LOG_INFO "Using 'viscous outflow boundary condition' at outer boundary.\n");
			break;	
		case boundary_condition_boundary_layer:
			logging::print_master(LOG_INFO "Using 'boundary layer boundary conditions' at outer boundary.\n");
			break;
        case boundary_condition_keplerian:
            logging::print_master(LOG_INFO "Using 'keplarian boundary conditions' at inner boundary.\n");
            break;
	}

	// Mass Transfer
	if (parameters::massoverflow) {
		logging::print_master(LOG_INFO "Mass Transfer of %g M_0/orbit will be spread on %i gridcells (sigma = %g).\n",parameters::mof_value, int(NAzimuthal*3.0*parameters::mof_sigma),mof_sigma);
	}

	// Boundary layer
	if( boundary_inner == boundary_condition_boundary_layer ) {
		logging::print_master(LOG_INFO "Boundary Layer: Radial velocity at inner boundary is %e * V_Kepler.\n", vrad_fraction_of_kepler);
		logging::print_master(LOG_INFO "Boundary Layer: Stellar rotation rate is %f * Om_Kepler.\n", stellar_rotation_rate);
	}
	if( boundary_outer == boundary_condition_boundary_layer ) {
		logging::print_master(LOG_INFO "Boundary Layer: Mass Accretion Rate is %g Solar Masses per Year.\n", mass_accretion_rate*units::mass.get_cgs_factor()/units::time.get_cgs_factor()*units::cgs_Year/units::cgs_Msol);
	}
	logging::print_master(LOG_INFO "Boundary Layer: Radial Viscosity is multiplied by a factor of %f.\n", radial_viscosity_factor);

	if (damping)
	{
		for(unsigned int i = 0; i < damping_vector.size(); ++i)
		{
			logging::print_master(LOG_INFO "%s\n", damping_vector[i].description_inner.c_str());
			logging::print_master(LOG_INFO "%s\n", damping_vector[i].description_outer.c_str());
		}
	}
	else
	{
		logging::print_master(LOG_INFO "Damping at boundaries is disabled.\n");
	}

	logging::print_master(LOG_INFO "Surface density factor: %g\n", density_factor);
	logging::print_master(LOG_INFO "Tau factor: %g\n", tau_factor);
	logging::print_master(LOG_INFO "Kappa factor: %g\n", kappa_factor);

    logging::print_master(LOG_INFO "Minimum temperature: %.5e\n", minimum_temperature);
    logging::print_master(LOG_INFO "Maximum temperature: %.5e\n", maximum_temperature);

	logging::print_master(LOG_INFO "Heating from star is %s. Using %s model with ramping time of %g and a total factor %g.\n", heating_star_enabled ? "enabled" : "disabled", heating_star_simple ? "simplified" : "advanced", heating_star_ramping_time, heating_star_factor);
	logging::print_master(LOG_INFO "Heating from viscous dissipation is %s. Using a total factor of %g.\n", heating_viscous_enabled ? "enabled" : "disabled", heating_viscous_factor);
	logging::print_master(LOG_INFO "Cooling (beta) is %s. Using beta = %g.\n", cooling_beta_enabled ? "enabled" : "disabled", cooling_beta);
	logging::print_master(LOG_INFO "Cooling (radiative) is %s. Using a total factor of %g.\n", cooling_radiative_enabled ? "enabled" : "disabled", cooling_radiative_factor);
	logging::print_master(LOG_INFO "Radiative diffusion is %s. Using %s omega = %lf with a maximum %u interations.\n", radiative_diffusion_enabled ? "enabled" : "disabled", radiative_diffusion_omega_auto_enabled ? "auto" : "fixed", radiative_diffusion_omega, radiative_diffusion_max_iterations);

	logging::print_master(LOG_INFO "CFL parameter: %g\n", CFL);

	switch (opacity) {
		case opacity_lin:
			logging::print_master(LOG_INFO "Opacity uses tables from Lin & Papaloizou, 1985\n");
			break;
		case opacity_bell:
			logging::print_master(LOG_INFO "Opacity uses tables from Bell & Lin, 1994\n");
			break;
		case opacity_zhu:
			logging::print_master(LOG_INFO "Opacity uses tables from Zhu et al., 2012\n");
			break;
		case opacity_kramers:
			logging::print_master(LOG_INFO "Kramers opacity and constant electron scattering (Thomson) used.\n");
			break;
		case opacity_const_op:
			logging::print_master(LOG_INFO "Using constant opacity kappa_R = %e.\n", kappa_const);
			break;
	}

	if (write_lightcurves) {
		char *buffer, *temp;
		unsigned int pos = 0;

		buffer = (char*)malloc((pos+1)*sizeof(char));
		for (unsigned int i = 0; i < lightcurves_radii.size(); ++i) {
			int size = asprintf(&temp, "%lf, ", lightcurves_radii[i]);
			buffer = (char*)realloc(buffer, (pos+size+1)*sizeof(char));
			strcpy(&buffer[pos], temp);
			free(temp);
			pos+=size;
		}
		buffer[pos-2] = '0';

		logging::print_master(LOG_INFO "Lightcurves radii are: %s\n", buffer);
	}

	// particles
	logging::print_master(LOG_INFO "Particles are %s.\n", integrate_particles ? "enabled" : "disabled");
	if (integrate_particles) {
		logging::print_master(LOG_INFO "Using %u particles with a radius of %g and a density of %g.\n", number_of_particles, particle_radius, particle_density);
		logging::print_master(LOG_INFO "Distributing particles with a r^%.2g profile from %g to %g.\n", particle_slope, particle_minimum_radius, particle_maximum_radius);
		logging::print_master(LOG_INFO "Particles are considered escaped from the system when they reach a distance of %g.\n", particle_escape_radius);
		logging::print_master(LOG_INFO "Particles gas drag is %s.\n", particle_gas_drag_enabled ? "enabled" : "disabled");
		logging::print_master(LOG_INFO "Particles disk gravity is %s.\n", particle_disk_gravity_enabled ? "enabled" : "disabled");
	}
}

void write_grid_data_to_file()
{
  /* Write a file containing the base units to the output folder. */

  FILE *fd = 0;
  char *fd_filename;

  if (CPU_Master) {
    if (asprintf(&fd_filename, "%s%s", OUTPUTDIR, "dimensions.dat") == -1) {
      logging::print_master(LOG_ERROR "Not enough memory for string buffer.\n");
      PersonalExit(1);
    }
    fd = fopen(fd_filename, "w");
    if (fd == NULL) {
      logging::print_master(LOG_ERROR
                            "Can't write 'dimensions.dat' file. Aborting.\n");
      PersonalExit(1);
    }

    free(fd_filename);


    //fprintf(fd, "#XMIN\tXMAX\tYMIN\tYMAX\tNX\tNY\tNGHX\tNGHY\tRadial_spacing\n");

    char radial_spacing_str[256];

    switch(radial_grid_type)
    {
        case(arithmetic_spacing):
        {
            strncpy(radial_spacing_str, "Arithmetic", 256);
            break;
        }
        case(logarithmic_spacing):
        {
            strncpy(radial_spacing_str, "Logarithmic", 256);
            break;
        }
        case(exponential_spacing):
        {
            strncpy(radial_spacing_str, "Exponential", 256);
            break;
        }
    }

    fprintf(fd, "#RMIN\tRMAX\tPHIMIN\tPHIMAX          \tNRAD\tNAZ\tNGHRAD\tNGHAZ\tRadial_spacing\n");
    fprintf(fd, "%.16g\t%.16g\t%.16g\t%.16g\t%d\t%d\t%d\t%d\t%s\n",
            RMIN, RMAX, 0.0, 2*PI, NRadial, NAzimuthal, 1, 1, radial_spacing_str);
    fclose(fd);
  }
  MPI_Barrier(MPI_COMM_WORLD);

}

} // namespace parameters
