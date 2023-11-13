/**
	\file parameters.cpp
	\author Tobias Mueller <Tobias_Mueller@twam.info>
*/
#include <algorithm>
#include <cmath>

#include <vector>

#include "LowTasks.h"
#include "boundary_conditions/boundary_conditions.h"
#include "config.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "output.h"
#include "parameters.h"
#include "config.h"
#include "units.h"
#include "config.h"
#include "fld.h"

#include <limits>
constexpr double DBL_EPSILON = std::numeric_limits<double>::epsilon();

namespace parameters
{

double aspectratio_ref;
int aspectratio_mode;
double constant_viscosity;
double viscous_alpha;
bool VISCOUS_ACCRETION;
double sigma_slope;
double OMEGAFRAME;
double IMPOSEDDISKDRIFT;
double flaring_index;
double ADIABATICINDEX; // Also used for polytropic energy equation
double POLYTROPIC_CONSTANT;

bool CartesianParticles;
bool ParticleGravityCalcInCartesian;


double monitor_timestep;
unsigned int Nmonitor;
unsigned int Nsnap;
double quantities_radius_limit;
double disk_radius_mass_fraction;


double cps;

int ShockTube = false;
bool SpreadingRing = false;

// energy equations
bool Adiabatic = false;
bool Polytropic = false;
bool Locally_Isothermal = false;

bool variableGamma = false;
double hydrogenMassFraction;


t_radial_grid radial_grid_type;
const char *radial_grid_names[] = {"logarithmic", "arithmetic", "exponential",
				   "custom"};
double exponential_cell_size_factor;

double minimum_temperature;
double maximum_temperature;
double MU;
bool heating_viscous_enabled;
double heating_viscous_factor;
bool heating_star_enabled;

double surface_cooling_factor;
bool cooling_surface_enabled;
bool cooling_beta_enabled;
double cooling_beta_ramp_up;
double cooling_beta;
bool cooling_beta_reference;
bool cooling_beta_model;
bool cooling_beta_floor;

bool cooling_scurve_enabled;
bool cooling_scurve_type;

t_initialize_condition sigma_initialize_condition;
std::string sigma_filename;
int random_seed;
bool sigma_randomize;
double sigma_random_factor;
double sigma_floor;
double sigma_feature_size;
bool sigma_adjust;
double sigma_diskmass;
double sigma0;
t_initialize_condition energy_initialize_condition;
std::string energy_filename;

bool keep_mass_constant;

t_artificial_viscosity artificial_viscosity;
double artificial_viscosity_factor;
bool artificial_viscosity_dissipation;

bool calculate_disk;

// control centering of frame
unsigned int n_bodies_for_hydroframe_center;
unsigned int corotation_reference_body;
bool corotating;


int AlphaMode;
double alphaCold;
double alphaHot;

bool cbd_ring;
double cbd_ring_position;
double cbd_ring_width;
double cbd_decay_width;
double cbd_decay_exponent;
double cbd_ring_enhancement_factor;
double center_mass_density_correction_factor;

bool profile_cutoff_outer;
double profile_cutoff_point_outer;
double profile_cutoff_width_outer;

bool profile_cutoff_inner;
double profile_cutoff_point_inner;
double profile_cutoff_width_inner;

bool disk_feedback;
bool accrete_without_disk_feedback;
bool fast_transport;
int hydro_integrator;
int indirect_term_mode;

bool planet_orbit_disk_test;

bool do_init_secondary_disk;

double density_factor;
double tau_factor;
double tau_min;
double kappa_factor;

bool v_azimuthal_with_quadropole_support;

bool self_gravity;
t_sg self_gravity_mode;
unsigned int self_gravity_steps_between_kernel_update;
double self_gravity_aspectratio_change_threshold;

bool body_force_from_potential;

bool write_torques;

bool bitwise_exact_restarting;
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
bool compatibility_smoothing_planetloc;
bool compatibility_no_star_smoothing;
bool correct_disk_selfgravity;

bool initialize_pure_keplerian;
bool initialize_vradial_zero;

double star_radius;
double star_temperature;

double radial_viscosity_factor;
double mass_accretion_rate;
double accretion_radius_fraction;
double klahr_smoothing_radius;

double visc_accret_massflow_test;

double CFL;
double CFL_max_var;
double HEATING_COOLING_CFL_LIMIT;

unsigned int number_of_particles;
bool integrate_particles;
double particle_radius;
unsigned int particle_species_number;
double particle_radius_increase_factor;
double particle_eccentricity;
double particle_density;
double particle_slope;
double particle_minimum_radius;
double particle_maximum_radius;
double particle_minimum_escape_radius;
double particle_maximum_escape_radius;
double particle_minimum_escape_radius_sq;
double particle_maximum_escape_radius_sq;
bool particle_gas_drag_enabled;
bool particle_disk_gravity_enabled;
bool particle_dust_diffusion;
t_particle_integrator particle_integrator;

// for constant opacity
double kappa_const = 1.0;

void exitOnDeprecatedSetting(std::string setting_name, std::string reason,
			     std::string instruction)
{
    if (config::cfg.contains(setting_name)) {
	std::string warn = "Deprecated setting '" + setting_name +
			   "' found in config file. " + reason + " " +
			   instruction + "\n";
	logging::print_master((LOG_ERROR + warn).c_str());
	die("Config Error");
    }
}

static void read_output_config(t_data &data) {
    // output settings
    bool do_write_1D = config::cfg.get_flag("DoWrite1DFiles", true);
    data[t_data::SIGMA].set_write(
	config::cfg.get_flag("WriteDensity", true), do_write_1D);
    data[t_data::V_RADIAL].set_write(
	config::cfg.get_flag("WriteVelocity", true), do_write_1D);
    data[t_data::V_AZIMUTHAL].set_write(
	config::cfg.get_flag("WriteVelocity", true), do_write_1D);
    data[t_data::ENERGY].set_write(
	config::cfg.get_flag("WriteEnergy", true), do_write_1D);
    data[t_data::TEMPERATURE].set_write(
	config::cfg.get_flag("WriteTemperature", false), do_write_1D);
    data[t_data::SOUNDSPEED].set_write(
	config::cfg.get_flag("WriteSoundSpeed", false), do_write_1D);
    data[t_data::GAMMAEFF].set_write(
	config::cfg.get_flag("WriteEffectiveGamma", false),
	do_write_1D);
    data[t_data::GAMMA1].set_write(
	config::cfg.get_flag("WriteFirstAdiabaticIndex", false),
	do_write_1D);
	data[t_data::ALPHA].set_write(
	config::cfg.get_flag("WriteAlpha", false),
	do_write_1D);
    data[t_data::MU].set_write(
	config::cfg.get_flag("WriteMeanMolecularWeight", false),
	do_write_1D);
    data[t_data::PRESSURE].set_write(
	config::cfg.get_flag("WritePressure", false), do_write_1D);
    data[t_data::TOOMRE].set_write(
	config::cfg.get_flag("WriteToomre", false), do_write_1D);
    data[t_data::ECCENTRICITY_X].set_write(
	config::cfg.get_flag("WriteEccentricity", false), do_write_1D);
    data[t_data::ECCENTRICITY_Y].set_write(
	config::cfg.get_flag("WriteEccentricity", false), do_write_1D);
    data[t_data::POTENTIAL].set_write(
	config::cfg.get_flag("WritePotential", false), do_write_1D);
    data[t_data::QPLUS].set_write(
	config::cfg.get_flag("WriteQPlus", false), do_write_1D);
    data[t_data::QMINUS].set_write(
	config::cfg.get_flag("WriteQMinus", false), do_write_1D);
    data[t_data::KAPPA].set_write(
	config::cfg.get_flag("WriteKappa", false), do_write_1D);
    data[t_data::TAU_COOL].set_write(
	config::cfg.get_flag("WriteTauCool", false), do_write_1D);
    data[t_data::ALPHA_GRAV].set_write(
	config::cfg.get_flag("WriteAlphaGrav", false), do_write_1D);
    data[t_data::ALPHA_GRAV_MEAN].set_write(
	config::cfg.get_flag("WriteAlphaGravMean", false),
	do_write_1D);
    data[t_data::ALPHA_REYNOLDS].set_write(
	config::cfg.get_flag("WriteAlphaReynolds", false),
	do_write_1D);
    data[t_data::ALPHA_REYNOLDS_MEAN].set_write(
	config::cfg.get_flag("WriteAlphaReynoldsMean", false),
	do_write_1D);
    data[t_data::VISCOSITY].set_write(
	config::cfg.get_flag("WriteViscosity", false), do_write_1D);
    data[t_data::DIV_V].set_write(
	config::cfg.get_flag("WriteDivV", false), do_write_1D);
    data[t_data::T_REYNOLDS].set_write(
	config::cfg.get_flag("WriteTReynolds", false), do_write_1D);
    data[t_data::T_GRAVITATIONAL].set_write(
	config::cfg.get_flag("WriteTGravitational", false),
	do_write_1D);
    data[t_data::ADVECTION_TORQUE].set_write(
	config::cfg.get_flag("WriteGasTorques", false), do_write_1D);
    data[t_data::GRAVITATIONAL_TORQUE_NOT_INTEGRATED].set_write(
	config::cfg.get_flag("WriteGasTorques", false), do_write_1D);
    data[t_data::VISCOUS_TORQUE].set_write(
	config::cfg.get_flag("WriteGasTorques", false), do_write_1D);
    data[t_data::P_DIVV].set_write(
	config::cfg.get_flag("WritepDV", false), do_write_1D);
    data[t_data::TAU].set_write(
	config::cfg.get_flag("WriteTau", false), do_write_1D);
    data[t_data::SCALE_HEIGHT].set_write(
	config::cfg.get_flag("WriteScaleHeight", false), do_write_1D);
    data[t_data::ASPECTRATIO].set_write(
	config::cfg.get_flag("WriteAspectratio", false), do_write_1D);
    data[t_data::VISIBILITY].set_write(
	config::cfg.get_flag("WriteVisibility", false), do_write_1D);
    data[t_data::LUMINOSITY_1D].set_write(
	config::cfg.get_flag("WriteRadialLuminosity", false));
    data[t_data::DISSIPATION_1D].set_write(
	config::cfg.get_flag("WriteRadialDissipation", false));
    data[t_data::TAU_EFF].set_write(
	config::cfg.get_flag("WriteVerticalOpticalDepth", false),
	do_write_1D);
	data[t_data::SG_ACCEL_RAD].set_write(
	config::cfg.get_flag("WriteSGAccelRad", false),do_write_1D);
	data[t_data::SG_ACCEL_AZI].set_write(
	config::cfg.get_flag("WriteSGAccelAzi", false),do_write_1D);

    write_torques = config::cfg.get_flag("WriteTorques", false);

    write_disk_quantities =
	config::cfg.get_flag("WriteDiskQuantities", true);
    write_at_every_timestep =
	config::cfg.get_flag("WriteAtEveryTimestep", true);
    write_lightcurves =
	config::cfg.get_flag("WriteLightCurves", false);

	bitwise_exact_restarting = config::cfg.get_flag("BitwiseExactRestarting", false);

    write_massflow = config::cfg.get_flag("WriteMassFlow", false);
    data[t_data::MASSFLOW].set_write(write_massflow, do_write_1D);

    log_after_steps = config::cfg.get<unsigned int>("LogAfterSteps", 0);
    log_after_real_seconds =
	config::cfg.get<double>("LogAfterRealSeconds", 600.0);


    // parse light curve radii
    if (config::cfg.contains("WriteLightCurvesRadii")) {
	// get light curves radii string

	std::string lightcurve_config =
	    config::cfg.get<std::string>("WriteLightCurvesRadii");

	char *lightcurves_radii_string =
	    new char[lightcurve_config.length() + 1];
	strcpy(lightcurves_radii_string, lightcurve_config.c_str());

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

	delete[] lightcurves_radii_string;

	// sort vector
	sort(lightcurves_radii.begin(), lightcurves_radii.end());
    }

}


static void read_scurve_config() {

	const std::string str = config::cfg.get_lowercase("ScurveType", "Kimura");

	if (str == "kimura") { // Kimura et al. 2020 (https://doi.org/10.1093/pasj/psz144)
		cooling_scurve_type = true;
	} else if (str == "ichikawa") { // Ichikawa & Osaki 1992 (https://ui.adsabs.harvard.edu/abs/1992PASJ...44...15I/abstract)
		cooling_scurve_type = false;
	} else {
		throw std::runtime_error("Invalid choice for scurve type: " + str);
	}
}


static void read_surface_cooling_config(){

	cooling_surface_enabled = false;
    cooling_scurve_enabled = false;
	
	const std::string surface_cooling = config::cfg.get_lowercase("SurfaceCooling", "No");
	if (surface_cooling == "no" || surface_cooling == "off" || surface_cooling == "false") {
		cooling_surface_enabled = false;
		cooling_scurve_enabled = false;
	} else if (surface_cooling == "thermal") {
		cooling_surface_enabled = true;
		cooling_scurve_enabled = false;
	} else if (surface_cooling == "scurve") {
		cooling_surface_enabled = false;
		cooling_scurve_enabled = true;
	} else {
		throw std::runtime_error("Invalid choice for surface cooling: " + surface_cooling);
	}

	read_scurve_config();

    surface_cooling_factor =
	config::cfg.get<double>("CoolingRadiativeFactor", 1.0);
}

static void read_opacity_config(){

	const std::string str = config::cfg.get_lowercase("Opacity", "Lin");

	if (str == "lin") {
		opacity = opacity_lin;
	} else if (str == "bell") {
		opacity = opacity_bell;
	} else if (str == "constant") {
		opacity = opacity_const_op;
	} else if (str == "simple") {
		opacity = opacity_simple; // see Gennaro D'Angelo et al. 2003
	} else {
		throw std::runtime_error("Invalid choice for opacity type: " + str);
	}

	// Get this always for setting default value.
    units::precise_unit L0 = units::L0;
    units::precise_unit M0 = units::M0;

	kappa_const = config::cfg.get<double>("KappaConst", 1.0, (L0*L0)/M0);
}


static void read_beta_cooling_config(){
	cooling_beta_enabled = config::cfg.get_flag("CoolingBetaLocal", "No");

    cooling_beta = config::cfg.get<double>("CoolingBeta", 1.0);

    units::precise_unit T0 = units::T0;
    cooling_beta_ramp_up =
	config::cfg.get<double>("CoolingBetaRampUp", 0.0, T0);

	cooling_beta_reference = false;
	cooling_beta_model = false;
	cooling_beta_floor = false;

	const std::string str = config::cfg.get_lowercase("CoolingBetaReference", "Zero");
	if (str == "zero") {
	} else if (str == "reference") {
		cooling_beta_reference = true;
	} else if (str == "diskmodel") {
		cooling_beta_model = true;
	} else if (str == "floor") {
		cooling_beta_floor = true;
	} else {
		throw std::runtime_error("Invalid choice for cooling beta reference: " + str);
	}

}

static void read_radiation_config() {

	// read all parameters to get default values
	read_surface_cooling_config();
	read_opacity_config();
	read_beta_cooling_config();
	fld::config();
}

    

void read(const std::string &filename, t_data &data)
{
	if (filename.compare(filename.size()-4,4,".par") == 0) {
		die("This version of fargo uses the new yaml config files.\n%s looks like a par file.\nUse Tools/ini2yml.py to convert your config file!", filename.c_str());
	}
	config::cfg.load_file(filename);
	logging::print_master(LOG_INFO "Using parameter file %s\n", filename.c_str());
    // units
	const std::string l0s = config::cfg.get<std::string>("l0", "1.0");
	const std::string m0s = config::cfg.get<std::string>("m0", "1.0");
	const std::string t0s = config::cfg.get<std::string>("t0", "1.0");
	const std::string temp0s = config::cfg.get<std::string>("temp0", "1.0");


	units::set_baseunits(l0s, m0s, t0s, temp0s);
	// units::set_baseunits(l0s, m0s);

	units::precise_unit L0 = units::L0;
	units::precise_unit M0 = units::M0;
	units::precise_unit T0 = units::T0;
	units::precise_unit Temp0 = units::Temp0;

	constants::initialize_constants();

    // now we now everything to compute unit factors
    units::calculate_unit_factors();

    /* grid */
    RMIN = config::cfg.get<double>("Rmin", L0);
    RMAX = config::cfg.get<double>("Rmax", L0);

    NRadial = config::cfg.get<unsigned int>("Nrad", 64);
    NAzimuthal = config::cfg.get<unsigned int>("Naz", 64);


    // Disk radius = radius at which disk_radius_mass_fraction percent
    // of the total mass inside the domain is contained
    // 0.99 is used in Kley et al. 2008 "Simulations of eccentric disks .."
    disk_radius_mass_fraction =
    config::cfg.get<double>("DiskRadiusMassFraction", 0.99);

    quantities_radius_limit =
	config::cfg.get<double>("QuantitiesRadiusLimit", 2.0 * RMAX, L0);

    if (quantities_radius_limit <= RMIN) {
	logging::print_master(LOG_INFO "QuantitiesRadiusLimit %.5e < %.5e (Rmin) too small, setting it to %.5e\n", quantities_radius_limit, RMIN, 2.0*RMAX);
	quantities_radius_limit = 2.0 * RMAX;
    }
	logging::print_master(LOG_INFO "Computing disk quantities within %.5e L0 from coordinate center\n", quantities_radius_limit);

    exponential_cell_size_factor =
	config::cfg.get<double>("ExponentialCellSizeFactor", 1.41);
    switch (config::cfg.get_first_letter_lowercase("RadialSpacing", "ARITHMETIC")) {
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

	read_output_config(data);

    // boundary conditions
	boundary_conditions::parse_config();

    calculate_disk = config::cfg.get_flag("Disk", "yes");
	
    corotation_reference_body =
	config::cfg.get<unsigned int>("CorotationReferenceBody", 1);

    // set number of bodies used to calculate barycenter in Interpret.cpp

    MU = config::cfg.get<double>("mu", 1.0);

    minimum_temperature =
	config::cfg.get<double>("MinimumTemperature", "3 K", Temp0);
    maximum_temperature =
	config::cfg.get<double>("MaximumTemperature", "1.0e300 K", Temp0);

    heating_viscous_enabled =
	config::cfg.get_flag("HeatingViscous", "Yes");
    heating_viscous_factor =
	config::cfg.get<double>("HeatingViscousFactor", 1.0);

	read_radiation_config();

    // initialisation
    initialize_pure_keplerian =
	config::cfg.get_flag("InitializePureKeplerian", "no");
    initialize_vradial_zero =
	config::cfg.get_flag("InitializeVradialZero", "no");

    switch (config::cfg.get_first_letter_lowercase("SigmaCondition", "Profile")) {
    case 'p': // Profile
	sigma_initialize_condition = initialize_condition_profile;
	break;
    case 'n': // Profile
	sigma_initialize_condition =
	    initialize_condition_profile_Nbody_centered;
	energy_initialize_condition =
	    initialize_condition_profile_Nbody_centered;
	break;
    case '1': // 1D
	sigma_initialize_condition = initialize_condition_read1D;
	break;
    case '2': // 2D
	sigma_initialize_condition = initialize_condition_read2D;
	break;
    default:
	die("Invalid setting for SigmaCondition: %s",
	    config::cfg.get<std::string>("SigmaCondition", "Profile").c_str());
    }

	sigma_filename = config::cfg.get<std::string>("SigmaFilename", "");
    
    switch (config::cfg.get_first_letter_lowercase("EnergyCondition", "Profile")) {
    case 'p': // Profile
	energy_initialize_condition = initialize_condition_profile;
	break;
    case 'n': // Profile
	sigma_initialize_condition =
	    initialize_condition_profile_Nbody_centered;
	energy_initialize_condition =
	    initialize_condition_profile_Nbody_centered;
	break;
    case '1': // 1D
	energy_initialize_condition = initialize_condition_read1D;
	break;
    case '2': // 2D
	energy_initialize_condition = initialize_condition_read2D;
	break;
    default:
	die("Invalid setting for EnergyCondition: %s",
	    config::cfg.get<std::string>("EnergyCondition", "Profile").c_str());
    }

	energy_filename = config::cfg.get<std::string>("EnergyFilename", "");

    random_seed = config::cfg.get<int>("RandomSeed", 0);
    sigma_randomize = config::cfg.get_flag("RandomSigma", "no");
    sigma_random_factor = config::cfg.get<double>("RandomFactor", 0.1);
    sigma_feature_size =
	config::cfg.get<double>("FeatureSize", (RMAX - RMIN) / 150, L0);
    sigma_floor = config::cfg.get<double>("SigmaFloor", 1e-9);
    sigma0 = config::cfg.get<double>("Sigma0", "173 g/cm2", M0/(L0*L0));
    sigma_adjust = config::cfg.get_flag("SetSigma0", "no");
    sigma_diskmass = config::cfg.get<double>("DiskMass", 0.01, M0);
    density_factor = config::cfg.get<double>("DensityFactor", std::sqrt(2.0 * M_PI));

    tau_factor = config::cfg.get<double>("TauFactor", 0.5);
    kappa_factor = config::cfg.get<double>("KappaFactor", 1.0);
	tau_min = config::cfg.get<double>("TauMin", 0.01);

    v_azimuthal_with_quadropole_support = config::cfg.get_flag("VazimuthalConsidersQuadropoleMoment", "no");

    // artificial visocisty
    switch (config::cfg.get_first_letter_lowercase("ArtificialViscosity", "SN")) {
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
	die("Invalid setting for ArtificialViscosity: %s",
	    config::cfg.get<std::string>("ArtificialViscosity", "SN").c_str());
    }
    artificial_viscosity_dissipation =
	config::cfg.get_flag("ArtificialViscosityDissipation", "yes");
    artificial_viscosity_factor =
	config::cfg.get<double>("ArtificialViscosityFactor", 1.41);
    // warning
    if (config::cfg.contains("CVNR")) {
	die("Parameter CVNR has been renamed to ArtificialViscosityFactor");
    }

	keep_mass_constant =
	config::cfg.get_flag("KeepDiskMassConstant", "no");
	if(keep_mass_constant){
		logging::print_master(LOG_INFO "Disk mass is kept constant at initial value.\n");
	}


    // self gravity
    self_gravity = config::cfg.get_flag("SelfGravity", "no");

	const std::string sgmode = config::cfg.get_lowercase("SelfGravityMode", "besselkernel");
	if (sgmode == "basic") {
		self_gravity_mode = t_sg::sg_B;
	} else if (sgmode == "symmetric") {
		self_gravity_mode = t_sg::sg_S;
	} else if (sgmode == "besselkernel") {
		self_gravity_mode = t_sg::sg_BK;
	} else {
		logging::print_master(LOG_ERROR "Selfgravity mode %s is not supported. Choices: b, s, BesselKernel\n", sgmode.c_str());
		die("Configuration error.");
	}

	self_gravity_steps_between_kernel_update = config::cfg.get<unsigned int>("SelfGravityStepsBetweenKernelUpdate", 20);
	self_gravity_aspectratio_change_threshold = config::cfg.get<double>("SelfGravityAspectRatioChangeThreshold", 0.001);

    if (self_gravity) {
		logging::print_master(LOG_INFO "Self gravity enabled. It uses the '%s' mode. The kernel is updated every %u steps and after aspect ratio changed by %f.\n", sgmode.c_str(), self_gravity_steps_between_kernel_update, self_gravity_aspectratio_change_threshold);
    }



    //
    thickness_smoothing =
	config::cfg.get<double>("ThicknessSmoothing", 0.6);
    thickness_smoothing_sg = config::cfg.get<double>(
	"ThicknessSmoothingSG", 1.2); // recommended value from Müller, Kley & Meru 2012
	compatibility_smoothing_planetloc = config::cfg.get_flag("CompatibilitySmoothingPlanetLoc", "no");
	compatibility_no_star_smoothing = config::cfg.get_flag("CompatibilityNoStarSmoothing", "no");

	correct_disk_selfgravity = config::cfg.get_flag("CorrectDiskSelfgravity", self_gravity ? "no" : "yes");
    do_init_secondary_disk = config::cfg.get_flag("SecondaryDisk", "no");


	//local alpha
	AlphaMode = config::cfg.get<int>("AlphaMode", 0);
	alphaCold = config::cfg.get<double>("AlphaCold", 0.01);
	alphaHot = config::cfg.get<double>("AlphaHot", 0.1);

	if(parameters::AlphaMode == SCURVE_ALPHA || parameters::AlphaMode == SCURVE_IONFRACTION){
        // already continously writes alpha
        data[t_data::ALPHA].set_do_before_write(nullptr);
    }


	cbd_ring =
	config::cfg.get_flag("CircumBinaryRing", "no");
	cbd_ring_position =
	config::cfg.get<double>("CircumBinaryRingPosition", 4.5, L0);
	cbd_ring_width =
	config::cfg.get<double>("CircumBinaryRingWidth", 0.6, L0);
	cbd_decay_width =
	config::cfg.get<double>("CircumBinaryDecayWidth", cbd_ring_width*1.4, L0);
	cbd_decay_exponent =
	config::cfg.get<double>("CircumBinaryDecayExponent", 0.75);
	cbd_ring_enhancement_factor =
	config::cfg.get<double>("CircumBinaryRingEnhancementFactor", 2.5);
	center_mass_density_correction_factor =
	config::cfg.get<double>("CenterProfileDensityCorrectionFactor", 1.0);

    // profile damping outer
    profile_cutoff_outer =
	config::cfg.get_flag("ProfileCutoffOuter", "no");
    profile_cutoff_point_outer =
	config::cfg.get<double>("ProfileCutoffPointOuter", 1.0e300, L0);
    profile_cutoff_width_outer =
	config::cfg.get<double>("ProfileCutoffWidthOuter", 1.0, L0);

    // profile damping inner
    profile_cutoff_inner =
	config::cfg.get_flag("ProfileCutoffInner", "no");
    profile_cutoff_point_inner =
	config::cfg.get<double>("ProfileCutoffPointInner", 0.0, L0);
    profile_cutoff_width_inner =
	config::cfg.get<double>("ProfileCutoffWidthInner", 1.0, L0);

    exitOnDeprecatedSetting(
	"FeelsDisk",
	"Replaced by parameter DiskFeedback for clarification since it affects more than just the star.",
	"Please replace 'FeelsDisk' by 'DiskFeedback'.");

	exitOnDeprecatedSetting(
	"DefaultStar",
	"The star parameters are now part of the planet system configuration.",
	"Please add the central object to the planet configuration.");

    disk_feedback = config::cfg.get_flag("DiskFeedback", "yes");
	accrete_without_disk_feedback = config::cfg.get_flag("AccreteWithoutDiskFeedback", "no");
	planet_orbit_disk_test = config::cfg.get_flag("PlanetOrbitDiskTest", "no");

	indirect_term_mode = config::cfg.get<int>("IndirectTermMode", INDIRECT_TERM_REBOUND);

	switch(indirect_term_mode){
		case INDIRECT_TERM_REBOUND:
		{
			logging::print_master(LOG_INFO "Indirect Term computed as effective Hydro center acceleratrion with shifting the Nbody system to the center.\n");
			break;
		}
		case INDIRECT_TERM_EULER:
		{
			logging::print_master(LOG_INFO "Indirect Term computed as current force (euler) on Hydro center with shifting the Nbody system to the center.\n");
			break;
		}
	}


    body_force_from_potential =
	config::cfg.get_flag("BodyForceFromPotential", "yes");
    if (body_force_from_potential) {
	logging::print_master(LOG_INFO
			      "Body force on gas computed via potential.\n");
    } else {
	logging::print_master(LOG_INFO
			      "Body force on gas computed via force.\n");
    }

    // boundary layer parameters
    radial_viscosity_factor =
	config::cfg.get<double>("RadialViscosityFactor", 1.);

    mass_accretion_rate =
	config::cfg.get<double>("MassAccretionRate", 1.e-9, M0/T0);
    accretion_radius_fraction =
	config::cfg.get<double>("MassAccretionRadius", 1.0);
    klahr_smoothing_radius = config::cfg.get<double>(
	"KlahrSmoothingRadius", accretion_radius_fraction);

    visc_accret_massflow_test =
        config::cfg.get_flag("ViscAccretMassflowTest", "no");

    CFL = config::cfg.get<double>("CFL", 0.5);
    HEATING_COOLING_CFL_LIMIT =
	config::cfg.get<double>("HeatingCoolingCFLlimit", 10.0);

    CFL_max_var = config::cfg.get<double>("CFLmaxVar", 1.1);

    // particles
    CartesianParticles =
	config::cfg.get_flag("CartesianParticles", "no");
    integrate_particles =
	config::cfg.get_flag("IntegrateParticles", "no");
    number_of_particles =
	config::cfg.get<unsigned int>("NumberOfParticles", 0);

    particle_radius = config::cfg.get<double>("ParticleRadius", "100.0 cm", L0);
    particle_species_number =
	config::cfg.get<unsigned int>("ParticleSpeciesNumber", 1);
    particle_radius_increase_factor =
	config::cfg.get<double>("ParticleRadiusIncreaseFactor", 10.0);
    particle_eccentricity =
	config::cfg.get<double>("ParticleEccentricity", 0.0);
    particle_density = config::cfg.get<double>("ParticleDensity", "2.65 g/cm3", M0/(L0*L0*L0));
    particle_slope = config::cfg.get<double>(
	"ParticleSurfaceDensitySlope", sigma_slope);
    particle_slope =
	-particle_slope; // particle distribution scales with  r^slope, so we
			 // introduces the minus here to make it r^-slope (same
			 // as for gas)
    particle_slope +=
	1.0; // particles are distributed over a whole simulation ring which
	     // introduces a factor 1/r for the particle surface density
    particle_minimum_radius =
	config::cfg.get<double>("ParticleMinimumRadius", RMIN, L0);
    particle_maximum_radius =
	config::cfg.get<double>("ParticleMaximumRadius", RMAX, L0);
    particle_minimum_escape_radius =
	config::cfg.get<double>("ParticleMinimumEscapeRadius", RMIN, L0);
    particle_maximum_escape_radius =
	config::cfg.get<double>("ParticleMaximumEscapeRadius", RMAX, L0);
    particle_gas_drag_enabled =
	config::cfg.get_flag("ParticleGasDragEnabled", true);
    particle_disk_gravity_enabled =
	config::cfg.get_flag("ParticleDiskGravityEnabled", false);
	particle_dust_diffusion =
	config::cfg.get_flag("ParticleDustDiffusion", false);



    // particle integrator
    switch (
	config::cfg.get_first_letter_lowercase("ParticleIntegrator", "m")) {
    case 'e': // Adaptive explicit
	particle_integrator = integrator_adaptive;
	break;
    case 'm': // exponential midpoint
	particle_integrator = integrator_exponential_midpoint;

	if (!particle_gas_drag_enabled) {
	    logging::print_master(
		LOG_ERROR
		"Do not use exponential midpoint particle integrator without gas drag, use the explicit integrator instead.\n");
	}

	break;
    default:
    die("Invalid setting for Particle Integrator: %s	with key %c",
	    config::cfg.get<std::string>("ParticleIntegrator", "s").c_str(),
	    config::cfg.get_first_letter_lowercase("ParticleIntegrator", "s"));
    }

    if (CartesianParticles && (particle_integrator == integrator_exponential_midpoint)) {
	// exponential midpoint integrator only implemented in polar
	// coordiantes, but forces can be calculated in cartesian coordinates
	CartesianParticles = false;
	ParticleGravityCalcInCartesian = true;
    }

    if (particle_disk_gravity_enabled && (!self_gravity)) {
	logging::print_master(
	    LOG_ERROR
	    "Cannot enable particle_disk_gravity_enabled while self_gravity is off!\n");
	PersonalExit(1);
    }

    // second and second last radius, so that partices stay out of the ghost
    // cells
    if (particle_minimum_escape_radius < RMIN) {
	logging::print_master(
	    LOG_WARNING
	    "particle_minimum_escape_radius can't be smaller than the inner radius of the domain. Setting particle_minimum_escape_radius to inner radius of the domain.\n");
	particle_minimum_escape_radius = RMIN;
    }

    if (particle_maximum_escape_radius > RMAX) {
	logging::print_master(
	    LOG_WARNING
	    "particle_maximum_escape_radius can't be larger than the outer radius of the domain. Setting particle_maximum_escape_radius to outer radius of the domain.\n");
	particle_maximum_escape_radius = RMAX;
    }

    particle_maximum_escape_radius_sq =
	std::pow(particle_maximum_escape_radius, 2) -
	DBL_EPSILON; // DBL for safety
    particle_minimum_escape_radius_sq =
	std::pow(particle_minimum_escape_radius, 2) + DBL_EPSILON;
}


void summarize_parameters()
{
    // artifical viscosity
    switch (artificial_viscosity) {
    case artificial_viscosity_none:
	logging::print_master(LOG_INFO "Using no artificial viscosity.\n");
	break;
    case artificial_viscosity_TW:
	logging::print_master(
	    LOG_INFO
	    "Using Tscharnuter-Winkler (1979) artificial viscosity with C = %lf.\n",
	    parameters::artificial_viscosity_factor);
	logging::print_master(
		LOG_INFO "Artificial viscosity is %s for dissipation.\n",
		parameters::artificial_viscosity_dissipation ? "used" : "not used");
		break;
    case artificial_viscosity_SN:
	logging::print_master(
	    LOG_INFO
	    "Using Stone-Norman (1991, ZEUS-2D) artificial viscosity with C = %lf.\n",
	    parameters::artificial_viscosity_factor);
	logging::print_master(
	    LOG_INFO "Artificial viscosity is %s for dissipation.\n",
	    parameters::artificial_viscosity_dissipation ? "used" : "not used");
	break;
    }

    // Mass Transfer
    if (boundary_conditions::rochelobe_overflow) {
	logging::print_master(
	    LOG_INFO
	    "Mass Transfer from planet #%d of %g M_sun/orbit with Ts = %g K and ramping time t_ramp = %g P_orb.\n",
	    boundary_conditions::rof_planet, boundary_conditions::rof_mdot, boundary_conditions::rof_temperature, boundary_conditions::rof_rampingtime);
    }
	if (boundary_conditions::rof_variableTransfer){
	logging::print_master(LOG_INFO
	"Mass Transfer is variable with gamma= %g and an averaging time of t_avg = %g P_orb. \n", boundary_conditions::rof_gamma, boundary_conditions::rof_averaging_time);
	}

	boundary_conditions::describe_damping();


    logging::print_master(LOG_INFO "Surface density factor: %g\n",
			  density_factor);
    logging::print_master(LOG_INFO "Tau factor: %g\n", tau_factor);
	logging::print_master(LOG_INFO "Tau min: %g\n", tau_min);
    logging::print_master(LOG_INFO "Kappa factor: %g\n", kappa_factor);

    logging::print_master(LOG_INFO "Minimum temperature: %.5e K = %.5e\n",
			  minimum_temperature, minimum_temperature*units::temperature);
    logging::print_master(LOG_INFO "Maximum temperature: %.5e K = %.5e\n",
			  maximum_temperature, maximum_temperature*units::temperature);

    logging::print_master(
	LOG_INFO
	"Heating from viscous dissipation is %s. Using a total factor of %g.\n",
	heating_viscous_enabled ? "enabled" : "disabled",
	heating_viscous_factor);
	logging::print_master(LOG_INFO "Cooling (beta) is %s and reference temperature is %s. Using beta = %g.\n",
			  cooling_beta_enabled ? "enabled" : "disabled",
			  cooling_beta_reference ? "the initial value" : cooling_beta_floor ? "floor" : "zero",
			  cooling_beta);

    logging::print_master(
	LOG_INFO "Cooling (radiative) is %s. Using a total factor of %g.\n",
	cooling_surface_enabled ? "enabled" : "disabled",
	surface_cooling_factor);

    logging::print_master(
        LOG_INFO "S-curve cooling is %s. \n",
        cooling_scurve_enabled ? "enabled" : "disabled");


    if (Adiabatic) {
	logging::print_master(
	    LOG_INFO
	    "CFL parameter: %g	heating/cooling (dT/T) limited to %g%% per hydro step\n",
	    CFL, 100.0 * HEATING_COOLING_CFL_LIMIT);
    } else {
	logging::print_master(LOG_INFO "CFL parameter: %g\n", CFL);
    }

    switch (opacity) {
    case opacity_lin:
	logging::print_master(
	    LOG_INFO "Opacity uses tables from Lin & Papaloizou, 1985\n");
	break;
    case opacity_bell:
	logging::print_master(LOG_INFO
			      "Opacity uses tables from Bell & Lin, 1994\n");
	break;
    case opacity_const_op:
	logging::print_master(LOG_INFO "Using constant opacity kappa_R = %e.\n",
			      kappa_const);
	break;
    case opacity_simple:
	logging::print_master(
	    LOG_INFO
	    "Using opacity from Gennaro D'Angelo et al. 2003 with kappa_0 = %e.\n",
	    kappa_const);
	break;
    }

    if (write_lightcurves) {

	std::string lightcurves_radii_string = "";
	for (unsigned int i = 0; i < lightcurves_radii.size(); ++i) {
		lightcurves_radii_string += std::to_string(lightcurves_radii[i]);
		if (i != lightcurves_radii.size()-1) {
			lightcurves_radii_string += ", ";
		}
	}
	logging::print_master(LOG_INFO "Lightcurves radii are: %s\n", lightcurves_radii_string.c_str());
    }

    // particles
    logging::print_master(LOG_INFO "Particles are %s.\n",
			  integrate_particles ? "enabled" : "disabled");
    if (integrate_particles) {
	logging::print_master(
	    LOG_INFO
	    "Using %u particles with a radius of %g and a density of %g.\n",
	    number_of_particles, particle_radius, particle_density);
	logging::print_master(
	    LOG_INFO
	    "Distributing particles with a r^%.2g profile from %g to %g with a eccentricity from 0.0 to %g.\n",
	    particle_slope, particle_minimum_radius, particle_maximum_radius,
	    particle_eccentricity);
	logging::print_master(
	    LOG_INFO
	    "Particles are considered escaped from the system when they reach a distance of %g or %g.\n",
	    particle_minimum_escape_radius, particle_maximum_escape_radius);
	logging::print_master(LOG_INFO "Particles gas drag is %s.\n",
			      particle_gas_drag_enabled ? "enabled"
							: "disabled");
	logging::print_master(LOG_INFO "Particles disk gravity is %s.\n",
			      particle_disk_gravity_enabled ? "enabled"
							    : "disabled");
	switch (particle_integrator) {
	case integrator_adaptive:
	    logging::print_master(
		LOG_INFO "Particles use the (explicit) adaptive integrator\n");
	    break;
	case integrator_exponential_midpoint:
	    logging::print_master(
		LOG_INFO "Particles use the exponential midpoint integrator\n");
	    break;
	default:
	    die("Invalid setting for Particle Integrator: %s",
		(config::cfg.get<std::string>("ParticleIntegrator", "s")).c_str());
	}
    }
}

void write_grid_data_to_file()
{
    /* Write a file containing the base units to the output folder. */

    FILE *fd = 0;
    
	if (CPU_Master) {

	const std::string filename = output::outdir + "dimensions.dat";

	fd = fopen(filename.c_str(), "w");
	if (fd == NULL) {
	    logging::print_master(
		LOG_ERROR "Can't write 'dimensions.dat' file. Aborting.\n");
	    PersonalExit(1);
	}

	// fprintf(fd,
	// "#XMIN\tXMAX\tYMIN\tYMAX\tNX\tNY\tNGHX\tNGHY\tRadial_spacing\n");

	char radial_spacing_str[256];

	switch (radial_grid_type) {
	case (arithmetic_spacing): {
	    strncpy(radial_spacing_str, "Arithmetic", 256);
	    break;
	}
	case (logarithmic_spacing): {
	    strncpy(radial_spacing_str, "Logarithmic", 256);
	    break;
	}
	case (exponential_spacing): {
	    strncpy(radial_spacing_str, "Exponential", 256);
	    break;
	}
	case (custom_spacing): {
	    strncpy(radial_spacing_str, "Custom", 256);
	    break;
	}
	}

	fprintf(
	    fd,
	    "#RMIN\tRMAX\tPHIMIN\tPHIMAX          \tNRAD\tNAZ\tNGHRAD\tNGHAZ\tRadial_spacing\n");
	fprintf(fd, "%.16g\t%.16g\t%.16g\t%.16g\t%d\t%d\t%d\t%d\t%s\n", RMIN,
		RMAX, 0.0, 2 * M_PI, NRadial, NAzimuthal, 1, 1,
		radial_spacing_str);
	fclose(fd);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace parameters
