/**
	\file parameters.cpp
	\author Tobias Mueller <Tobias_Mueller@twam.info>
*/
#include <algorithm>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <vector>

#include "LowTasks.h"
#include "boundary_conditions.h"
#include "config.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "output.h"
#include "parameters.h"
#include "string.h"
#include "util.h"
#include "viscosity.h"
#include "config.h"
#include "units.h"

#include <limits>
constexpr double DBL_EPSILON = std::numeric_limits<double>::epsilon();

namespace parameters
{

double ASPECTRATIO_REF;
int ASPECTRATIO_MODE;
int EXPLICIT_VISCOSITY;
double STS_NU;
double VISCOSITY;
double ALPHAVISCOSITY;
int VISCOUS_ACCRETION;
double SIGMASLOPE;
double OMEGAFRAME;
double IMPOSEDDISKDRIFT;
double FLARINGINDEX;
double ADIABATICINDEX; // Also used for polytropic energy equation
double POLYTROPIC_CONSTANT;

bool CartesianParticles;
bool ParticlesInCartesian;


double DT;
unsigned int NINTERM;
unsigned int NTOT;
double quantities_radius_limit;


int ShockTube = false;
bool SpreadingRing = false;

// energy equations
bool Adiabatic = false;
bool Polytropic = false;
bool Locally_Isothermal = false;

bool variableGamma = false;

t_radial_grid radial_grid_type;
const char *radial_grid_names[] = {"logarithmic", "arithmetic", "exponential",
				   "custom"};
double exponential_cell_size_factor;

std::string PRESCRIBED_BOUNDARY_OUTER_FILE;

t_boundary_condition boundary_inner;
t_boundary_condition boundary_outer;
bool domegadr_zero;

double viscous_outflow_speed;

bool damping;
bool is_damping_initial = false;
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
bool heating_star_simple;

double cooling_radiative_factor;
bool cooling_radiative_enabled;
bool cooling_beta_enabled;
double cooling_beta_ramp_up;
double cooling_beta;
bool cooling_beta_initial;
bool cooling_beta_aspect_ratio;


bool radiative_diffusion_enabled;
double radiative_diffusion_omega;
bool radiative_diffusion_omega_auto_enabled;
unsigned int radiative_diffusion_max_iterations;

t_initialize_condition sigma_initialize_condition;
std::string sigma_filename;
int random_seed;
bool sigma_randomize;
double sigma_random_factor;
double sigma_floor;
double sigma_feature_size;
bool sigma_adjust;
bool sigma0_in_code_units;
double sigma_discmass;
double sigma0;
t_initialize_condition energy_initialize_condition;
std::string energy_filename;

t_artificial_viscosity artificial_viscosity;
double artificial_viscosity_factor;
bool artificial_viscosity_dissipation;

bool calculate_disk;

// control centering of frame
unsigned int n_bodies_for_hydroframe_center;
unsigned int corotation_reference_body;

bool massoverflow;
unsigned int mof_planet;
double mof_temperature;
double mof_value;
double mof_rampingtime;

bool profile_cutoff_outer;
double profile_cutoff_point_outer;
double profile_cutoff_width_outer;

bool profile_cutoff_inner;
double profile_cutoff_point_inner;
double profile_cutoff_width_inner;

bool disk_feedback;

bool integrate_planets;
bool do_init_secondary_disk;

double density_factor;
double tau_factor;
double kappa_factor;

bool self_gravity;
bool body_force_from_potential;

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
double accretion_radius_fraction;
double klahr_smoothing_radius;

unsigned int zbuffer_size;
double zbuffer_maxangle;

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
t_particle_integrator integrator;

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

static t_DampingType write_damping_type(t_damping_type type_inner,
					t_damping_type type_outer,
					t_data::t_polargrid_type quantity,
					t_data::t_polargrid_type quantity0,
					std::string description)
{
    t_DampingType damping_type;
    damping_type.array_to_damp = quantity;
    damping_type.array_with_damping_values = quantity0;
    damping_type.type_inner = type_inner;
    damping_type.type_outer = type_outer;
    std::string description_inner;
    std::string description_outer;

    switch (type_inner) {
    case damping_none:
	damping_type.inner_damping_function = nullptr;
	description_inner =
	    "Damping " + description + " is disabled at inner boundary.";
	break;
    case damping_mean:
	damping_type.inner_damping_function =
	    &boundary_conditions::damping_single_inner_mean;
	description_inner =
	    "Damping " + description + " to mean value at inner boundary.";
	break;
    case damping_initial:
	damping_type.inner_damping_function =
	    &boundary_conditions::damping_single_inner;
	description_inner =
	    "Damping " + description + " to initial value at inner boundary.";
	break;
    case damping_zero:
	damping_type.inner_damping_function =
	    &boundary_conditions::damping_single_inner_zero;
	description_inner =
	    "Damping " + description + " to zero at inner boundary.";
	break;
    case damping_visc:
	damping_type.inner_damping_function =
	    &boundary_conditions::damping_vradial_inner_visc;
	description_inner = "Damping " + description +
			    " to viscous radial speed at inner boundary.";
	break;
    }

    switch (type_outer) {
    case damping_none:
	damping_type.outer_damping_function = nullptr;
	description_outer =
	    "Damping " + description + " is disabled at outer boundary.";
	break;
    case damping_mean:
	damping_type.outer_damping_function =
	    &boundary_conditions::damping_single_outer_mean;
	description_outer =
	    "Damping " + description + " to mean value at outer boundary.";
	break;
    case damping_initial:
	damping_type.outer_damping_function =
	    &boundary_conditions::damping_single_outer;
	description_outer =
	    "Damping " + description + " to initial value at outer boundary.";
	break;
    case damping_zero:
	damping_type.outer_damping_function =
	    &boundary_conditions::damping_single_outer_zero;
	description_outer =
	    "Damping " + description + " to zero at outer boundary.";
	break;
    case damping_visc:
	die(("Damping " + description +
	     " to viscous radial speed at outer boundary not implemented!\n")
		.c_str());
	break;
    }

    logging::print_master(LOG_INFO "%s\n", description_inner.c_str());
    logging::print_master(LOG_INFO "%s\n", description_outer.c_str());

    return damping_type;
}

/**
	Get a value as t_boundary_condition to a corresponding key  if
   available, else set to default

	\param key key
	\param defvalue default value
	\returns t_boundary_condition
*/
static t_damping_type value_as_boudary_damping_default(const char *key,
						const char *defvalue)
{
    const std::string ret = config::cfg.get<std::string>(key, defvalue);

    t_damping_type boundary_condition;
    switch (tolower(ret[0])) {
    case 'n':
	boundary_condition = parameters::t_damping_type::damping_none;
	break;
    case 'i':
	boundary_condition = parameters::t_damping_type::damping_initial;
	break;
    case 'y': // for legacy compatibility
	boundary_condition = parameters::t_damping_type::damping_initial;
	break;
    case 'm':
	boundary_condition = parameters::t_damping_type::damping_mean;
	break;
    case 'z':
	boundary_condition = parameters::t_damping_type::damping_zero;
	break;
    case 'v':
	boundary_condition = parameters::t_damping_type::damping_visc;
	break;
    default:
	boundary_condition = parameters::t_damping_type::damping_none;
    }
    return boundary_condition;

}

void read(const std::string &filename, t_data &data)
{
	if (filename.compare(filename.size()-4,4,".par") == 0) {
		die("This version of fargo uses the new yaml config files.\n%s looks like a par file.\nUse Tools/ini2yml.py to convert your config file!", filename.c_str());
	}
	config::cfg.load_file(filename);

    // units
	const std::string l0s = config::cfg.get<std::string>("l0", "");
	const std::string m0s = config::cfg.get<std::string>("m0", "");

	units::set_baseunits(l0s, m0s);

	units::precise_unit L0 = units::L0;
	units::precise_unit M0 = units::M0;
	units::precise_unit T0 = units::T0;
	units::precise_unit Temp0 = units::Temp0;

	constants::initialize_constants();

    // now we now everything to compute unit factors
    units::calculate_unit_factors();

    /* grid */
    NRadial = config::cfg.get<unsigned int>("NRAD", 64);
    NAzimuthal = config::cfg.get<unsigned int>("NSEC", 64);
    RMIN = config::cfg.get<double>("RMIN", L0);
    RMAX = config::cfg.get<double>("RMAX", L0);

    quantities_radius_limit =
	config::cfg.get<double>("QUANTITIESRADIUSLIMIT", 2.0 * RMAX, L0);

    if (quantities_radius_limit == 0.0) {
	quantities_radius_limit = 2.0 * RMAX;
    }

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
    data[t_data::MU].set_write(
	config::cfg.get_flag("WriteMeanMolecularWeight", false),
	do_write_1D);
    data[t_data::PRESSURE].set_write(
	config::cfg.get_flag("WritePressure", false), do_write_1D);
    data[t_data::TOOMRE].set_write(
	config::cfg.get_flag("WriteToomre", false), do_write_1D);
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
    data[t_data::ECCENTRICITY].set_write(
	config::cfg.get_flag("WriteEccentricity", false), do_write_1D);
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

    write_torques = config::cfg.get_flag("WriteTorques", false);

    write_disk_quantities =
	config::cfg.get_flag("WriteDiskQuantities", true);
    write_at_every_timestep =
	config::cfg.get_flag("WriteAtEveryTimestep", true);
    write_lightcurves =
	config::cfg.get_flag("WriteLightCurves", false);

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

    // boundary conditions
    switch (
	config::cfg.get_first_letter_lowercase("InnerBoundary", "Open")) {
    case 'o':
	boundary_inner = boundary_condition_open;
	break;
    case 'n':
	boundary_inner = boundary_condition_nonreflecting;
	break;
    case 'c':
	boundary_inner = boundary_condition_center_of_mass_initial;
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
    case 'j':
	boundary_inner = boundary_condition_jibin_spreading_ring;
	break;
    case 'b':
	boundary_inner = boundary_condition_boundary_layer;
	break;
    case 'k':
	boundary_inner = boundary_condition_keplerian;
	break;
    case 'p':
	boundary_inner = boundary_condition_precribed_time_variable;
	break;
    case 'z':
	boundary_inner = boundary_condition_zero_gradient;
	break;
    default:
	die("Invalid setting for InnerBoundary: %s",
	    config::cfg.get<std::string>("InnerBoundary", "Open").c_str());
    }

    switch (
		config::cfg.get_first_letter_lowercase("OuterBoundary", "Open")) {
    case 'o':
	boundary_outer = boundary_condition_open;
	break;
    case 'n':
	boundary_outer = boundary_condition_nonreflecting;
	break;
    case 'c':
	boundary_outer = boundary_condition_center_of_mass_initial;
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
    case 'j':
	boundary_outer = boundary_condition_jibin_spreading_ring;
	break;
    case 'b':
	boundary_outer = boundary_condition_boundary_layer;
	break;
    case 'k':
	boundary_outer = boundary_condition_keplerian;
	break;
    case 'p':
	boundary_outer = boundary_condition_precribed_time_variable;
	break;
    case 'z':
	boundary_outer = boundary_condition_zero_gradient;
	break;
    default:
	die("Invalid setting for OuterBoundary: %s",
	    config::cfg.get<std::string>("OuterBoundary", "Open").c_str());
    }

    // check if file for prescribed time variable boundary exists
    if (config::cfg.contains("PRESCRIBEDBOUNDARYFILEOUTER")) {
	if (config::cfg.get<std::string>("PRESCRIBEDBOUNDARYFILEOUTER").length() >
	    0) {
			PRESCRIBED_BOUNDARY_OUTER_FILE = config::cfg.get<std::string>("PRESCRIBEDBOUNDARYFILEOUTER") ;
	} else {
	    die("Error looking for data for the prescribed time variable boundary condition. Path could not be read!\n");
	}
    } else {
	PRESCRIBED_BOUNDARY_OUTER_FILE = "";
    }

    domegadr_zero = config::cfg.get_flag("DomegaDrZero", false);

    if (domegadr_zero)
	logging::print_master(
	    LOG_INFO "Using zero torque condition at outer boundary\n");

    viscous_outflow_speed =
	config::cfg.get<double>("ViscousOutflowSpeed", 1.0);

    damping = config::cfg.get_flag("Damping", false);

    damping_inner_limit =
	config::cfg.get<double>("DampingInnerLimit", 1.05);
    if (damping_inner_limit < 1) {
	die("DampingInnerLimit must not be <1\n");
    }
    damping_outer_limit =
	config::cfg.get<double>("DampingOuterLimit", 0.95);
    if (damping_outer_limit > 1) {
	die("DampingOuterLimit must not be >1\n");
    }
    damping_time_factor =
	config::cfg.get<double>("DampingTimeFactor", 1.0);

    t_damping_type tmp_damping_inner;
    t_damping_type tmp_damping_outer;

    if (config::cfg.contains("DampingVRadial"))
	die("DampingVRadial flag is decrepated used DampingVRadialInner and DampingVRadialOuter instead!");

    tmp_damping_inner =
	value_as_boudary_damping_default("DampingVRadialInner", "None");
    tmp_damping_outer =
	value_as_boudary_damping_default("DampingVRadialOuter", "None");

    if (tmp_damping_inner == parameters::t_damping_type::damping_visc) {
	damping_vector.push_back(
	    write_damping_type(tmp_damping_inner, tmp_damping_outer,
			       t_data::V_RADIAL, t_data::VISCOSITY, "VRadial"));
    } else {
	damping_vector.push_back(
	    write_damping_type(tmp_damping_inner, tmp_damping_outer,
			       t_data::V_RADIAL, t_data::V_RADIAL0, "VRadial"));
    }

    if (config::cfg.contains("DampingVAzimuthal"))
	die("DampingVRadial flag is decrepated used DampingVAzimuthalInner and DampingVAzimuthalOuter instead!");

    tmp_damping_inner = value_as_boudary_damping_default(
	"DampingVAzimuthalInner", "None");
    tmp_damping_outer = value_as_boudary_damping_default(
	"DampingVAzimuthalOuter", "None");

    damping_vector.push_back(write_damping_type(
	tmp_damping_inner, tmp_damping_outer, t_data::V_AZIMUTHAL,
	t_data::V_AZIMUTHAL0, "VAzimuthal"));

    if (config::cfg.contains("DampingSurfaceDensity"))
	die("DampingSurfaceDensity flag is decrepated used DampingSurfaceDensityInner and DampingSurfaceDensityOuter instead!");

    tmp_damping_inner = value_as_boudary_damping_default(
	"DampingSurfaceDensityInner", "None");
    tmp_damping_outer = value_as_boudary_damping_default(
	"DampingSurfaceDensityOuter", "None");

    damping_vector.push_back(
	write_damping_type(tmp_damping_inner, tmp_damping_outer, t_data::SIGMA,
			   t_data::SIGMA0, "SurfaceDensity"));

    if (config::cfg.contains("DampingEnergy"))
	die("DampingEnergy flag is decrepated used DampingEnergyInner and DampingEnergyOuter instead!");

    tmp_damping_inner =
	value_as_boudary_damping_default("DampingEnergyInner", "None");
    tmp_damping_outer =
	value_as_boudary_damping_default("DampingEnergyOuter", "None");

    damping_vector.push_back(
	write_damping_type(tmp_damping_inner, tmp_damping_outer, t_data::ENERGY,
			   t_data::ENERGY0, "Energy"));
    damping_energy_id = (int)damping_vector.size() - 1;

    calculate_disk = config::cfg.get_flag("DISK", "yes");
	
    corotation_reference_body =
	config::cfg.get<unsigned int>("CorotationReferenceBody", 1);

    // set number of bodies used to calculate barycenter in Interpret.cpp

    MU = config::cfg.get<double>("mu", 1.0);
    minimum_temperature =
	config::cfg.get<double>("MinimumTemperature", "3 K", Temp0);
    maximum_temperature =
	config::cfg.get<double>("MaximumTemperature", "1.0e300 K", Temp0);

    heating_viscous_enabled =
	config::cfg.get_flag("HeatingViscous", "no");
    heating_viscous_factor =
	config::cfg.get<double>("HeatingViscousFactor", 1.0);
    
	heating_star_enabled = config::cfg.get_flag("HeatingStar", "no");
	heating_star_factor =
	config::cfg.get<double>("HeatingStarFactor", 1.0);
    heating_star_simple =
	config::cfg.get_flag("HeatingStarSimple", "yes");

    radiative_diffusion_enabled =
	config::cfg.get_flag("RadiativeDiffusion", "no");
    radiative_diffusion_omega =
	config::cfg.get<double>("RadiativeDiffusionOmega", 1.5);
    radiative_diffusion_omega_auto_enabled =
	config::cfg.get_flag("RadiativeDiffusionAutoOmega", "no");
    radiative_diffusion_max_iterations = config::cfg.get<unsigned int>(
	"RadiativeDiffusionMaxIterations", 50000);

    zbuffer_size = config::cfg.get<unsigned int>("zbufferSize", 100);
    zbuffer_maxangle =
	config::cfg.get<double>("zbufferMaxAngle", 10.0 / 180.0 * M_PI);

    cooling_radiative_factor =
	config::cfg.get<double>("CoolingRadiativeFactor", 1.0);
    cooling_radiative_enabled =
	config::cfg.get_flag("CoolingRadiativeLocal", "no");
    cooling_beta_enabled =
	config::cfg.get_flag("CoolingBetaLocal", "no");
    cooling_beta = config::cfg.get<double>("CoolingBeta", 1.0);
    cooling_beta_ramp_up =
	config::cfg.get<double>("CoolingBetaRampUp", 0.0, T0);
	cooling_beta_aspect_ratio = false;
	cooling_beta_initial = false;
	switch (config::cfg.get_first_letter_lowercase("CoolingBetaReference", "Zero")) {
    case 'z': // Zero
	break;
    case 'i': // Initial
	cooling_beta_initial = true;
	break;
    case 'a': // AspectRatio
	cooling_beta_aspect_ratio = true;
	die("Not implemented yet: CoolingBetaReference: AspectRatio");
	break;
    default:
	die("Invalid setting for CoolingBetaReference: %s",
	    config::cfg.get<std::string>("CoolingBetaReference", "Zero").c_str());
    }

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
    case 's': // Boundary Layer: initialize w/ SS-73-Standard-Solution
	sigma_initialize_condition = initialize_condition_shakura_sunyaev;
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
    case 's': // Boundary Layer: initialize w/ SS-73-Standard-Solution
	energy_initialize_condition = initialize_condition_shakura_sunyaev;
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
    sigma0 = config::cfg.get<double>("SIGMA0", 173., M0/(L0*L0));
    sigma_adjust = config::cfg.get_flag("SetSigma0", "no");
    sigma_discmass = config::cfg.get<double>("discmass", 0.01, M0);
    density_factor = config::cfg.get<double>("DensityFactor", std::sqrt(2.0 * M_PI));

    tau_factor = config::cfg.get<double>("TauFactor", 0.5);
    kappa_factor = config::cfg.get<double>("KappaFactor", 1.0);

    EXPLICIT_VISCOSITY =
	config::cfg.get_flag("ExplicitViscosity", "yes");

    if (EXPLICIT_VISCOSITY) {
	logging::print_master(LOG_INFO "Using EXPLICIT VISCOSITY\n");
    } else {
	logging::print_master(LOG_INFO "Using SUPER TIMESTEPPINGG VISCOSITY\n");
    }

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
	config::cfg.get_flag("ArtificialViscosityDissipation", true);
    artificial_viscosity_factor =
	config::cfg.get<double>("ArtificialViscosityFactor", 1.41);
    // warning
    if (config::cfg.contains("CVNR")) {
	die("Parameter CVNR has been renamed to ArtificialViscosityFactor");
    }

    //
    thickness_smoothing =
	config::cfg.get<double>("ThicknessSmoothing", 0.6);
    thickness_smoothing_sg = config::cfg.get<double>(
	"ThicknessSmoothingSG", thickness_smoothing);
    integrate_planets = config::cfg.get_flag("IntegratePlanets", "yes");
    do_init_secondary_disk =
	config::cfg.get_flag("SecondaryDisk", "no");

    // mass overflow
    massoverflow = config::cfg.get_flag("massoverflow", "no");
    mof_planet = config::cfg.get<int>("mofplanet", 1);
    mof_temperature = config::cfg.get<double>("moftemperature", "1000.0 K", Temp0);
    mof_value = config::cfg.get<double>("mofvalue", 10E-9, M0/T0);
    mof_rampingtime = config::cfg.get<double>("moframpingtime", 30.0);

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

    // self gravity
    self_gravity = config::cfg.get_flag("SelfGravity", "no");

    if (self_gravity) {
	logging::print_master(LOG_INFO "Self gravity enabled.\n");
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

    // opacity
    switch (config::cfg.get_first_letter_lowercase("Opacity", "Lin")) {
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
	kappa_const = config::cfg.get<double>("KappaConst", 1.0);
	break;
    case 's': // simple, see Gennaro D'Angelo et al. 2003
	opacity = opacity_simple;
	kappa_const = config::cfg.get<double>("KappaConst", 1.0);
	break;
    default:
	die("Invalid setting for Opacity: %s",
	    config::cfg.get<std::string>("Opacity", "Lin").c_str());
    }

    if (heating_star_enabled) {
	if (star_radius*L0 < 0.1*units::solar_radius_in_au*units::au) {
	    die("Star radius is smaller than Jupiter with %.3e [R_sol]. This cannot be an active star\n",
		star_radius);
	}
    }

    // boundary layer parameters
    radial_viscosity_factor =
	config::cfg.get<double>("RadialViscosityFactor", 1.);
    vrad_fraction_of_kepler = config::cfg.get<double>("VRadIn", 1.6e-3);
    stellar_rotation_rate =
	config::cfg.get<double>("StellarRotation", 0.1);
    mass_accretion_rate =
	config::cfg.get<double>("MassAccretionRate", 1.e-9, M0/T0);
    accretion_radius_fraction =
	config::cfg.get<double>("MassAccretionRadius", 1.0);
    klahr_smoothing_radius = config::cfg.get<double>(
	"KlahrSmoothingRadius", accretion_radius_fraction);

    CFL = config::cfg.get<double>("CFL", 0.5);
    HEATING_COOLING_CFL_LIMIT =
	config::cfg.get<double>("HeatingCoolingCFLlimit", 5.0e-2);

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
	config::cfg.get<unsigned int>("ParticleSpeciesNumber", 5);
    particle_radius_increase_factor =
	config::cfg.get<double>("ParticleRadiusIncreaseFactor", 10.0);
    particle_eccentricity =
	config::cfg.get<double>("ParticleEccentricity", 0.0);
    particle_density = config::cfg.get<double>("ParticleDensity", "2.65 g/cm3", M0/(L0*L0*L0));
    particle_slope = config::cfg.get<double>(
	"ParticleSurfaceDensitySlope", SIGMASLOPE);
    particle_slope =
	-particle_slope; // particle distribution scales with  r^slope, so we
			 // introduces the minus here to make it r^-slope (same
			 // as for gas)
    particle_slope +=
	1.0; // particles are distributed over a whole simulation ring which
	     // introduces a factor 1/r for the particle surface density
    particle_minimum_radius =
	config::cfg.get<double>("ParticleMinimumRadius", RMIN);
    particle_maximum_radius =
	config::cfg.get<double>("ParticleMaximumRadius", RMAX);
    particle_minimum_escape_radius =
	config::cfg.get<double>("ParticleMinimumEscapeRadius", RMIN);
    particle_maximum_escape_radius =
	config::cfg.get<double>("ParticleMaximumEscapeRadius", RMAX);
    particle_gas_drag_enabled =
	config::cfg.get_flag("ParticleGasDragEnabled", true);
    particle_disk_gravity_enabled =
	config::cfg.get_flag("ParticleDiskGravityEnabled", false);
	particle_dust_diffusion =
	config::cfg.get_flag("ParticleDustDiffusion", false);



    // particle integrator
    switch (
	config::cfg.get_first_letter_lowercase("ParticleIntegrator", "s")) {
    case 'e': // Explicit
	integrator = integrator_explicit;
	break;
    case 'a': // Adaptive
	integrator = integrator_adaptive;
	break;
    case 's': // Semi-implicit
	integrator = integrator_semiimplicit;

	if (!particle_gas_drag_enabled) {
	    logging::print_master(
		LOG_ERROR
		"Do not use semi-implicit particle integrator without gas drag, use the explicit integrator instead.\n");
	}

	break;
    case 'm': // exponential midpoint
	integrator = integrator_exponential_midpoint;

	if (!particle_gas_drag_enabled) {
	    logging::print_master(
		LOG_ERROR
		"Do not use exponential midpoint particle integrator without gas drag, use the explicit integrator instead.\n");
	}

	break;
    case 'i': // Implicit
	integrator = integrator_implicit;

	if (!particle_gas_drag_enabled) {
	    logging::print_master(
		LOG_ERROR
		"Do not use implicit particle integrator without gas drag, use the explicit integrator instead.\n");
	}

	break;
    default:
	die("Invalid setting for Particle Integrator: %s	with key %s",
	    config::cfg.get<std::string>("ParticleIntegrator", "s").c_str(),
	    config::cfg.get_first_letter_lowercase("ParticleIntegrator", "s"));
    }

    if (CartesianParticles && ((integrator == integrator_implicit) ||
			       integrator == integrator_semiimplicit ||
			       integrator == integrator_exponential_midpoint)) {
	// implicit and semiimplicit integrator only implemented in polar
	// coordiantes, but forces can be calculated in cartesian coordinates
	CartesianParticles = false;
	ParticlesInCartesian = true;
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

    // boundary conditions
    switch (boundary_inner) {
    case boundary_condition_open:
	logging::print_master(
	    LOG_INFO "Using 'open boundary condition' at inner boundary.\n");
	break;
    case boundary_condition_reflecting:
	logging::print_master(
	    LOG_INFO
	    "Using 'reflecting boundary condition' at inner boundary.\n");
	break;
    case boundary_condition_zero_gradient:
	logging::print_master(
	    LOG_INFO
	    "Using 'zero gradient boundary condition' at inner boundary.\n");
	break;
    case boundary_condition_nonreflecting:
	logging::print_master(
	    LOG_INFO
	    "Using 'nonreflecting boundary condition' at inner boundary.\n");
	break;
    case boundary_condition_evanescent:
	logging::print_master(
	    LOG_INFO
	    "Using 'evanescent boundary condition' at inner boundary.\n");
	break;
    case boundary_condition_viscous_outflow:
	logging::print_master(
	    LOG_INFO
	    "Using 'viscous outflow boundary condition' at inner boundary.\n");
	break;
    case boundary_condition_jibin_spreading_ring:
	logging::print_master(
	    LOG_INFO
	    "Using 'boundary_condition_jibin_spreading_ring' at inner boundary.\n");
	break;
    case boundary_condition_boundary_layer:
	logging::print_master(
	    LOG_INFO
	    "Using 'boundary layer boundary conditions' at inner boundary.\n");
	break;
    case boundary_condition_keplerian:
	logging::print_master(
	    LOG_INFO
	    "Using 'keplarian boundary conditions' at inner boundary.\n");
	break;
    case boundary_condition_precribed_time_variable:
	die("Inner precribed time variable boundary condition is not implemented yet!\n");
	break;
    case boundary_condition_center_of_mass_initial:
	die("Inner boundary initial condition in center of mass is not implemented yet!\n");
	break;
    }

    switch (boundary_outer) {
    case boundary_condition_open:
	logging::print_master(
	    LOG_INFO "Using 'open boundary condition' at outer boundary.\n");
	break;
    case boundary_condition_reflecting:
	logging::print_master(
	    LOG_INFO
	    "Using 'reflecting boundary condition' at outer boundary.\n");
	break;
    case boundary_condition_zero_gradient:
	logging::print_master(
	    LOG_INFO
	    "Using 'zero gradient boundary condition' at outer boundary.\n");
	break;
    case boundary_condition_nonreflecting:
	logging::print_master(
	    LOG_INFO
	    "Using 'nonreflecting boundary condition' at outer boundary.\n");
	break;
    case boundary_condition_evanescent:
	if (domegadr_zero) {
	    die("domegadr_zero = true and evanescent outer boundary condition is not allowed!\n");
	}
	logging::print_master(
	    LOG_INFO
	    "Using 'evanescent boundary condition' at outer boundary.\n");
	break;
    case boundary_condition_viscous_outflow:
	logging::print_master(
	    LOG_INFO
	    "Using 'viscous outflow boundary condition' at outer boundary.\n");
	break;
    case boundary_condition_jibin_spreading_ring:
	logging::print_master(
	    LOG_INFO
	    "Using 'boundary_condition_jibin_spreading_ring' at outer boundary.\n");
	break;
    case boundary_condition_boundary_layer:
	if (domegadr_zero) {
	    die("domegadr_zero = true and boundary layer outer boundary condition is not allowed!\n");
	}
	logging::print_master(
	    LOG_INFO
	    "Using 'boundary layer boundary conditions' at outer boundary.\n");
	break;
    case boundary_condition_keplerian:
	logging::print_master(
	    LOG_INFO
	    "Using 'keplarian boundary conditions' at inner boundary.\n");
	break;
    case boundary_condition_precribed_time_variable:
	if (domegadr_zero) {
	    die("domegadr_zero = true and prescribed time variable outer boundary condition is not allowed!\n");
	}
	logging::print_master(
	    LOG_INFO
	    "Using 'time variable boundary conditions' at inner boundary.\n");
	break;
    case boundary_condition_center_of_mass_initial:
	logging::print_master(
	    LOG_INFO
	    "Using 'initial boundary in center of mass frame' at outer boundary.\n");
	break;
    }

    // Mass Transfer
    if (parameters::massoverflow) {
	logging::print_master(
	    LOG_INFO
	    "Mass Transfer from planet #%d of %g M_sun/orbit with Ts = %g K and ramping time t_ramp = %g P_orb.\n",
	    mof_planet, mof_value, mof_temperature, mof_rampingtime);
    }

    // Boundary layer
    if (boundary_inner == boundary_condition_boundary_layer) {
	logging::print_master(
	    LOG_INFO
	    "Boundary Layer: Radial velocity at inner boundary is %e * V_Kepler.\n",
	    vrad_fraction_of_kepler);
	logging::print_master(
	    LOG_INFO
	    "Boundary Layer: Stellar rotation rate is %f * Om_Kepler.\n",
	    stellar_rotation_rate);
    }
    if (boundary_outer == boundary_condition_boundary_layer) {
	logging::print_master(
	    LOG_INFO
	    "Boundary Layer: Mass Accretion Rate is %g Solar Masses per Year.\n",
	    mass_accretion_rate * units::mass.get_cgs_factor() /
		units::time.get_cgs_factor() * units::cgs_Year /
		units::cgs_Msol);
    }
    logging::print_master(
	LOG_INFO
	"Boundary Layer: Radial Viscosity is multiplied by a factor of %f.\n",
	radial_viscosity_factor);

    if (!damping) {
	logging::print_master(LOG_INFO "Damping at boundaries is disabled.\n");
	is_damping_initial = false;
    } else {
	is_damping_initial = false;
	for (unsigned int i = 0; i < damping_vector.size(); ++i) {
	    is_damping_initial =
		is_damping_initial ||
		(damping_vector[i].type_inner == damping_initial);
	    is_damping_initial =
		is_damping_initial ||
		(damping_vector[i].type_outer == damping_initial);
	}
    }

    logging::print_master(LOG_INFO "Surface density factor: %g\n",
			  density_factor);
    logging::print_master(LOG_INFO "Tau factor: %g\n", tau_factor);
    logging::print_master(LOG_INFO "Kappa factor: %g\n", kappa_factor);

    logging::print_master(LOG_INFO "Minimum temperature: %.5e K = %.5e\n",
			  minimum_temperature*units::temperature, minimum_temperature);
    logging::print_master(LOG_INFO "Maximum temperature: %.5e K = %.5e\n",
			  maximum_temperature*units::temperature, maximum_temperature);

    logging::print_master(
	LOG_INFO
	"Heating from star is %s. Using %s model with ramping time of %g and a total factor %g.\n",
	heating_star_enabled ? "enabled" : "disabled",
	heating_star_simple ? "simplified" : "advanced",
	heating_star_ramping_time, heating_star_factor);
    logging::print_master(
	LOG_INFO
	"Heating from viscous dissipation is %s. Using a total factor of %g.\n",
	heating_viscous_enabled ? "enabled" : "disabled",
	heating_viscous_factor);
	logging::print_master(LOG_INFO "Cooling (beta) is %s and reference temperature is %s. Using beta = %g.\n",
			  cooling_beta_enabled ? "enabled" : "disabled",
			  cooling_beta_initial ? "the initial value" : "zero",
			  cooling_beta);

    logging::print_master(
	LOG_INFO "Cooling (radiative) is %s. Using a total factor of %g.\n",
	cooling_radiative_enabled ? "enabled" : "disabled",
	cooling_radiative_factor);
    logging::print_master(
	LOG_INFO
	"Radiative diffusion is %s. Using %s omega = %lf with a maximum %u interations.\n",
	radiative_diffusion_enabled ? "enabled" : "disabled",
	radiative_diffusion_omega_auto_enabled ? "auto" : "fixed",
	radiative_diffusion_omega, radiative_diffusion_max_iterations);

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
    case opacity_zhu:
	logging::print_master(LOG_INFO
			      "Opacity uses tables from Zhu et al., 2012\n");
	break;
    case opacity_kramers:
	logging::print_master(
	    LOG_INFO
	    "Kramers opacity and constant electron scattering (Thomson) used.\n");
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
	char *buffer, *temp;
	unsigned int pos = 0;

	buffer = (char *)malloc((pos + 1) * sizeof(char));
	for (unsigned int i = 0; i < lightcurves_radii.size(); ++i) {
	    int size = asprintf(&temp, "%lf, ", lightcurves_radii[i]);
	    buffer = (char *)realloc(buffer, (pos + size + 1) * sizeof(char));
	    strncpy(&buffer[pos], temp, strlen(&buffer[pos]));
	    free(temp);
	    pos += size;
	}
	buffer[pos - 2] = '0';

	logging::print_master(LOG_INFO "Lightcurves radii are: %s\n", buffer);
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
	switch (integrator) {
	case integrator_explicit:
	    logging::print_master(LOG_INFO
				  "Particles use the explicit integrator\n");
	    break;
	case integrator_adaptive:
	    logging::print_master(
		LOG_INFO "Particles use the (explicit) adaptive integrator\n");
	    break;
	case integrator_semiimplicit: // Semi-implicit
	    logging::print_master(
		LOG_INFO "Particles use the semiimplicit integrator\n");
	    break;
	case integrator_exponential_midpoint:
	    logging::print_master(
		LOG_INFO "Particles use the exponential midpoint integrator\n");
	    break;
	case integrator_implicit: // Implicit
	    logging::print_master(LOG_INFO
				  "Particles use the implicit integrator\n");
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
