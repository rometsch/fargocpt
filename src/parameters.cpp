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

#include <limits>
constexpr double DBL_EPSILON = std::numeric_limits<double>::epsilon();

namespace parameters
{

// energy euations
bool Adiabatic = false;
bool Polytropic = false;
bool Locally_Isothermal = false;

t_radial_grid radial_grid_type;
const char *radial_grid_names[] = {"logarithmic", "arithmetic", "exponential",
				   "custom"};
double exponential_cell_size_factor;

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
double cooling_beta_ramp_up;
double cooling_beta;

bool radiative_diffusion_enabled;
double radiative_diffusion_omega;
bool radiative_diffusion_omega_auto_enabled;
unsigned int radiative_diffusion_max_iterations;

bool centrifugal_balance;
t_initialize_condition sigma_initialize_condition;
std::string sigma_filename;
int random_seed;
bool sigma_randomize;
double sigma_random_factor;
double sigma_floor;
double sigma_feature_size;
bool sigma_adjust;
double sigma_discmass;
double sigma0;
t_initialize_condition energy_initialize_condition;
std::string energy_filename;

t_artificial_viscosity artificial_viscosity;
double artificial_viscosity_factor;
bool artificial_viscosity_dissipation;

bool calculate_disk;

// control centering of frame
bool default_star;
unsigned int n_bodies_for_hydroframe_center;
unsigned int corotation_reference_body;

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
double roche_smoothing_factor;
double thickness_smoothing_sg;

bool roche_smoothing_enabled;
bool thickness_smoothing_at_cell;

bool exclude_hill;

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
bool sloppy_cfl;

double L0;
double M0;

unsigned int number_of_particles;
bool integrate_particles;
double particle_radius;
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
t_particle_integrator integrator;

// for constant opacity
double kappa_const = 1.0;

t_DampingType write_damping_type(t_damping_type type_inner,
				 t_damping_type type_outer,
				 t_data::t_polargrid_type quantity,
				 t_data::t_polargrid_type quantity0,
				 std::string description)
{
    t_DampingType damping_type;
    damping_type.array_to_damp = quantity;
    damping_type.array_with_damping_values = quantity0;
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
    }

    damping_type.description_inner = description_inner;

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
    }

    damping_type.description_outer = description_outer;

    return damping_type;
}

/**
	Get a value as t_boundary_condition to a corresponding key  if
   available, else set to default

	\param key key
	\param defvalue default value
	\returns t_boundary_condition
*/
t_damping_type value_as_boudary_damping_default(const char *key,
						const char *defvalue)
{
    const std::string ret = config.get<std::string>(key, defvalue);
    const char *string_key = ret.c_str();

    t_damping_type boundary_condition;
    switch (tolower(*string_key)) {
    case 'n':
	boundary_condition = t_damping_type::damping_none;
	break;
    case 'i':
	boundary_condition = t_damping_type::damping_initial;
	break;
    case 'y': // for legacy compatibility
	boundary_condition = t_damping_type::damping_initial;
	break;
    case 'm':
	boundary_condition = t_damping_type::damping_mean;
	break;
    case 'z':
	boundary_condition = t_damping_type::damping_zero;
	break;
    default:
	boundary_condition = t_damping_type::damping_none;
    }
    return boundary_condition;
}

void read(char *filename, t_data &data)
{
    config.load_file(filename);

    // units
    L0 = config.get<double>("l0", 1.0);
    M0 = config.get<double>("m0", 1.0);

    // star parameters
    star_temperature = config.get<double>("StarTemperature", 5778);
    star_radius = config.get<double>("StarRadius", 0.009304813);

    parse_grid_config();

    parse_output_config(data);

    parse_boundary_config();

    parse_damping_config();

    parse_nbody_config();

    parse_disk_config();

    parse_initialization_config();

    parse_viscosity_config();

    parse_opacity_config();

    parse_massoverflow_config();

    parse_hydrosolver_config();

    parse_boundarylayer_config();

    parse_particle_config();
}

void apply_units()
{
    star_temperature /= units::temperature.get_cgs_factor();
    mass_accretion_rate = mass_accretion_rate *
			  (units::cgs_Msol / units::cgs_Year) * 1. /
			  units::mass_accretion_rate.get_cgs_factor();
    sigma0 /= units::surface_density.get_cgs_factor();
    particle_radius /= units::length.get_cgs_factor();
    particle_density /= units::density.get_cgs_factor();
}

void summarize_parameters()
{
    // artifical viscosity
    switch (artificial_viscosity) {
    case artificial_viscosity_none:
	logging::info_master("Using no artificial viscosity.\n");
	break;
    case artificial_viscosity_TW:
	logging::info_master(
	    "Using Tscharnuter-Winkler (1979) artificial viscosity with C = %lf.\n",
	    parameters::artificial_viscosity_factor);
	logging::info_master(
	    "Artificial viscosity is %s for dissipation.\n",
	    parameters::artificial_viscosity_dissipation ? "used" : "not used");
	break;
    case artificial_viscosity_SN:
	logging::info_master(
	    "Using Stone-Norman (1991, ZEUS-2D) artificial viscosity with C = %lf.\n",
	    parameters::artificial_viscosity_factor);
	logging::info_master(
	    "Artificial viscosity is %s for dissipation.\n",
	    parameters::artificial_viscosity_dissipation ? "used" : "not used");
	break;
    }

    // boundary conditions
    switch (boundary_inner) {
    case boundary_condition_open:
	logging::info_master(
	    "Using 'open boundary condition' at inner boundary.\n");
	break;
    case boundary_condition_reflecting:
	logging::info_master(
	    "Using 'reflecting boundary condition' at inner boundary.\n");
	break;
    case boundary_condition_nonreflecting:
	logging::info_master(
	    "Using 'nonreflecting boundary condition' at inner boundary.\n");
	break;
    case boundary_condition_evanescent:
	logging::info_master(
	    "Using 'evanescent boundary condition' at inner boundary.\n");
	break;
    case boundary_condition_viscous_outflow:
	logging::info_master(
	    "Using 'viscous outflow boundary condition' at inner boundary.\n");
	break;
    case boundary_condition_boundary_layer:
	logging::info_master(
	    "Using 'boundary layer boundary conditions' at inner boundary.\n");
	break;
    case boundary_condition_keplerian:
	logging::info_master(
	    "Using 'keplarian boundary conditions' at inner boundary.\n");
	break;
    }

    switch (boundary_outer) {
    case boundary_condition_open:
	logging::info_master(
	    "Using 'open boundary condition' at outer boundary.\n");
	break;
    case boundary_condition_reflecting:
	logging::info_master(
	    "Using 'reflecting boundary condition' at outer boundary.\n");
	break;
    case boundary_condition_nonreflecting:
	logging::info_master(
	    "Using 'nonreflecting boundary condition' at outer boundary.\n");
	break;
    case boundary_condition_evanescent:
	logging::info_master(
	    "Using 'evanescent boundary condition' at outer boundary.\n");
	break;
    case boundary_condition_viscous_outflow:
	logging::info_master(
	    "Using 'viscous outflow boundary condition' at outer boundary.\n");
	break;
    case boundary_condition_boundary_layer:
	logging::info_master(
	    "Using 'boundary layer boundary conditions' at outer boundary.\n");
	break;
    case boundary_condition_keplerian:
	logging::info_master(
	    "Using 'keplarian boundary conditions' at inner boundary.\n");
	break;
    }

    // Mass Transfer
    if (parameters::massoverflow) {
	logging::info_master(
	    "Mass Transfer of %g M_0/orbit will be spread on %i gridcells (sigma = %g).\n",
	    parameters::mof_value,
	    int(NAzimuthal * 3.0 * parameters::mof_sigma), mof_sigma);
    }

    // Boundary layer
    if (boundary_inner == boundary_condition_boundary_layer) {
	logging::info_master(
	    "Boundary Layer: Radial velocity at inner boundary is %e * V_Kepler.\n",
	    vrad_fraction_of_kepler);
	logging::info_master(
	    "Boundary Layer: Stellar rotation rate is %f * Om_Kepler.\n",
	    stellar_rotation_rate);
    }
    if (boundary_outer == boundary_condition_boundary_layer) {
	logging::info_master(
	    "Boundary Layer: Mass Accretion Rate is %g Solar Masses per Year.\n",
	    mass_accretion_rate * units::mass.get_cgs_factor() /
		units::time.get_cgs_factor() * units::cgs_Year /
		units::cgs_Msol);
    }
    logging::info_master(
	"Boundary Layer: Radial Viscosity is multiplied by a factor of %f.\n",
	radial_viscosity_factor);

    if (damping) {
	for (unsigned int i = 0; i < damping_vector.size(); ++i) {
	    logging::info_master("%s\n",
				 damping_vector[i].description_inner.c_str());
	    logging::info_master("%s\n",
				 damping_vector[i].description_outer.c_str());
	}
    } else {
	logging::info_master("Damping at boundaries is disabled.\n");
    }

    logging::info_master("Surface density factor: %g\n", density_factor);
    logging::info_master("Tau factor: %g\n", tau_factor);
    logging::info_master("Kappa factor: %g\n", kappa_factor);

    logging::info_master("Minimum temperature: %.5e\n", minimum_temperature);
    logging::info_master("Maximum temperature: %.5e\n", maximum_temperature);

    logging::info_master(
	"Heating from star is %s. Using %s model with ramping time of %g and a total factor %g.\n",
	heating_star_enabled ? "enabled" : "disabled",
	heating_star_simple ? "simplified" : "advanced",
	heating_star_ramping_time, heating_star_factor);
    logging::info_master(
	"Heating from viscous dissipation is %s. Using a total factor of %g.\n",
	heating_viscous_enabled ? "enabled" : "disabled",
	heating_viscous_factor);
    logging::info_master("Cooling (beta) is %s. Using beta = %g.\n",
			 cooling_beta_enabled ? "enabled" : "disabled",
			 cooling_beta);
    logging::info_master(
	"Cooling (radiative) is %s. Using a total factor of %g.\n",
	cooling_radiative_enabled ? "enabled" : "disabled",
	cooling_radiative_factor);
    logging::info_master(
	"Radiative diffusion is %s. Using %s omega = %lf with a maximum %u interations.\n",
	radiative_diffusion_enabled ? "enabled" : "disabled",
	radiative_diffusion_omega_auto_enabled ? "auto" : "fixed",
	radiative_diffusion_omega, radiative_diffusion_max_iterations);

    logging::info_master("CFL parameter: %g\n", CFL);

    switch (opacity) {
    case opacity_lin:
	logging::info_master(
	    "Opacity uses tables from Lin & Papaloizou, 1985\n");
	break;
    case opacity_bell:
	logging::info_master("Opacity uses tables from Bell & Lin, 1994\n");
	break;
    case opacity_zhu:
	logging::info_master("Opacity uses tables from Zhu et al., 2012\n");
	break;
    case opacity_kramers:
	logging::info_master(
	    "Kramers opacity and constant electron scattering (Thomson) used.\n");
	break;
    case opacity_const_op:
	logging::info_master("Using constant opacity kappa_R = %e.\n",
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
	    strcpy(&buffer[pos], temp);
	    free(temp);
	    pos += size;
	}
	buffer[pos - 2] = '0';

	logging::info_master("Lightcurves radii are: %s\n", buffer);
    }

    // particles
    logging::info_master("Particles are %s.\n",
			 integrate_particles ? "enabled" : "disabled");
    if (integrate_particles) {
	logging::info_master(
	    "Using %u particles with a radius of %g and a density of %g.\n",
	    number_of_particles, particle_radius, particle_density);
	logging::info_master(
	    "Distributing particles with a r^%.2g profile from %g to %g with a eccentricity from 0.0 to %g.\n",
	    particle_slope, particle_minimum_radius, particle_maximum_radius,
	    particle_eccentricity);
	logging::info_master(
	    "Particles are considered escaped from the system when they reach a distance of %g or %g.\n",
	    particle_minimum_escape_radius, particle_maximum_escape_radius);
	logging::info_master("Particles gas drag is %s.\n",
			     particle_gas_drag_enabled ? "enabled"
						       : "disabled");
	logging::info_master("Particles disk gravity is %s.\n",
			     particle_disk_gravity_enabled ? "enabled"
							   : "disabled");
	switch (integrator) {
	case integrator_explicit:
	    logging::info_master("Particles use the explicit integrator\n");
	    break;
	case integrator_adaptive:
	    logging::info_master(
		"Particles use the (explicit) adaptive integrator\n");
	    break;
	case integrator_semiimplicit: // Semi-implicit
	    logging::info_master("Particles use the semiimplicit integrator\n");
	    break;
	case integrator_implicit: // Implicit
	    logging::info_master("Particles use the implicit integrator\n");
	    break;
	default:
	    die("Invalid setting for Particle Integrator: %s",
		config.get<std::string>("ParticleIntegrator", "s").c_str());
	}
    }
}

void write_grid_data_to_file()
{
    /* Write a file containing the base units to the output folder. */

    FILE *fd = 0;
    char *fd_filename;

    if (CPU_Master) {
	if (asprintf(&fd_filename, "%s%s", OUTPUTDIR.c_str(),
		     "dimensions.dat") == -1) {
	    logging::error_master("Not enough memory for string buffer.\n");
	    PersonalExit(1);
	}
	fd = fopen(fd_filename, "w");
	if (fd == NULL) {
	    logging::error_master(
		"Can't write 'dimensions.dat' file. Aborting.\n");
	    PersonalExit(1);
	}

	free(fd_filename);

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
		RMAX, 0.0, 2 * PI, NRadial, NAzimuthal, 1, 1,
		radial_spacing_str);
	fclose(fd);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void parse_grid_config()
{
    /* grid */
    NRadial = config.get<unsigned int>("NRAD", 64);
    NAzimuthal = config.get<unsigned int>("NSEC", 64);
    RMIN = config.get<double>("RMIN", 1.0);
    RMAX = config.get<double>("RMAX", 1.0);

    exponential_cell_size_factor =
	config.get<double>("ExponentialCellSizeFactor", 1.41);
    switch (
	tolower(config.get<std::string>("RadialSpacing", "ARITHMETIC")[0])) {
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
}

void parse_output_config(t_data &data)
{
    data[t_data::DENSITY].set_write(config.get_flag("WriteDensity", true));
    data[t_data::V_RADIAL].set_write(config.get_flag("WriteVelocity", true));
    data[t_data::V_AZIMUTHAL].set_write(config.get_flag("WriteVelocity", true));
    data[t_data::ENERGY].set_write(config.get_flag("WriteEnergy", true));
    data[t_data::TEMPERATURE].set_write(
	config.get_flag("WriteTemperature", false));
    data[t_data::SOUNDSPEED].set_write(
	config.get_flag("WriteSoundSpeed", false));
    data[t_data::PRESSURE].set_write(config.get_flag("WritePressure", false));
    data[t_data::TOOMRE].set_write(config.get_flag("WriteToomre", false));
    data[t_data::TOOMRE_1D].set_write(
	config.get_flag("WriteRadialToomre", false));
    data[t_data::QPLUS].set_write(config.get_flag("WriteQPlus", false));
    data[t_data::QMINUS].set_write(config.get_flag("WriteQMinus", false));
    data[t_data::KAPPA].set_write(config.get_flag("WriteKappa", false));
    data[t_data::TAU_COOL].set_write(config.get_flag("WriteTauCool", false));
    data[t_data::ALPHA_GRAV].set_write(
	config.get_flag("WriteAlphaGrav", false));
    data[t_data::ALPHA_GRAV_1D].set_write(
	config.get_flag("WriteRadialAlphaGrav", false));
    data[t_data::ALPHA_GRAV_MEAN].set_write(
	config.get_flag("WriteAlphaGravMean", false));
    data[t_data::ALPHA_GRAV_MEAN_1D].set_write(
	config.get_flag("WriteRadialAlphaGravMean", false));
    data[t_data::ALPHA_REYNOLDS].set_write(
	config.get_flag("WriteAlphaReynolds", false));
    data[t_data::ALPHA_REYNOLDS_1D].set_write(
	config.get_flag("WriteRadialAlphaReynolds", false));
    data[t_data::ALPHA_REYNOLDS_MEAN].set_write(
	config.get_flag("WriteAlphaReynoldsMean", false));
    data[t_data::ALPHA_REYNOLDS_MEAN_1D].set_write(
	config.get_flag("WriteRadialAlphaReynoldsMean", false));
    data[t_data::VISCOSITY].set_write(config.get_flag("WriteViscosity", false));
    data[t_data::DIV_V].set_write(config.get_flag("WriteDivV", false));
    data[t_data::ECCENTRICITY].set_write(
	config.get_flag("WriteEccentricity", false));
    data[t_data::T_REYNOLDS].set_write(
	config.get_flag("WriteTReynolds", false));
    data[t_data::T_GRAVITATIONAL].set_write(
	config.get_flag("WriteTGravitational", false));
    data[t_data::P_DIVV].set_write(config.get_flag("WritepDV", false));
    data[t_data::TAU].set_write(config.get_flag("WriteTau", false));
    data[t_data::ASPECTRATIO].set_write(
	config.get_flag("WriteAspectRatio", false));
    data[t_data::VISIBILITY].set_write(
	config.get_flag("WriteVisibility", false));
    data[t_data::LUMINOSITY_1D].set_write(
	config.get_flag("WriteRadialLuminosity", false));
    data[t_data::DISSIPATION_1D].set_write(
	config.get_flag("WriteRadialDissipation", false));
    data[t_data::TAU_EFF].set_write(
	config.get_flag("WriteVerticalOpticalDepth", false));

    write_torques = config.get_flag("WriteTorques", false);

    write_disk_quantities = config.get_flag("WriteDiskQuantities", false);
    write_at_every_timestep = config.get_flag("WriteAtEveryTimestep", false);
    write_lightcurves = config.get_flag("WriteLightCurves", false);

    write_massflow = config.get_flag("WriteMassFlow", false);

    log_after_steps = config.get<double>("LogAfterSteps", 0);
    log_after_real_seconds = config.get<double>("LogAfterRealSeconds", 600.0);

    // parse light curve radii
    if (config.contains("WriteLightCurvesRadii")) {
	// get light curves radii string

	std::string lightcurve_config =
	    config.get<std::string>("WriteLightCurvesRadii");

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

void parse_boundary_config()
{
    // boundary conditions
    switch (tolower(config.get<std::string>("InnerBoundary", "Open")[0])) {
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
	die("Invalid setting for InnerBoundary: %s",
	    config.get<std::string>("InnerBoundary", "Open").c_str());
    }

    switch (tolower(config.get<std::string>("OuterBoundary", "Open")[0])) {
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
	die("Invalid setting for OuterBoundary: %s",
	    config.get<std::string>("OuterBoundary", "Open").c_str());
    }

    domegadr_zero = config.get_flag("DomegaDrZero", false);

    if (!config.contains("OuterBoundary")) {
	logging::error_master("OuterBoundary doesn't exist. Old .par file?\n");
	die("died for convenience ;)");
    }
}

void parse_damping_config()
{
    damping = config.get_flag("Damping", false);

    damping_inner_limit = config.get<double>("DampingInnerLimit", 1.05);
    if (damping_inner_limit < 1) {
	die("DampingInnerLimit must not be <1\n");
    }
    damping_outer_limit = config.get<double>("DampingOuterLimit", 0.95);
    if (damping_outer_limit > 1) {
	die("DampingOuterLimit must not be >1\n");
    }
    damping_time_factor = config.get<double>("DampingTimeFactor", 1.0);

    t_damping_type tmp_damping_inner;
    t_damping_type tmp_damping_outer;

    if (config.contains("DampingVRadial"))
	die("DampingVRadial flag is decrepated used DampingVRadialInner and DampingVRadialOuter instead!");

    tmp_damping_inner =
	value_as_boudary_damping_default("DampingVRadialInner", "None");
    tmp_damping_outer =
	value_as_boudary_damping_default("DampingVRadialOuter", "None");

    damping_vector.push_back(
	write_damping_type(tmp_damping_inner, tmp_damping_outer,
			   t_data::V_RADIAL, t_data::V_RADIAL0, "VRadial"));

    if (config.contains("DampingVAzimuthal"))
	die("DampingVRadial flag is decrepated used DampingVAzimuthalInner and DampingVAzimuthalOuter instead!");

    tmp_damping_inner =
	value_as_boudary_damping_default("DampingVAzimuthalInner", "None");
    tmp_damping_outer =
	value_as_boudary_damping_default("DampingVAzimuthalOuter", "None");

    damping_vector.push_back(write_damping_type(
	tmp_damping_inner, tmp_damping_outer, t_data::V_AZIMUTHAL,
	t_data::V_AZIMUTHAL0, "VAzimuthal"));

    if (config.contains("DampingSurfaceDensity"))
	die("DampingSurfaceDensity flag is decrepated used DampingSurfaceDensityInner and DampingSurfaceDensityOuter instead!");

    tmp_damping_inner =
	value_as_boudary_damping_default("DampingSurfaceDensityInner", "None");
    tmp_damping_outer =
	value_as_boudary_damping_default("DampingSurfaceDensityOuter", "None");

    damping_vector.push_back(write_damping_type(
	tmp_damping_inner, tmp_damping_outer, t_data::DENSITY, t_data::DENSITY0,
	"SurfaceDensity"));

    tmp_damping_inner =
	value_as_boudary_damping_default("DampingEnergy", "None");
    if (config.contains("DampingEnergy"))
	die("DampingEnergy flag is decrepated used DampingEnergyInner and DampingEnergyOuter instead!");

    tmp_damping_inner =
	value_as_boudary_damping_default("DampingEnergyInner", "None");
    tmp_damping_outer =
	value_as_boudary_damping_default("DampingEnergyOuter", "None");

    damping_vector.push_back(
	write_damping_type(tmp_damping_inner, tmp_damping_outer, t_data::ENERGY,
			   t_data::ENERGY0, "Energy"));
    damping_energy_id = (int)damping_vector.size() - 1;
}

void parse_nbody_config()
{
    default_star = config.get_flag("DefaultStar", true);
    corotation_reference_body =
	config.get<unsigned int>("CorotationReferenceBody", 1);
    thickness_smoothing = config.get<double>("ThicknessSmoothing", 0.0);
    integrate_planets = config.get_flag("IntegratePlanets", true);
	exclude_hill = config.get_flag("ExcludeHILL", false);
}

void parse_disk_config()
{
    calculate_disk = config.get_flag("DISK", true);
    disk_feedback = config.get_flag("DiskFeedback", true);

    MU = config.get<double>("mu", 1.0);
    minimum_temperature = config.get<double>("MinimumTemperature", 3);
    maximum_temperature = config.get<double>(
	"MaximumTemperature", std::numeric_limits<double>::max());
    if (maximum_temperature < 0) {
	maximum_temperature = std::numeric_limits<double>::max();
    }

    // TODO: remove temporary warning
    if (config.contains("HeatingViscous") == false) {
	die("please specify HeatingViscous in config file");
    }
    heating_viscous_enabled = config.get_flag("HeatingViscous", false);
    heating_viscous_factor = config.get<double>("HeatingViscousFactor", 1.0);
    heating_star_enabled = config.get_flag("HeatingStar", false);
    heating_star_factor = config.get<double>("HeatingStarFactor", 1.0);
    heating_star_ramping_time =
	config.get<double>("HeatingStarRampingTime", 0.0);
    heating_star_simple = config.get_flag("HeatingStarSimple", false);

    radiative_diffusion_enabled = config.get_flag("RadiativeDiffusion", false);
    radiative_diffusion_omega =
	config.get<double>("RadiativeDiffusionOmega", 1.5);
    radiative_diffusion_omega_auto_enabled =
	config.get_flag("RadiativeDiffusionAutoOmega", false);
    radiative_diffusion_max_iterations =
	config.get<unsigned int>("RadiativeDiffusionMaxIterations", 50000);

    zbuffer_size = config.get<unsigned int>("zbufferSize", 100);
    zbuffer_maxangle = config.get<double>("zbufferMaxAngle", 10.0 / 180.0 * PI);

    cooling_radiative_factor =
	config.get<double>("CoolingRadiativeFactor", 1.0);
    cooling_radiative_enabled = config.get_flag("CoolingRadiativeLocal", false);
    cooling_beta_enabled = config.get_flag("CoolingBetaLocal", false);
    cooling_beta = config.get<double>("CoolingBeta", 1.0);
    cooling_beta_ramp_up = config.get<double>("CoolingBetaRampUp", 0.0);
}

void parse_initialization_config()
{
    // initialisation
    initialize_pure_keplerian =
	config.get_flag("InitializePureKeplerian", false);
    initialize_vradial_zero = config.get_flag("InitializeVradialZero", false);
    centrifugal_balance = config.get_flag("CentrifugalBalance", false);

    switch (tolower(config.get<std::string>("SigmaCondition", "Profile")[0])) {
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
	die("Invalid setting for SigmaCondition: %s",
	    config.get<std::string>("SigmaCondition", "Profile").c_str());
    }

    if (config.contains("SigmaFilename")) {
	sigma_filename = config.get<std::string>("SigmaFilename", "");
    }

    switch (tolower(config.get<std::string>("EnergyCondition", "Profile")[0])) {
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
	die("Invalid setting for EnergyCondition: %s",
	    config.get<std::string>("EnergyCondition", "Profile").c_str());
    }

    if (config.contains("EnergyFilename")) {
	energy_filename = config.get<std::string>("EnergyFilename", "");
    }

    random_seed = config.get<int>("RandomSeed", 0);
    sigma_randomize = config.get_flag("RandomSigma", false);
    sigma_random_factor = config.get<double>("RandomFactor", 0.1);
    sigma_feature_size = config.get<double>("FeatureSize", (RMAX - RMIN) / 150);
    sigma_floor = config.get<double>("SigmaFloor", 1e-9);
    sigma0 = config.get<double>("SIGMA0", 173.);
    sigma_adjust = config.get_flag("SetSigma0", false);
    sigma_discmass = config.get<double>("discmass", 0.01);
    density_factor = config.get<double>("DensityFactor", 2.0);

    // profile damping
    profile_damping = config.get_flag("ProfileDamping", false);
    profile_damping_point = config.get<double>("ProfileDampingPoint", 0.0);
    profile_damping_width = config.get<double>("ProfileDampingWidth", 1.0);
}

void parse_viscosity_config()
{
    // artificial visocisty
    switch (tolower(config.get<std::string>("ArtificialViscosity", "SN")[0])) {
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
	    config.get<std::string>("ArtificialViscosity", "SN").c_str());
    }
    artificial_viscosity_dissipation =
	config.get_flag("ArtificialViscosityDissipation", true);
    artificial_viscosity_factor =
	config.get<double>("ArtificialViscosityFactor", 1.41);
    // warning
    if (config.contains("CVNR")) {
	die("Parameter CVNR has been renamed to ArtificialViscosityFactor");
    }
}

void parse_selfgravity_config()
{
    thickness_smoothing_sg =
	config.get<double>("ThicknessSmoothingSG", thickness_smoothing);
    self_gravity = config.get_flag("SelfGravity", false);
}

void parse_opacity_config()
{
    tau_factor = config.get<double>("TauFactor", 1.0);
    kappa_factor = config.get<double>("KappaFactor", 1.0);

    // opacity
    switch (tolower(config.get<std::string>("Opacity", "Lin")[0])) {
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
	kappa_const = config.get<double>("KappaConst", 1.0);
	break;
    default:
	die("Invalid setting for Opacity: %s",
	    config.get<std::string>("Opacity", "Lin").c_str());
    }
}

void parse_massoverflow_config()
{
    massoverflow = config.get_flag("massoverflow", false);
    mof_planet = config.get<unsigned int>("mofplanet", 0);
    mof_sigma = config.get<double>("mofsigma", 1.0);
    mof_value = config.get<double>("mofvalue", 10E-9);
}

void parse_hydrosolver_config()
{
    CFL = config.get<double>("CFL", 0.5);
    sloppy_cfl = config.get_flag("SloppyCFL", false);
}

void parse_boundarylayer_config()
{
    radial_viscosity_factor = config.get<double>("RadialViscosityFactor", 1.);
    vrad_fraction_of_kepler = config.get<double>("VRadIn", 1.6e-3);
    stellar_rotation_rate = config.get<double>("StellarRotation", 0.1);
    mass_accretion_rate = config.get<double>("MassAccretionRate", 1.e-9);
}

void parse_particle_config()
{
    CartesianParticles = config.get_flag("CartesianParticles", false);
    integrate_particles = config.get_flag("IntegrateParticles", false);
    number_of_particles = config.get<unsigned int>("NumberOfParticles", 0);
    particle_radius = config.get<double>("ParticleRadius", 100.0);
    particle_eccentricity = config.get<double>("ParticleEccentricity", 0.0);
    particle_density = config.get<double>("ParticleDensity", 2.65);
    particle_slope =
	config.get<double>("ParticleSurfaceDensitySlope", SIGMASLOPE);
    particle_slope =
	-particle_slope; // particle distribution scales with  r^slope, so we
			 // introduces the minus here to make it r^-slope (same
			 // as for gas)
    particle_slope +=
	1.0; // particles are distributed over a whole simulation ring which
	     // introduces a factor 1/r for the particle surface density
    particle_minimum_radius = config.get<double>("ParticleMinimumRadius", RMIN);
    particle_maximum_radius = config.get<double>("ParticleMaximumRadius", RMAX);
    particle_minimum_escape_radius =
	config.get<double>("ParticleMinimumEscapeRadius", RMIN);
    particle_maximum_escape_radius =
	config.get<double>("ParticleMaximumEscapeRadius", RMAX);
    particle_gas_drag_enabled = config.get_flag("ParticleGasDragEnabled", true);
    particle_disk_gravity_enabled =
	config.get_flag("ParticleDiskGravityEnabled", false);
    // particle integrator
    switch (tolower(config.get<std::string>("ParticleIntegrator", "s")[0])) {
    case 'e': // Explicit
	integrator = integrator_explicit;
	break;
    case 'a': // Adaptive
	integrator = integrator_adaptive;
	break;
    case 's': // Semi-implicit
	integrator = integrator_semiimplicit;

	if (!particle_gas_drag_enabled) {
	    logging::error_master(
		"Do not use semi-implicit particle integrator without gas drag, use the explicit integrator instead.\n");
	}

	break;
    case 'i': // Implicit
	integrator = integrator_implicit;

	if (!particle_gas_drag_enabled) {
	    logging::error_master(
		"Do not use implicit particle integrator without gas drag, use the explicit integrator instead.\n");
	}

	break;
    default:
	die("Invalid setting for Particle Integrator: %s	with key %s",
	    config.get<std::string>("ParticleIntegrator", "s").c_str(),
	    tolower(config.get<std::string>("ParticleIntegrator", "s")[0]));
    }

    if (CartesianParticles && ((integrator == integrator_implicit) ||
			       integrator == integrator_semiimplicit)) {
	// implicit and semiimplicit integrator only implemented in polar
	// coordiantes, but forces can be calculated in cartesian coordinates
	CartesianParticles = false;
	ParticlesInCartesian = true;
    }

    if (particle_disk_gravity_enabled && (!self_gravity)) {
	logging::error_master(
	    "Cannot enable particle_disk_gravity_enabled while self_gravity is off!\n");
	PersonalExit(1);
    }

    // second and second last radius, so that partices stay out of the ghost
    // cells
    if (particle_minimum_escape_radius < RMIN) {
	logging::warning_master(
	    "particle_minimum_escape_radius can't be smaller than the inner radius of the domain. Setting particle_minimum_escape_radius to inner radius of the domain.\n");
	particle_minimum_escape_radius = RMIN;
    }

    if (particle_maximum_escape_radius > RMAX) {
	logging::warning_master(
	    "particle_maximum_escape_radius can't be larger than the outer radius of the domain. Setting particle_maximum_escape_radius to outer radius of the domain.\n");
	particle_maximum_escape_radius = RMAX;
    }

    particle_maximum_escape_radius_sq =
	pow2(particle_maximum_escape_radius) - DBL_EPSILON; // DBL for safety
    particle_minimum_escape_radius_sq =
	pow2(particle_minimum_escape_radius) + DBL_EPSILON;
}

} // namespace parameters
