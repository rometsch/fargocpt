#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>

#include "Interpret.h"
#include "LowTasks.h"
#include "config.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "parameters.h"
#include "start_mode.h"
#include "units.h"

extern int damping_energy_id;
extern std::vector<parameters::t_DampingType> damping_vector;

// frame
int Corotating, GuidingCenter;

int FastTransport;
int OuterSourceMass, CICPlanet;

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include "options.h"
#include <sys/stat.h>

// used for calling python script needed for getting the polytropic constants
std::string exec(const char *cmd)
{
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
	throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
	result += buffer.data();
    }
    return result;
}

void get_polytropic_constants(double &K, double &gamma)
{
	// P_poly = K * Sigma**gamma
	// P_poly = K * Sigma0**gamma * r**(-p * gamma)
	// P_iso = Sigma * C_s**2
	// with Sigma = Sigma0 * r**(-p); C_s = h * vk * r**(F)
	// P_iso = Sigma0 * h**2 * G*M * r**(-1) * r**(2*F)
	// P_iso = Sigma0 * h**2 * G*M * r**(-1 - p + 2*F)
	// through comparisson of coefficients, we find the following relations:
	const double p = SIGMASLOPE;
	const double F = FLARINGINDEX;
	const double h = ASPECTRATIO_REF;
	gamma = (-1.0 - p + 2.0*F)/(-p);
	K = std::pow(h, 2) * std::pow(parameters::sigma0, 1.0-gamma);
}

std::string getFileName(const std::string &s)
{

    char sep = '/';

#ifdef _WIN32
    sep = '\\';
#endif

    size_t i = s.rfind(sep, s.length());
    if (i != std::string::npos) {
	return (s.substr(i + 1, s.length() - i));
    } else {
	return s;
    }

    return ("");
}

static void _mkdir(const char *dir, mode_t mode)
{
    // from
    // https://stackoverflow.com/questions/2336242/recursive-mkdir-system-call-on-unix
    char tmp[256];
    char *p = NULL;
    size_t len;

    snprintf(tmp, sizeof(tmp), "%s", dir);
    len = strlen(tmp);
    if (tmp[len - 1] == '/')
	tmp[len - 1] = 0;
    for (p = tmp + 1; *p; p++)
	if (*p == '/') {
	    *p = 0;
	    mkdir(tmp, mode);
	    *p = '/';
	}
    mkdir(tmp, S_IRWXU);
}

void ReadVariables(char *filename, t_data &data, int argc, char **argv)
{
    // read config from
    // config::read_config_from_file(filename);
    parameters::read(filename, data);

	constants::initialize_constants();

	// now we now everything to compute unit factors
	units::calculate_unit_factors();

	// TODO: This should definitely done in parameters.cpp, where values are
	// read, but parameters::read() is called before
	// units::calculate_unit_factors() so it is not possible. Moving the read()
	// call causes an error.
	parameters::apply_units();

	parameters::ShockTube = config::value_as_bool_default("ShockTube", 0);
	parameters::SpreadingRing =
	config::value_as_bool_default("SpreadingRing", NO);

    SIGMASLOPE = config::value_as_double_default("SIGMASLOPE", 0.0);
    IMPOSEDDISKDRIFT = config::value_as_double_default("IMPOSEDDISKDRIFT", 0.0);

    FLARINGINDEX = config::value_as_double_default("FLARINGINDEX", 0.0);

    std::string par_file_name = getFileName(filename);
    par_file_name = par_file_name.substr(0, par_file_name.size() - 4) + "/";
    if (asprintf(&OUTPUTDIR, "%s",
		 config::value_as_string_default("OUTPUTDIR",
						 par_file_name.c_str())) < 0) {
	logging::print_master(LOG_ERROR "Not enough memory!\n");
    }

    // Create output directory if it doesn't exist
    if (CPU_Master) {
	struct stat buffer;
	if (stat(OUTPUTDIR, &buffer)) {
	    _mkdir(OUTPUTDIR, 0700);
	}
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // check if planet config exists
    if ((config::key_exists("PLANETCONFIG")) &&
	(strlen(config::value_as_string("PLANETCONFIG")) > 0)) {
	if (asprintf(&PLANETCONFIG, "%s",
		     config::value_as_string("PLANETCONFIG")) < 0) {
	    logging::print_master(LOG_ERROR "Not enough memory!\n");
	}
    } else {
	PLANETCONFIG = NULL;
    }

    start_mode::configure_start_mode();

    if (CPU_Master) {
	// copy setup files into the output folder
	std::string output_folder = std::string(OUTPUTDIR);
	std::string par_file = getFileName(filename);
	if (output_folder.back() != '/') {
	    output_folder += "/";
	}
	std::string par_filename = output_folder + par_file;

	if (start_mode::mode == start_mode::mode_restart) {
	    char str[12];
	    sprintf(str, "%d", start_mode::restart_from);

	    par_filename += "_restart_";
	    par_filename += str;
	}

	std::remove(par_filename.c_str());
	std::ofstream new_par_file;
	new_par_file.open(par_filename.c_str(),
			  std::ofstream::out | std::ofstream::trunc);

	new_par_file << "###  Used launch options:";
	for (int i = 0; i < argc; i++) {
	    new_par_file << " " << argv[i];
	}
	new_par_file << "\n\n\n";
	new_par_file.close();

	std::filebuf old_par_file, append_new_par_file;
	old_par_file.open(filename, std::ios::in);
	append_new_par_file.open(par_filename.c_str(),
				 std::ios::app | std::ios::out);

	std::copy(std::istreambuf_iterator<char>(&old_par_file), {},
		  std::ostreambuf_iterator<char>(&append_new_par_file));
	append_new_par_file.close();
	old_par_file.close();

	if (PLANETCONFIG != NULL) {
	    std::string planet_file = std::string(PLANETCONFIG);
	    std::string planet_filename =
		output_folder + getFileName(planet_file);
	    if (start_mode::mode == start_mode::mode_restart) {
		char str[12];
		sprintf(str, "%d", start_mode::restart_from);
		planet_filename += "_restart_";
		planet_filename += str;
	    }

	    std::filebuf old_planet_file, append_new_planet_file;
	    old_planet_file.open(planet_file.c_str(), std::ios::in);
	    append_new_planet_file.open(planet_filename.c_str(),
					std::ios::trunc | std::ios::out);
	    std::copy(std::istreambuf_iterator<char>(&old_planet_file), {},
		      std::ostreambuf_iterator<char>(&append_new_planet_file));
	    append_new_planet_file.close();
	    old_planet_file.close();
	}
    }
    MPI_Barrier(MPI_COMM_WORLD);
    OuterSourceMass = config::value_as_bool_default("OUTERSOURCEMASS", 0);

    switch (tolower(*config::value_as_string_default("TRANSPORT", "Fast"))) {
    case 'f':
	FastTransport = 1;
	break;
    case 's':
	FastTransport = 0;
	break;
    default:
	die("Invalid setting for Transport");
    }

    // time settings
    NTOT = config::value_as_unsigned_int_default("NTOT", 1000);
    NINTERM = config::value_as_unsigned_int_default("NINTERM", 10);
    DT = config::value_as_double_default("DT", 1.0);

    if ((parameters::radial_grid_type == parameters::logarithmic_spacing) ||
	(parameters::radial_grid_type == parameters::exponential_spacing)) {
	double c = log(RMAX / RMIN);
	double optimal_N_azimuthal =
			M_PI / ((exp(c / NRadial) - 1.0) / (exp(c / NRadial) + 1.0));

	// check if optimal azimuthal cell number differs from actual azimuthal
	// cell number by more than 10%
	if (fabs(((double)NAzimuthal - optimal_N_azimuthal) /
		 (double)NAzimuthal) > 0.1) {
	    logging::print_master(
		LOG_WARNING
		"You have %u cells in azimuthal direction. This should be %u cells to have quadratic cells!\n",
		NAzimuthal, lround(optimal_N_azimuthal));
	}
    }

	dphi = 2.0 * M_PI / (double)NAzimuthal;
	invdphi = (double)NAzimuthal / (2.0 * M_PI);

    // disc
    ASPECTRATIO_REF = config::value_as_double_default("ASPECTRATIO", 0.05);

    if (!config::key_exists("OuterBoundary")) {
	logging::print_master(LOG_ERROR
			      "OuterBoundary doesn't exist. Old .par file?\n");
	die("died for convenience ;)");
    }

    // Frame settings
    Corotating = 0;
    GuidingCenter = 0;
    switch (tolower(*config::value_as_string_default("Frame", "Fixed"))) {
    case 'f': // Fixed
	break;
    case 'c': // Corotating
	Corotating = 1;
	break;
    case 'g': // Guiding-Center
	Corotating = 1;
	GuidingCenter = 1;
	break;
    default:
	die("Invalid setting for Frame");
    }
    OMEGAFRAME = config::value_as_double_default("OMEGAFRAME", 0);

    // Barycenter mode
    switch (tolower(
	*config::value_as_string_default("HydroFrameCenter", "primary"))) {
    case 'p': // primary
	parameters::n_bodies_for_hydroframe_center = 1;
	break;
    case 'b': // binary
	parameters::n_bodies_for_hydroframe_center = 2;
	break;
    case 't': // tertiary
	parameters::n_bodies_for_hydroframe_center = 3;
	break;
    case 'q': // quaternary
	parameters::n_bodies_for_hydroframe_center = 4;
	break;
    case 'a': // all
	parameters::n_bodies_for_hydroframe_center = 0;
	// will be set to number of bodies when loading planet system
	break;
    default:
	die("Invalid setting for HydroFrameCenter: %s",
	    config::value_as_string_default("HydroFrameCenter", "primary"));
    }

    parameters::exitOnDeprecatedSetting(
	"IndirectTerm",
	"Indirect terms are now handled automatically to avoid unphysical settings.",
	"Please remove the setting.");
    parameters::exitOnDeprecatedSetting(
	"IndirectTermPlanet",
	"Indirect terms are now handled automatically to avoid unphysical settings.",
	"Please remove the setting.");

    if (!parameters::Adiabatic) // if energy is not needed, delete the energy
				// damping boundary conditions
    {
	parameters::damping_vector.erase(parameters::damping_vector.begin() +
					 parameters::damping_energy_id);
    }

    // delete unneeded calls to damping functions
    // this must be performed after deleting the energy damping boundary,
    // otherwise damping_energy_id would be incorrect.
    auto delete_damping_condition =
	[&](const parameters::t_DampingType damper) {
	    return damper.inner_damping_function == nullptr &&
			   damper.outer_damping_function == nullptr
		       ? true
		       : false;
	};
    parameters::damping_vector.erase(
	std::remove_if(parameters::damping_vector.begin(),
		       parameters::damping_vector.end(),
		       delete_damping_condition),
	parameters::damping_vector.end());

    CICPlanet = config::value_as_bool_default("CICPLANET", 0);

    ALPHAVISCOSITY = config::value_as_double_default("ALPHAVISCOSITY", 0.0);
    VISCOSITY = config::value_as_double_default("VISCOSITY", 0.0);

	if (!EXPLICIT_VISCOSITY && ALPHAVISCOSITY == 0.0 &&
	(parameters::artificial_viscosity_factor == 0.0 ||
	 parameters::artificial_viscosity ==
		 parameters::artificial_viscosity_none) &&
	VISCOSITY == 0.0) {
	logging::print_master(
		LOG_ERROR
		"You cannot use super time-stepping without any viscosity!\n");
	PersonalExit(1);
	}

	STS_NU = config::value_as_double_default("STSNU", 0.01);

    if ((ALPHAVISCOSITY != 0.0) && (VISCOSITY != 0.0)) {
	logging::print_master(LOG_ERROR "You cannot use at the same time\n");
	logging::print_master(LOG_ERROR "VISCOSITY and ALPHAVISCOSITY.\n");
	logging::print_master(LOG_ERROR
			      "Edit the parameter file so as to remove\n");
	logging::print_master(LOG_ERROR
			      "one of these variables and run again.\n");
	PersonalExit(1);
    }

    if (ALPHAVISCOSITY != 0.0) {
	ViscosityAlpha = YES;
	logging::print_master(LOG_INFO "Viscosity is of alpha type\n");
    }

    if (parameters::thickness_smoothing <= 0.0) {
	logging::print_master(
	    LOG_ERROR
	    "A non-vanishing potential smoothing length is required.\n");
	}

    // Add a trailing slash to OUTPUTDIR if needed
    if (OUTPUTDIR[strlen(OUTPUTDIR) - 1] != '/') {
	unsigned int size = strlen(OUTPUTDIR);
	OUTPUTDIR = (char *)realloc(OUTPUTDIR, size + 2);
	OUTPUTDIR[size] = '/';
	OUTPUTDIR[size + 1] = 0;
    }

	const double T0 = config::value_as_double_default("TemperatureCGS0", 0.0);
	if (T0 != 0.0) // rescale ASPECTRATIO_REF according to cgs Temperature
	ASPECTRATIO_REF =
		sqrt(T0 * units::temperature.get_inverse_cgs_factor() *
		 constants::R / parameters::MU);

	const bool VISCOSITY_in_CGS =
	config::value_as_bool_default("VISCOSITYINCGS", false);
	if (VISCOSITY_in_CGS) {
	VISCOSITY =
		VISCOSITY * units::kinematic_viscosity.get_inverse_cgs_factor();
	}

	/// EoS can only be read after apply_units(), because fitting
	/// polytropic constants requires parameters::sigma0 to be in code units.
	// Energy equation / Adiabatic
	char Adiabatic_deprecated =
	tolower(*config::value_as_string_default("Adiabatic", "false"));

	if (Adiabatic_deprecated == 'n') {
	logging::print_master(
		LOG_INFO
		"Warning : Setting the isothermal equation of state with the flag 'Adiabatic   NO' is deprecated. Use 'EquationOfState   Isothermal' instead.\n");
	}
	if (Adiabatic_deprecated == 'y') {
	parameters::Adiabatic = true;
	logging::print_master(
		LOG_INFO
		"Warning : Setting the ideal equation of state with the flag 'Adiabatic    YES' is deprecated. Use 'EquationOfState   Adiabatic' instead.\n");

	ADIABATICINDEX =
		config::value_as_double_default("AdiabaticIndex", 7.0 / 5.0);
	if ((parameters::Adiabatic) && (ADIABATICINDEX == 1)) {
		logging::print_master(
		LOG_WARNING
		"You cannot have Adiabatic=true and AdiabatcIndex = 1. I decided to put Adiabatic=false, to  simulate a locally isothermal equation of state. Please check that it what you really wanted to do!\n");
		parameters::Adiabatic = false;
	}
	} else {
	char eos_string[512];
	strncpy(
		eos_string,
		config::value_as_string_default("EquationOfState", "Isothermal"),
		256); // same as MAXNAME from config.cpp
	for (char *t = eos_string; *t != '\0'; ++t) {
		*t = tolower(*t);
	}

	bool could_read_eos = false;
	if (strcmp(eos_string, "isothermal") == 0 ||
		strcmp(eos_string, "iso") == 0) {
		could_read_eos = true;
		parameters::Adiabatic = false;
		parameters::Polytropic = false;
		parameters::Locally_Isothermal = true;
		logging::print_master(LOG_INFO
				  "Using isothermal equation of state.\n");
	}
	if (strcmp(eos_string, "adiabatic") == 0 ||
		strcmp(eos_string, "ideal") == 0) {
		could_read_eos = true;

		// Energy equation / Adiabatic
		parameters::Adiabatic = true;

		char ADIABATICINDEX_string[512];
		strncpy(
		ADIABATICINDEX_string,
		config::value_as_string_default("AdiabaticIndex", "7.0/5.0"),
		256); // same as MAXNAME from config.cpp
		for (char *t = ADIABATICINDEX_string; *t != '\0'; ++t) {
		*t = tolower(*t);
		}

		if (strcmp(ADIABATICINDEX_string, "fit_isothermal") == 0 ||
		strcmp(ADIABATICINDEX_string, "fit isothermal") == 0) {
		logging::print_master(
			LOG_ERROR
			"Automatic AdiabatcIndex determination only available for polytropic equation of state\n");
		PersonalExit(1);
		} else {
		ADIABATICINDEX = config::value_as_double_default(
			"AdiabaticIndex", 7.0 / 5.0);
		}

		if ((parameters::Adiabatic) && (ADIABATICINDEX == 1)) {
		logging::print_master(
			LOG_WARNING
			"You cannot have Adiabatic=true and AdiabatcIndex = 1. I decided to put Adiabatic=false, to simulate a locally isothermal equation of state. Please check that it what you really wanted to do!\n");
		parameters::Adiabatic = false;
		}
		logging::print_master(LOG_INFO "Using ideal equation of state.\n");
	}

	if (strcmp(eos_string, "polytropic") == 0 ||
		strcmp(eos_string, "polytrop") == 0 ||
		strcmp(eos_string, "poly") == 0) {
		could_read_eos = true;

		// Equation of state / Polytropic
		parameters::Polytropic = true;
		double K = 0.0;
		double gamma = 0.0;

		char ADIABATICINDEX_string[512];
		strncpy(ADIABATICINDEX_string,
			config::value_as_string_default("AdiabaticIndex", "2.0"),
			256); // same as MAXNAME from config.cpp
		for (char *t = ADIABATICINDEX_string; *t != '\0'; ++t) {
		*t = tolower(*t);
		}

		if (strcmp(ADIABATICINDEX_string, "fit_isothermal") == 0 ||
		strcmp(ADIABATICINDEX_string, "fit isothermal") == 0) {
		get_polytropic_constants(K, gamma);
		ADIABATICINDEX = gamma;
		} else {
		ADIABATICINDEX =
			config::value_as_double_default("AdiabaticIndex", 2.0);
		}

		char POLYTROPIC_CONSTANT_string[512];
		strncpy(
		POLYTROPIC_CONSTANT_string,
		config::value_as_string_default("PolytropicConstant", "12.753"),
		256); // same as MAXNAME from config.cpp
		for (char *t = POLYTROPIC_CONSTANT_string; *t != '\0'; ++t) {
		*t = tolower(*t);
		}

		if (strcmp(POLYTROPIC_CONSTANT_string, "fit_isothermal") == 0 ||
		strcmp(POLYTROPIC_CONSTANT_string, "fit isothermal") == 0) {
		if (K == 0.0) // Call script only if needed
		{
			get_polytropic_constants(K, gamma);
		}
		POLYTROPIC_CONSTANT = K;
		} else {
		POLYTROPIC_CONSTANT = config::value_as_double_default(
			"PolytropicConstant", 12.753);
		}

		if ((parameters::Polytropic) && (ADIABATICINDEX == 1)) {
		logging::print_master(
			LOG_WARNING
			"You cannot have Polytropic=true and AdiabatcIndex = 1. I decided to put Polytropic=false, to simulate a locally isothermal equation of state. Please check that it what you really wanted to do!\n");
		parameters::Polytropic = false;
		}
		logging::print_master(LOG_INFO
				  "Using polytropic equation of state.\n");
	}

	if (!could_read_eos)
		die("Invalid setting for Energy Equation:   %s\n",
		config::value_as_string("EquationOfState"));
	}
}

void PrintUsage(char *execname)
{
    logging::print_master(
	LOG_ERROR
	"Usage : %s [-abcdeimnptvz] [-(0-9)] [-s number] [-f scaling] parameters file\n",
	execname);
    logging::print_master(
	LOG_ERROR
	"\n-a : Monitor mass and angular momentum at each timestep\n");
    logging::print_master(
	LOG_ERROR
	"-b : Adjust azimuthal velocity to impose strict centrifugal balance at t=0\n");
    logging::print_master(
	LOG_ERROR
	"-c : Sloppy CFL condition (checked at each DT, not at each timestep)\n");
    logging::print_master(
	LOG_ERROR
	"-d : Print some debugging information on 'stdout' at each timestep\n");
    logging::print_master(LOG_ERROR
			  "-e : Activate EU test problem torque file output\n");
    logging::print_master(
	LOG_ERROR
	"-f : Scale density array by 'scaling'. Useful to increase/decrease\n");
    logging::print_master(
	LOG_ERROR
	"     disk surface density after a restart, for instance.            \n");
    logging::print_master(
	LOG_ERROR "-i : tabulate Sigma profile as given by restart files\n");
    logging::print_master(
	LOG_ERROR
	"-n : Disable simulation. The program just reads parameters file\n");
    logging::print_master(LOG_ERROR
			  "-o : Overrides output directory of input file.\n");
    logging::print_master(
	LOG_ERROR "-p : Give profiling information at each time step\n");
    logging::print_master(
	LOG_ERROR
	"-s : Restart simulation, taking #'number' files as initial conditions\n");
    logging::print_master(
	LOG_ERROR
	"-v : Verbose mode. Tells everything about parameters file\n");
    logging::print_master(
	LOG_ERROR
	"-z : fake sequential built when evaluating sums on HD meshes\n");
    logging::print_master(
	LOG_ERROR "-(0-9) : only write initial (or restart) HD meshes,\n");
    logging::print_master(LOG_ERROR
			  "     proceed to the next nth output and exit\n");
    logging::print_master(
	LOG_ERROR
	"     This option must stand alone on one switch (-va -4 is legal, -v4a is not)\n");
    PersonalExit(1);
}

double TellNbOrbits(double time)
{
    return time / 2.0 / M_PI * sqrt(constants::G * 1.0 / 1.0 / 1.0 / 1.0);
}

double TellNbOutputs(double time) { return (time / DT / NINTERM); }

void TellEverything()
{
    double temp;

    if (!CPU_Master)
	return;

    logging::print_master(LOG_VERBOSE "Disc properties:\n");
    logging::print_master(LOG_VERBOSE "----------------\n");
    logging::print_master(LOG_VERBOSE "Inner Radius          : %g\n", RMIN);
    logging::print_master(LOG_VERBOSE "Outer Radius          : %g\n", RMAX);
    logging::print_master(LOG_VERBOSE "Aspect Ratio          : %g\n",
			  ASPECTRATIO_REF);
    logging::print_master(LOG_VERBOSE "VKep at inner edge    : %.3g\n",
			  sqrt(constants::G * 1.0 * (1. - 0.0) / RMIN));
    logging::print_master(LOG_VERBOSE "VKep at outer edge    : %.3g\n",
			  sqrt(constants::G * 1.0 / RMAX));
    /*
    logging::print_master(LOG_VERBOSE "boundary_inner        : %i\n",
    parameters::boundary_inner); logging::print_master(LOG_VERBOSE
    "boundary_outer        : %i\n", parameters::boundary_outer);
    */
    // temp=2.0*PI*parameters::sigma0/(2.0-SIGMASLOPE)*(pow(RMAX,2.0-SIGMASLOPE)
    // - pow(RMIN,2.0-SIGMASLOPE));	/* correct this and what follows... */
    // logging::print_master(LOG_VERBOSE "Initial Disk Mass             : %g\n",
    // temp); temp=2.0*PI*parameters::sigma0/(2.0-SIGMASLOPE)*(1.0 -
    // pow(RMIN,2.0-SIGMASLOPE)); logging::print_master(LOG_VERBOSE "Initial
    // Mass inner to r=1.0  : %g \n", temp);
    // temp=2.0*PI*parameters::sigma0/(2.0-SIGMASLOPE)*(pow(RMAX,2.0-SIGMASLOPE)
    // - 1.0); logging::print_master(LOG_VERBOSE "Initial Mass outer to r=1.0  :
    // %g \n", temp);
    logging::print_master(LOG_VERBOSE
			  "Travelling time for acoustic density waves :\n");
    temp = 2.0 / 3.0 / ASPECTRATIO_REF * (pow(RMAX, 1.5) - pow(RMIN, 1.5));
    logging::print_master(
	LOG_VERBOSE
	" * From Rmin to Rmax  : %.2g = %.2f orbits ~ %.1f outputs\n",
	temp, TellNbOrbits(temp), TellNbOutputs(temp));
    temp = 2.0 / 3.0 / ASPECTRATIO_REF * (pow(RMAX, 1.5) - pow(1.0, 1.5));
    logging::print_master(
	LOG_VERBOSE
	" * From r=1.0 to Rmax: %.2g = %.2f orbits ~ %.1f outputs\n",
	temp, TellNbOrbits(temp), TellNbOutputs(temp));
    temp = 2.0 / 3.0 / ASPECTRATIO_REF * (pow(1.0, 1.5) - pow(RMIN, 1.5));
    logging::print_master(
	LOG_VERBOSE
	" * From r=1.0 to Rmin: %.2g = %.2f orbits ~ %.1f outputs\n",
	temp, TellNbOrbits(temp), TellNbOutputs(temp));
    temp = 2.0 * M_PI * sqrt(RMIN * RMIN * RMIN / constants::G / 1.0);
    logging::print_master(LOG_VERBOSE
			  "Orbital time at Rmin  : %.3g ~ %.2f outputs\n",
			  temp, TellNbOutputs(temp));
    temp = 2.0 * M_PI * sqrt(RMAX * RMAX * RMAX / constants::G / 1.0);
    logging::print_master(LOG_VERBOSE
			  "Orbital time at Rmax  : %.3g ~ %.2f outputs\n",
			  temp, TellNbOutputs(temp));
    logging::print_master(LOG_VERBOSE "Sound speed :\n");
    logging::print_master(LOG_VERBOSE " * At unit radius     : %.3g\n",
			  ASPECTRATIO_REF * sqrt(constants::G * 1.0));
    logging::print_master(LOG_VERBOSE " * At outer edge      : %.3g\n",
			  ASPECTRATIO_REF * sqrt(constants::G * 1.0 / RMAX));
    logging::print_master(LOG_VERBOSE " * At inner edge      : %.3g\n",
			  ASPECTRATIO_REF * sqrt(constants::G * 1.0 / RMIN));
    logging::print_master(LOG_VERBOSE "Grid properties:\n");
    logging::print_master(LOG_VERBOSE "----------------\n");
    logging::print_master(LOG_VERBOSE "Number of (local) rings  : %d\n",
			  NRadial);
    logging::print_master(LOG_VERBOSE "Number of (global) rings : %d\n",
			  GlobalNRadial);
    logging::print_master(LOG_VERBOSE "Number of sectors        : %d\n",
			  NAzimuthal);
    logging::print_master(LOG_VERBOSE "Total (local) cells      : %d\n",
			  NRadial * NAzimuthal);
    logging::print_master(LOG_VERBOSE "Total (gobal) cells      : %d\n",
			  GlobalNRadial * NAzimuthal);
    logging::print_master(LOG_VERBOSE "Outputs properties:\n");
    logging::print_master(LOG_VERBOSE "-------------------\n");
    logging::print_master(
	LOG_VERBOSE "Time increment between outputs : %.3f = %.3f orbits\n",
	NINTERM * DT, TellNbOrbits(NINTERM * DT));
}
