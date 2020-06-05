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
int OuterSourceMass;

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

void get_polytropic_constants(char *filename, double &K, double &gamma)
{
    char command[1024];
    strcpy(command, "python Tools/get_polytropic_constants.py ");
    strcat(command, filename);
    strcat(command, " 0");
    // Python script fits the parameters ADIABATICINDEX and POLYTROPIC_CONSTANT,
    //      so that the pressure is equal to the pressure generated by the
    //      isothermal equation of state.
    std::string polytropic_string = exec(command);
    std::stringstream ss(polytropic_string);

    ss >> K;
    ss >> gamma;
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

static char lowercase_first_letter(const std::string &s)
{
    return tolower(s[0]);
}

void ReadVariables(char *filename, t_data &data, int argc, char **argv)
{
    // read config from
    // config::read_config_from_file(filename);
    parameters::read(filename, data);

    MASSTAPER = config.get("MASSTAPER", 0.0000001);
    ROCHESMOOTHING = config.get("ROCHESMOOTHING", 0.0);
    SIGMASLOPE = config.get("SIGMASLOPE", 0.0);
    IMPOSEDDISKDRIFT = config.get("IMPOSEDDISKDRIFT", 0.0);

    FLARINGINDEX = config.get("FLARINGINDEX", 0.0);

    std::string par_file_name = getFileName(filename);
    par_file_name = par_file_name.substr(0, par_file_name.size() - 4) + "/";
    OUTPUTDIR = config.get<std::string>("OUTPUTDIR", par_file_name);

    // Create output directory if it doesn't exist
    if (CPU_Master) {
	struct stat buffer;
	if (stat(OUTPUTDIR.c_str(), &buffer)) {
	    _mkdir(OUTPUTDIR.c_str(), 0700);
	}
    }

    MPI_Barrier(MPI_COMM_WORLD);

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
    }

    MPI_Barrier(MPI_COMM_WORLD);
    OuterSourceMass = config.get_flag("OUTERSOURCEMASS", false);

    switch (tolower(config.get<std::string>("TRANSPORT", "Fast")[0])) {
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
    NTOT = config.get("NTOT", 1000);
    NINTERM = config.get("NINTERM", 10);
    DT = config.get("DT", 1.0);

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

    // disc
    ASPECTRATIO_REF = config.get("ASPECTRATIO", 0.05);

    if (!config.contains("OuterBoundary")) {
	logging::print_master(LOG_ERROR
			      "OuterBoundary doesn't exist. Old .par file?\n");
	die("died for convenience ;)");
    }

    // Frame settings
    Corotating = 0;
    GuidingCenter = 0;
    switch (tolower(config.get<std::string>("Frame", "Fixed")[0])) {
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
    OMEGAFRAME = config.get("OMEGAFRAME", double(0));

    // Barycenter mode
    switch (lowercase_first_letter(config.get<std::string>("HydroFrameCenter", "primary"))) {
    case 'p': // primarys
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
	    config.get<std::string>("HydroFrameCenter", "primary").c_str());
    }

    if (config.contains("IndirectTerm")) {
	die("Indirect terms are now handled automatically to avoid unphysical settings.",
	    "Please remove the setting.");
    }

    if (config.contains("IndirectTermPlanet")) {
	die("Indirect terms are now handled automatically to avoid unphysical settings.",
	    "Please remove the setting.");
    }

    // Energy equation / Adiabatic
    char Adiabatic_deprecated =
	lowercase_first_letter(config.get<std::string>("Adiabatic", "false"));

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

	ADIABATICINDEX = config.get("AdiabaticIndex", 7.0 / 5.0);
	if ((parameters::Adiabatic) && (ADIABATICINDEX == 1)) {
	    logging::print_master(
		LOG_WARNING
		"You cannot have Adiabatic=true and AdiabatcIndex = 1. I decided to put Adiabatic=false, to  simulate a locally isothermal equation of state. Please check that it what you really wanted to do!\n");
	    parameters::Adiabatic = false;
	}
    } else {

	const std::string eos_string_tmp =
	    lowercase(config.get<std::string>("EquationOfState", "Isothermal"));
	const char *eos_string = eos_string_tmp.c_str();

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

	    const std::string ADIABATICINDEX_string_tmp =
		lowercase(config.get<std::string>("AdiabaticIndex", "7/5"));
	    const char *ADIABATICINDEX_string =
		ADIABATICINDEX_string_tmp.c_str();

	    if (strcmp(ADIABATICINDEX_string, "fit_isothermal") == 0 ||
		strcmp(ADIABATICINDEX_string, "fit isothermal") == 0) {
		logging::print_master(
		    LOG_ERROR
		    "Automatic AdiabatcIndex determination only available for polytropic equation of state\n");
		PersonalExit(1);
	    } else {
		ADIABATICINDEX = config.get("AdiabaticIndex", 7.0 / 5.0);
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

	    const std::string ADIABATICINDEX_string_tmp =
		lowercase(config.get<std::string>("AdiabaticIndex", "2.0"));
	    const char *ADIABATICINDEX_string =
		ADIABATICINDEX_string_tmp.c_str();

	    if (strcmp(ADIABATICINDEX_string, "fit_isothermal") == 0 ||
		strcmp(ADIABATICINDEX_string, "fit isothermal") == 0) {
		get_polytropic_constants(filename, K, gamma);
		ADIABATICINDEX = gamma;
	    } else {
		ADIABATICINDEX = config.get("AdiabaticIndex", 2.0);
	    }

	    const std::string POLYTROPIC_CONSTANT_string_tmp =
		lowercase(config.get<std::string>("AdiabaticIndex", "7/5"));
	    const char *POLYTROPIC_CONSTANT_string =
		POLYTROPIC_CONSTANT_string_tmp.c_str();

	    if (strcmp(POLYTROPIC_CONSTANT_string, "fit_isothermal") == 0 ||
		strcmp(POLYTROPIC_CONSTANT_string, "fit isothermal") == 0) {
		if (K == 0.0) // Call script only if needed
		{
		    get_polytropic_constants(filename, K, gamma);
		}
		POLYTROPIC_CONSTANT = K;
	    } else {
		POLYTROPIC_CONSTANT = config.get("PolytropicConstant", 12.753);
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

	if (!could_read_eos) {
	    const std::string eos_str = config.get("EquationOfState");
	    die("Invalid setting for Energy Equation:   %s\n", eos_str.c_str());
	}
    }

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

    ExcludeHill = config.get_flag("EXCLUDEHILL", false);
    ALPHAVISCOSITY = config.get("ALPHAVISCOSITY", 0.0);
    VISCOSITY = config.get("VISCOSITY", 0.0);

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

    if ((parameters::thickness_smoothing != 0.0) && (ROCHESMOOTHING != 0.0)) {
	logging::print_master(LOG_ERROR "You cannot use at the same time\n");
	logging::print_master(LOG_ERROR
			      "`ThicknessSmoothing' and `RocheSmoothing'.\n");
	logging::print_master(LOG_ERROR
			      "Edit the parameter file so as to remove\n");
	logging::print_master(LOG_ERROR
			      "one of these variables and run again.\n");
	PersonalExit(1);
    }

    if ((parameters::thickness_smoothing <= 0.0) && (ROCHESMOOTHING <= 0.0)) {
	logging::print_master(
	    LOG_ERROR
	    "A non-vanishing potential smoothing length is required.\n");
	logging::print_master(
	    LOG_ERROR "Please use either of the following variables:\n");
	logging::print_master(LOG_ERROR
			      "`ThicknessSmoothing' *or* `RocheSmoothing'.\n");
	logging::print_master(LOG_ERROR "before launching the run again.\n");
	PersonalExit(1);
    }

    if (ROCHESMOOTHING != 0.0) {
	RocheSmoothing = YES;
	logging::print_master(
	    LOG_INFO
	    "Planet potential smoothing scales with their Hill sphere.\n");
    } else if (config.get_flag("ThicknessSmoothingAtPlanet", false)) {
	ThicknessSmoothingAtCell = NO;
	ThicknessSmoothingAtPlanet = YES;
	logging::print_master(
	    LOG_INFO
	    "Planet potential smoothing uses disk scale height at planet location (bad choice!).\n");
    } else {
	ThicknessSmoothingAtCell = YES;
	ThicknessSmoothingAtPlanet = NO;
	logging::print_master(
	    LOG_INFO
	    "Planet potential smoothing uses disk scale height at gas cell location.\n");
    }

    // Add a trailing slash to OUTPUTDIR if needed
    if (OUTPUTDIR[OUTPUTDIR.length() - 1] != '/') {
	OUTPUTDIR.append("/");
    }

    constants::initialize_constants();

    // now we now everything to compute unit factors
    units::calculate_unit_factors();

    // TODO: This should definitely done in parameters.cpp, where values are
    // read, but parameters::read() is called before
    // units::calculate_unit_factors() so it is not possible. Moving the read()
    // call causes an error.
    parameters::apply_units();
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
    return time / 2.0 / PI * sqrt(constants::G * 1.0 / 1.0 / 1.0 / 1.0);
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
    temp = 2.0 * PI * sqrt(RMIN * RMIN * RMIN / constants::G / 1.0);
    logging::print_master(LOG_VERBOSE
			  "Orbital time at Rmin  : %.3g ~ %.2f outputs\n",
			  temp, TellNbOutputs(temp));
    temp = 2.0 * PI * sqrt(RMAX * RMAX * RMAX / constants::G / 1.0);
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
