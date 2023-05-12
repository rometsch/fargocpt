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
#include "config.h"
#include "output.h"
#include "simulation.h"

extern int damping_energy_id;
extern std::vector<parameters::t_DampingType> damping_vector;

// frame
int GuidingCenter;

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
#include <filesystem>

static void get_polytropic_constants(double &K, double &gamma)
{
    // P_poly = K * Sigma**gamma
    // P_poly = K * Sigma0**gamma * r**(-p * gamma)
    // P_iso = Sigma * C_s**2
    // with Sigma = Sigma0 * r**(-p); C_s = h * vk * r**(F)
    // P_iso = Sigma0 * h**2 * G*M * r**(-1) * r**(2*F)
    // P_iso = Sigma0 * h**2 * G*M * r**(-1 - p + 2*F)
    // through comparisson of coefficients, we find the following relations:
    const double p = parameters::SIGMASLOPE;
    const double F = parameters::FLARINGINDEX;
    const double h = parameters::ASPECTRATIO_REF;
    gamma = (-1.0 - p + 2.0 * F) / (-p);
    K = std::pow(h, 2) * std::pow(parameters::sigma0, 1.0 - gamma);
}

static std::string getFileName(const std::string &s)
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

void ReadVariables(const std::string &filename, t_data &data, int argc, char **argv)
{
    // read config from
	if (!std::filesystem::exists(filename)) {
		die("Can not find config file %s!\n", filename.c_str());
	}
    parameters::read(filename, data);
	auto &cfg = config::cfg;

    parameters::ShockTube = cfg.get<int>("ShockTube", 0);
    parameters::SpreadingRing =
	cfg.get_flag("SpreadingRing", false);

    sim::last_dt = cfg.get<double>("FirstDT", 1e-9, units::T0);

    parameters::SIGMASLOPE = cfg.get<double>("SIGMASLOPE", 0.0);
    parameters::IMPOSEDDISKDRIFT = cfg.get<double>("IMPOSEDDISKDRIFT", 0.0);

    parameters::FLARINGINDEX = cfg.get<double>("FLARINGINDEX", 0.0);

    std::string setup_name = getFileName(filename);
    setup_name = setup_name.substr(0, setup_name.size() - 4) + "/";

	output::outdir = cfg.get<std::string>("OUTPUTDIR", setup_name);
	if (output::outdir[output::outdir.length()-1] != '/') {
		output::outdir += "/";
	}
    
    ensure_directory_exists(output::outdir);
    MPI_Barrier(MPI_COMM_WORLD);

    start_mode::configure_start_mode();

    ensure_directory_exists(output::outdir + "snapshots/");
    ensure_directory_exists(output::outdir + "parameters/");
	ensure_directory_exists(output::outdir + "monitor/");
    MPI_Barrier(MPI_COMM_WORLD);

    if (CPU_Master) {

	// copy setup files into the output folder
	std::string output_folder = output::outdir + "parameters";
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

	// delete file with 'filename', do not confuse with overloaded 'remove' function that removes a part of the string
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
    OuterSourceMass = cfg.get_flag("OUTERSOURCEMASS", "no");

    switch (cfg.get_first_letter_lowercase("TRANSPORT", "Fast")) {
    case 'f':
	parameters::fast_transport = true;
	logging::print_master(LOG_INFO
	"Using FARGO algorithm for azimuthal advection.\n");
	break;
    case 's':
	parameters::fast_transport = false;
	logging::print_master(LOG_INFO
	"Using standard advection WITHOUT FARGO algorithm.\n");
	break;
    default:
	die("Invalid setting for Transport");
    }
    
	switch (cfg.get_first_letter_lowercase("Integrator", "Euler")) {
    case 'e':
	parameters::hydro_integrator = EULER_INTEGRATOR;
	logging::print_master(LOG_INFO
	"Using standard forward euler scheme for source terms.\n");
	break;
	case 'k':
	parameters::hydro_integrator = LEAPFROG_KICK_DRIFT_KICK;
	logging::print_master(LOG_INFO
	"Using leapfrog scheme for source terms.\n");
	break;
	case 'l':
	parameters::hydro_integrator = LEAPFROG_DRIFT_KICK_DRIFT;
	logging::print_master(LOG_INFO
	"Using leapfrog (gas: kick drift kick) (nbody: drift kick drift) scheme for source terms.\n");
	break;
    default:
	die("Invalid setting for Integrator (valid: Euler, Leapfrog");
    }
    


    // disc
    parameters::ASPECTRATIO_REF = cfg.get<double>("ASPECTRATIO", 0.05);
    parameters::ASPECTRATIO_MODE = cfg.get<int>("AspectRatioMode", 0);

    const double T0 = cfg.get<double>("Temperature0", "-1", units::Temp0);
	if (T0 > 0.0){ // rescale parameters::ASPECTRATIO_REF according to Temperature
	parameters::ASPECTRATIO_REF = sqrt(T0 * constants::R / parameters::MU);
	}

    // time settings
    parameters::NTOT = cfg.get<unsigned int>("NTOT", 1000);
    parameters::NINTERM = cfg.get<unsigned int>("NINTERM", 10);
    parameters::DT = cfg.get<double>("DT", 1.0);



	parameters::cps = config::cfg.get<double>("cps", -1.0);
	
	const double H = parameters::ASPECTRATIO_REF; // H(r=1)
	if (parameters::cps > 0) {
		if (config::cfg.contains("Nrad") || config::cfg.contains("Nsec")) {
			logging::print_master(LOG_INFO "Cps is set, overwriting Nrad and Nsec!\n");
		}
		const double cps = parameters::cps;
		switch (parameters::radial_grid_type) {
			case parameters::arithmetic_spacing:
				NRadial = std::round(cps*(RMAX - RMIN) / H);
				NAzimuthal = std::round(2*M_PI/(RMAX-RMIN)*NRadial);
				break;
			case parameters::logarithmic_spacing:
				NRadial = std::round(std::log(RMAX/(double )RMIN) / std::log(1 + H/cps));
				NAzimuthal = std::round(2*M_PI/(std::pow(RMAX/(double )RMIN, 1.0/(double) NRadial) - 1));
				break;
			default:
				die("Setting resolution is not supported for the selected radial grid spacing.");
		}
		logging::print_master( LOG_INFO
		"Grid resolution set using cps = %f\n", cps);
	}

	dphi = 2.0 * M_PI / (double)NAzimuthal;
    invdphi = (double)NAzimuthal / (2.0 * M_PI);

	double cpsrad;
	switch (parameters::radial_grid_type) {
	case parameters::arithmetic_spacing:
		cpsrad = H/ ((RMAX-RMIN)/(double) NRadial);
		break;
	case parameters::logarithmic_spacing:
		cpsrad = H/(std::pow(RMAX/((double) RMIN), 1.0/((double) NRadial)) - 1);
		break;
	default:
		cpsrad = -1;
	}
	const double cpsaz = NAzimuthal * H / (2*M_PI);


	logging::print_master(
		LOG_INFO
		"The grid has (Nrad, Naz) = (%u, %u) cells with (%f, %f) cps.\n",
		NRadial, NAzimuthal, cpsrad, cpsaz);


	if(parameters::klahr_smoothing_radius == 0.0){ // otherwise Ofast will optimize out zero checks and then divide by zero
		parameters::klahr_smoothing_radius = 0.001 * (RMAX - RMIN)/(double)NRadial;
	}



    if ((parameters::radial_grid_type == parameters::logarithmic_spacing) ||
	(parameters::radial_grid_type == parameters::exponential_spacing)) {
	double c = log(RMAX / RMIN);
	double optimal_N_azimuthal = M_PI / ((std::exp(c / NRadial) - 1.0) /
					     (std::exp(c / NRadial) + 1.0));

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

    switch (parameters::ASPECTRATIO_MODE) {
    case 0:
	logging::print_master(
	    LOG_INFO
	    "Computing scale height with respect to primary object.\n");
	break;
    case 1:
	logging::print_master(
	    LOG_INFO "Computing scale height with respect to nbody system.\n");
	break;
    case 2:
	logging::print_master(
	    LOG_INFO
	    "Computing scale height with respect to center of mass.\n");
	break;
    default:
	logging::print_master(
	    LOG_INFO
	    "Computing scale height with respect to primary object.\n");
    }

    if (!cfg.contains("OuterBoundary")) {
	logging::print_master(LOG_ERROR
			      "OuterBoundary doesn't exist. Old parameter file?\n");
	die("died for convenience ;)");
    }

	ECC_GROWTH_MONITOR = config::cfg.get_flag("WriteEccentricityChange", "no");
	ecc_old = 0.0;
	peri_old = 0.0;
	delta_ecc_source = 0.0;
	delta_peri_source = 0.0;
	delta_ecc_art_visc = 0.0;
	delta_peri_art_visc = 0.0;
	delta_ecc_visc = 0.0;
	delta_peri_visc = 0.0;
	delta_ecc_transport = 0.0;
	delta_peri_transport = 0.0;
	delta_ecc_damp = 0.0;
	delta_peri_damp = 0.0;


    // Frame settings
    parameters::corotating = false;
    GuidingCenter = false;
    switch (cfg.get_first_letter_lowercase("Frame", "Fixed")) {
    case 'f': // Fixed
	break;
    case 'c': // Corotating
	parameters::corotating = true;
	break;
    case 'g': // Guiding-Center
	parameters::corotating = false;
	GuidingCenter = false;
	break;
    default:
	die("Invalid setting for Frame");
    }
    parameters::OMEGAFRAME = cfg.get<double>("OMEGAFRAME", 0);

    // Barycenter mode
    switch (cfg.get_first_letter_lowercase("HydroFrameCenter", "primary")) {
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
	    cfg.get<std::string>("HydroFrameCenter", "primary"));
    }

    if (parameters::n_bodies_for_hydroframe_center != 1 &&
	parameters::ASPECTRATIO_MODE == 1) {
	logging::print_master(
	    LOG_INFO
	    "WARNING: MORE THAN 1 CENTRAL OBJECT AND NBODY ASPECTRATIO IS NOT TESTED OR DEBUGGED!\n");
    }

    parameters::exitOnDeprecatedSetting(
	"IndirectTerm",
	"Indirect terms are now handled automatically to avoid unphysical settings.",
	"Please remove the setting.");
    parameters::exitOnDeprecatedSetting(
	"IndirectTermPlanet",
	"Indirect terms are now handled automatically to avoid unphysical settings.",
	"Please remove the setting.");

    /// EoS can only be read after apply_units(), because fitting
    /// polytropic constants requires parameters::sigma0 to be in code units.
    // Energy equation / Adiabatic
    char Adiabatic_deprecated =
	cfg.get_first_letter_lowercase("Adiabatic", "false");

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

	parameters::ADIABATICINDEX =
	    cfg.get<double>("AdiabaticIndex", 7.0 / 5.0);
	if ((parameters::Adiabatic) && (parameters::ADIABATICINDEX == 1)) {
	    logging::print_master(
		LOG_WARNING
		"You cannot have Adiabatic=true and AdiabatcIndex = 1. I decided to put Adiabatic=false, to  simulate a locally isothermal equation of state. Please check that it what you really wanted to do!\n");
	    parameters::Adiabatic = false;
	}
    } else {
	char eos_string[512];
	strncpy(
	    eos_string,
	    cfg.get<std::string>("EquationOfState", "Isothermal").c_str(),
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
	    parameters::ADIABATICINDEX =
		cfg.get<double>("AdiabaticIndex", 7.0 / 5.0);
	    logging::print_master(
		LOG_INFO
		"Using isothermal equation of state. AdiabaticIndex = %.3f.\n",
		parameters::ADIABATICINDEX);
	}
	if (strcmp(eos_string, "adiabatic") == 0 ||
	    strcmp(eos_string, "ideal") == 0) {
	    could_read_eos = true;

	    // Energy equation / Adiabatic
	    parameters::Adiabatic = true;

	    char ADIABATICINDEX_string[512];
	    strncpy(
		ADIABATICINDEX_string,
		cfg.get<std::string>("AdiabaticIndex", "7.0/5.0").c_str(),
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
		parameters::ADIABATICINDEX = cfg.get<double>(
		    "AdiabaticIndex", 7.0 / 5.0);
	    }

	    if ((parameters::Adiabatic) && (parameters::ADIABATICINDEX == 1)) {
		logging::print_master(
		    LOG_WARNING
		    "You cannot have Adiabatic=true and AdiabatcIndex = 1. I decided to put Adiabatic=false, to simulate a locally isothermal equation of state. Please check that it what you really wanted to do!\n");
		parameters::Adiabatic = false;
	    }
	    logging::print_master(
		LOG_INFO
		"Using ideal equation of state. AdiabaticIndex = %.3f.\n",
		parameters::ADIABATICINDEX);
	}

	if (strcmp(eos_string, "pvtelaw") == 0 ||
	    strcmp(eos_string, "pvte") == 0) {
	    could_read_eos = true;

	    // Energy equation / Adiabatic
	    parameters::Adiabatic = true;
	    parameters::variableGamma = true;

	    char ADIABATICINDEX_string[512];
	    strncpy(
		ADIABATICINDEX_string,
		cfg.get<std::string>("AdiabaticIndex", "7.0/5.0").c_str(),
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
		parameters::ADIABATICINDEX = cfg.get<double>(
		    "AdiabaticIndex", 7.0 / 5.0);
	    }

	    if (parameters::ADIABATICINDEX == 1) {
		logging::print_master(
		    LOG_WARNING
		    "You cannot have PVTE EoS and AdiabatcIndex = 1.! Setting AdiabaticIndex to 7/5.\n");
		parameters::ADIABATICINDEX = 7.0 / 5.0;
	    }
	    logging::print_master(
		LOG_INFO
		"PVTE EoS: Using ideal equation of state with a variable AdiabaticIndex. Init Gamma = %g!\n",
		parameters::ADIABATICINDEX);
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
		    cfg.get<std::string>("AdiabaticIndex", "2.0").c_str(),
		    256); // same as MAXNAME from config.cpp
	    for (char *t = ADIABATICINDEX_string; *t != '\0'; ++t) {
		*t = tolower(*t);
	    }

	    if (strcmp(ADIABATICINDEX_string, "fit_isothermal") == 0 ||
		strcmp(ADIABATICINDEX_string, "fit isothermal") == 0) {
		get_polytropic_constants(K, gamma);
		parameters::ADIABATICINDEX = gamma;
	    } else {
		parameters::ADIABATICINDEX =
		    cfg.get<double>("AdiabaticIndex", 2.0);
	    }

	    char POLYTROPIC_CONSTANT_string[512];
	    strncpy(
		POLYTROPIC_CONSTANT_string,
		cfg.get<std::string>("PolytropicConstant", "12.753").c_str(),
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
		parameters::POLYTROPIC_CONSTANT = K;
	    } else {
		parameters::POLYTROPIC_CONSTANT = cfg.get<double>(
		    "PolytropicConstant", 12.753);
	    }

	    if ((parameters::Polytropic) && (parameters::ADIABATICINDEX == 1)) {
		logging::print_master(
		    LOG_WARNING
		    "You cannot have Polytropic=true and AdiabatcIndex = 1. I decided to put Polytropic=false, to simulate a locally isothermal equation of state. Please check that it what you really wanted to do!\n");
		parameters::Polytropic = false;
	    }
	    logging::print_master(
		LOG_INFO
		"Using polytropic equation of state. AdiabaticIndex = %.3f.\n",
		parameters::ADIABATICINDEX);
	}

	if (!could_read_eos)
	    die("Invalid setting for Energy Equation:   %s\n",
		cfg.get<std::string>("EquationOfState").c_str());
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

    CICPlanet = cfg.get_flag("CICPLANET", "no");

    parameters::ALPHAVISCOSITY = cfg.get<double>("ALPHAVISCOSITY", 0.0);
	parameters::VISCOSITY = cfg.get<double>("VISCOSITY", 0.0, units::L0*units::L0/units::T0);

    if (!parameters::EXPLICIT_VISCOSITY && parameters::ALPHAVISCOSITY == 0.0 &&
	(parameters::artificial_viscosity_factor == 0.0 ||
	 parameters::artificial_viscosity ==
	     parameters::artificial_viscosity_none) &&
	parameters::VISCOSITY == 0.0) {
	logging::print_master(
	    LOG_ERROR
	    "You cannot use super time-stepping without any viscosity!\n");
	PersonalExit(1);
    }

    parameters::STS_NU = cfg.get<double>("STSNU", 0.01);

    if ((parameters::ALPHAVISCOSITY != 0.0) && (parameters::VISCOSITY != 0.0)) {
	logging::print_master(LOG_ERROR "You cannot use at the same time\n");
	logging::print_master(LOG_ERROR "VISCOSITY and ALPHAVISCOSITY.\n");
	logging::print_master(LOG_ERROR
			      "Edit the parameter file so as to remove\n");
	logging::print_master(LOG_ERROR
			      "one of these variables and run again.\n");
	PersonalExit(1);
    }

    if (parameters::ALPHAVISCOSITY > 0.0 || parameters::AlphaMode > 0) {
		if (parameters::AlphaMode == 1){
			logging::print_master(LOG_INFO
			      "Smooth local alpha with a threshold of T = %.3e and alpha_cold = %.3e , alpha_Hot = %.3e\n", 
				  parameters::localAlphaThreshold, parameters::alphaCold, parameters::alphaHot);
		}
		else if (parameters::AlphaMode == 2){
			logging::print_master(LOG_INFO
			      "Switched local alpha with a threshold of T = %.3e and alpha_cold = %.3e , alpha_Hot = %.3e\n", 
				  parameters::localAlphaThreshold, parameters::alphaCold, parameters::alphaHot);
		}
		else if (parameters::AlphaMode == 3){
			logging::print_master(LOG_INFO
			      "Local alpha, based on hydrogen ionization with threshold = %.3e, alpha_cold = %.3e , alpha_Hot = %.3e\n", 
				  parameters::localAlphaThreshold, parameters::alphaCold, parameters::alphaHot);
		}
		else if (parameters::AlphaMode == 4){
			logging::print_master(LOG_INFO
			      "Local alpha based on hydrogen ionization (switched) threshold = %.3e, alpha_cold = %.3e , alpha_Hot = %.3e\n", 
				  parameters::localAlphaThreshold, parameters::alphaCold, parameters::alphaHot);
		}
		else if (parameters::AlphaMode == 5){
			logging::print_master(LOG_INFO
			      "Using local alpha after Coleman 2016. \n");
		}
		else if (parameters::AlphaMode == ALPHA_STAR_DIST_DEPENDEND){
			logging::print_master(LOG_INFO
				  "Using star dist alpha, scaling from %.3e to %.3e\n", parameters::alphaCold, parameters::alphaHot);
		}
		else {
			logging::print_master(LOG_INFO
			      "Viscosity is of alpha type with alpha = %.3e\n",
			      parameters::ALPHAVISCOSITY);
		}
    }

    if (parameters::thickness_smoothing <= 0.0) {
	logging::print_master(
	    LOG_ERROR
	    "A non-vanishing potential smoothing length is required.\n");
    }


	{
	/// Read Flux Limiter

	
	std::string flux_limiter_text = config::cfg.get<std::string>("FluxLimiter", "VanLeer");

	bool could_read_flux_limiter = false;
	if (strcmp(flux_limiter_text.c_str(), "vanleer") == 0 ||
		strcmp(flux_limiter_text.c_str(), "van") == 0 ||
		strcmp(flux_limiter_text.c_str(), "leer") == 0 ||
		strcmp(flux_limiter_text.c_str(), "vl") == 0  ||
		strcmp(flux_limiter_text.c_str(), "v") == 0) {
		logging::print_master(LOG_INFO "Using VanLeer flux limiter\n");
		flux_limiter_type = 0;
		could_read_flux_limiter = true;
	}

	if (strcmp(flux_limiter_text.c_str(), "mc") == 0 ||
		strcmp(flux_limiter_text.c_str(), "m") == 0) {
		logging::print_master(LOG_INFO "Using MC flux limiter\n");
		flux_limiter_type = 1;
		could_read_flux_limiter = true;
	}

	if(!could_read_flux_limiter){
		logging::print_master(LOG_INFO "Defaulting to VanLeer flux limiter\n");
		flux_limiter_type = 0;
	}
	}

	/// Read Viscosity Stuff
    StabilizeViscosity = cfg.get<int>("STABILIZEVISCOSITY", 0);

    if (StabilizeViscosity == 1) {
	logging::print_master(
	    LOG_INFO
	    "Using pseudo implicit viscosity to limit the viscosity update step\n");
    }
    if (StabilizeViscosity == 2) {
	logging::print_master(
	    LOG_INFO
	    "Using pseudo implicit viscosity to limit the time step size\n");
    }

	/// Read Viscosity Stuff
	StabilizeArtViscosity = config::cfg.get<double>("STABILIZEARTVISCOSITY", 0);

	if (StabilizeArtViscosity == 1) {
	logging::print_master(
		LOG_INFO
		"(WT Only) Using pseudo implicit artificial viscosity to limit the viscosity update step\n");
	}
	if (StabilizeArtViscosity == 2) {
	logging::print_master(
		LOG_INFO
		"(WT Only) Using pseudo implicit artificial viscosity to limit the time step size\n");
	}


    if (parameters::VISCOSITY != 0.0) {
	logging::print_master(
	    LOG_INFO "Viscosity is kinematic viscosity with nu = %.3e\n",
	    parameters::VISCOSITY);
    }


	if (parameters::VISCOUS_ACCRETION) {
	logging::print_master(
	    LOG_INFO
	    "VISCOUS_ACCRETION is true, recomputing viscosity before accreting mass.\n");
    }

	if (config::cfg.get_flag("WriteDefaultValues", "no")) {
		config::cfg.write_default(output::outdir + "default_config.yml");
	}

	// TODO: check with Lucas why needed?
	if(parameters::heating_star_enabled){
		data[t_data::ASPECTRATIO].set_do_before_write(nullptr);
	}

    parameters::VISCOUS_ACCRETION = false;
    if (parameters::boundary_inner == parameters::boundary_condition_viscous_outflow) {
		parameters::VISCOUS_ACCRETION = true;
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

double TellNbOutputs(double time) { return (time / parameters::DT / parameters::NINTERM); }

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
			  parameters::ASPECTRATIO_REF);
    logging::print_master(LOG_VERBOSE "VKep at inner edge    : %.3g\n",
			  std::sqrt(constants::G * 1.0 * (1. - 0.0) / RMIN));
    logging::print_master(LOG_VERBOSE "VKep at outer edge    : %.3g\n",
			  std::sqrt(constants::G * 1.0 / RMAX));
    /*
    logging::print_master(LOG_VERBOSE "boundary_inner        : %i\n",
    parameters::boundary_inner); logging::print_master(LOG_VERBOSE
    "boundary_outer        : %i\n", parameters::boundary_outer);
    */
    // temp=2.0*PI*parameters::sigma0/(2.0-parameters::SIGMASLOPE)*(pow(RMAX,2.0-parameters::SIGMASLOPE)
    // - pow(RMIN,2.0-parameters::SIGMASLOPE));	/* correct this and what follows... */
    // logging::print_master(LOG_VERBOSE "Initial Disk Mass             : %g\n",
    // temp); temp=2.0*PI*parameters::sigma0/(2.0-parameters::SIGMASLOPE)*(1.0 -
    // pow(RMIN,2.0-parameters::SIGMASLOPE)); logging::print_master(LOG_VERBOSE "Initial
    // Mass inner to r=1.0  : %g \n", temp);
    // temp=2.0*PI*parameters::sigma0/(2.0-parameters::SIGMASLOPE)*(pow(RMAX,2.0-parameters::SIGMASLOPE)
    // - 1.0); logging::print_master(LOG_VERBOSE "Initial Mass outer to r=1.0  :
    // %g \n", temp);
    logging::print_master(LOG_VERBOSE
			  "Travelling time for acoustic density waves :\n");
    temp = 2.0 / 3.0 / parameters::ASPECTRATIO_REF * (pow(RMAX, 1.5) - pow(RMIN, 1.5));
    logging::print_master(
	LOG_VERBOSE
	" * From Rmin to Rmax  : %.2g = %.2f orbits ~ %.1f outputs\n",
	temp, TellNbOrbits(temp), TellNbOutputs(temp));
    temp = 2.0 / 3.0 / parameters::ASPECTRATIO_REF * (pow(RMAX, 1.5) - pow(1.0, 1.5));
    logging::print_master(
	LOG_VERBOSE
	" * From r=1.0 to Rmax: %.2g = %.2f orbits ~ %.1f outputs\n",
	temp, TellNbOrbits(temp), TellNbOutputs(temp));
    temp = 2.0 / 3.0 / parameters::ASPECTRATIO_REF * (pow(1.0, 1.5) - pow(RMIN, 1.5));
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
			  parameters::ASPECTRATIO_REF * sqrt(constants::G * 1.0));
    logging::print_master(LOG_VERBOSE " * At outer edge      : %.3g\n",
			  parameters::ASPECTRATIO_REF * sqrt(constants::G * 1.0 / RMAX));
    logging::print_master(LOG_VERBOSE " * At inner edge      : %.3g\n",
			  parameters::ASPECTRATIO_REF * sqrt(constants::G * 1.0 / RMIN));
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
	parameters::NINTERM * parameters::DT, TellNbOrbits(parameters::NINTERM * parameters::DT));
}
