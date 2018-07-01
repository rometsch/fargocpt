#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "config.h"
#include "Interpret.h"
#include "LowTasks.h"
#include "constants.h"
#include "logging.h"
#include "global.h"
#include "units.h"
#include "parameters.h"

// energy euations
bool Adiabatic = false;
// frame
int Corotating, GuidingCenter;

int FastTransport, Indirect_Term, Indirect_Term_Planet;
int OuterSourceMass, CICPlanet;

void ReadVariables(char* filename, t_data &data)
{
	// read config from
	//config::read_config_from_file(filename);
	parameters::read(filename, data);


	MASSTAPER = config::value_as_double_default("MASSTAPER", 0.0000001);
	ROCHESMOOTHING = config::value_as_double_default("ROCHESMOOTHING", 0.0);
	SIGMASLOPE = config::value_as_double_default("SIGMASLOPE", 0.0);
	IMPOSEDDISKDRIFT = config::value_as_double_default("IMPOSEDDISKDRIFT", 0.0);

	FLARINGINDEX = config::value_as_double_default("FLARINGINDEX", 0.0);

	if (asprintf(&OUTPUTDIR, "%s", config::value_as_string_default("OUTPUTDIR", "out/")) < 0) {
		logging::print_master(LOG_ERROR "Not enough memory!\n");
	}

	// check if planet config exists
	if ((config::key_exists("PLANETCONFIG")) && (strlen(config::value_as_string("PLANETCONFIG"))>0)) {
		if (asprintf(&PLANETCONFIG, "%s", config::value_as_string("PLANETCONFIG")) < 0) {
			logging::print_master(LOG_ERROR "Not enough memory!\n");
		}
	} else {
		PLANETCONFIG = NULL;
	}

	OuterSourceMass = config::value_as_bool_default("OUTERSOURCEMASS", 0);

	switch (tolower(*config::value_as_string_default("TRANSPORT","Fast"))) {
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
	NTOT = config::value_as_unsigned_int_default("NTOT",1000);
	NINTERM = config::value_as_unsigned_int_default("NINTERM", 10);
	DT = config::value_as_double_default("DT",1.0);

	if ( (parameters::radial_grid_type == parameters::logarithmic_spacing) || (parameters::radial_grid_type == parameters::exponential_spacing) ) {
		double c = log(RMAX/RMIN);
		double optimal_N_azimuthal = M_PI/((exp(c/NRadial)-1.0)/(exp(c/NRadial)+1.0));

		// check if optimal azimuthal cell number differs from actual azimuthal cell number by more than 10%
		if (fabs(((double)NAzimuthal-optimal_N_azimuthal)/(double)NAzimuthal) > 0.1) {
			logging::print_master(LOG_WARNING "You have %u cells in azimuthal direction. This should be %u cells to have quadratic cells!\n",NAzimuthal,lround(optimal_N_azimuthal));
		}
	}

	// disc
	ASPECTRATIO = config::value_as_double_default("ASPECTRATIO", 0.05);

	if (!config::key_exists("OuterBoundary")) {
		logging::print_master(LOG_ERROR "OuterBoundary doesn't exist. Old .par file?\n");
		die("died for convenience ;)");
	}

	// Frame settings
	Corotating = 0;
	GuidingCenter = 0;
	switch (tolower(*config::value_as_string_default("Frame","Fixed"))) {
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
	OMEGAFRAME = config::value_as_double_default("OMEGAFRAME",0);

	Indirect_Term = config::value_as_bool_default("IndirectTerm", 1);
	Indirect_Term_Planet = config::value_as_bool_default("IndirectTermPlanet", 1);

	// Energy equation / Adiabatic
	Adiabatic = config::value_as_bool_default("Adiabatic", 0);
	ADIABATICINDEX = config::value_as_double_default("AdiabaticIndex", 7.0/5.0);
	if ( (Adiabatic) && (ADIABATICINDEX == 1) ) {
		logging::print_master(LOG_WARNING "You cannot have Adiabatic=true and AdiabatcIndex = 1. I decided to put Adiabatic=false, to simulate a locally isothermal equation of state. Please check that it what you really wanted to do!\n");
		Adiabatic = 0;
	}

	ExcludeHill = config::value_as_bool_default("EXCLUDEHILL", 0);
	CICPlanet = config::value_as_bool_default("CICPLANET", 0);

	ALPHAVISCOSITY = config::value_as_double_default("ALPHAVISCOSITY", 0.0);
	VISCOSITY = config::value_as_double_default("VISCOSITY", 0.0);

	if ((ALPHAVISCOSITY != 0.0) && (VISCOSITY != 0.0)) {
		logging::print_master(LOG_ERROR "You cannot use at the same time\n");
		logging::print_master(LOG_ERROR "VISCOSITY and ALPHAVISCOSITY.\n");
		logging::print_master(LOG_ERROR "Edit the parameter file so as to remove\n");
		logging::print_master(LOG_ERROR "one of these variables and run again.\n");
		PersonalExit(1);
	}

	if (ALPHAVISCOSITY != 0.0) {
		ViscosityAlpha = YES;
		logging::print_master(LOG_INFO "Viscosity is of alpha type\n");
	}

	if ((parameters::thickness_smoothing != 0.0) && (ROCHESMOOTHING != 0.0)) {
		logging::print_master(LOG_ERROR "You cannot use at the same time\n");
		logging::print_master(LOG_ERROR "`ThicknessSmoothing' and `RocheSmoothing'.\n");
		logging::print_master(LOG_ERROR "Edit the parameter file so as to remove\n");
		logging::print_master(LOG_ERROR "one of these variables and run again.\n");
		PersonalExit(1);
	}

	if ((parameters::thickness_smoothing <= 0.0) && (ROCHESMOOTHING <= 0.0)) {
		logging::print_master(LOG_ERROR "A non-vanishing potential smoothing length is required.\n");
		logging::print_master(LOG_ERROR "Please use either of the following variables:\n");
		logging::print_master(LOG_ERROR "`ThicknessSmoothing' *or* `RocheSmoothing'.\n");
		logging::print_master(LOG_ERROR "before launching the run again.\n");
		PersonalExit(1);
	}

	if (ROCHESMOOTHING != 0.0) {
		RocheSmoothing = YES;
		logging::print_master(LOG_INFO "Planet potential smoothing scales with their Hill sphere.\n");
	}

	// Add a trailing slash to OUTPUTDIR if needed
	if (OUTPUTDIR[strlen(OUTPUTDIR)-1] != '/') {
		unsigned int size = strlen(OUTPUTDIR);
		OUTPUTDIR = (char*)realloc(OUTPUTDIR, size+2);
		OUTPUTDIR[size] = '/';
		OUTPUTDIR[size+1] = 0;
	}

	constants::initialize_constants();

	// now we now everything to compute unit factors
	units::calculate_unit_factors();

	// TODO: This should definitely done in parameters.cpp, where values are read, but parameters::read() is called before units::calculate_unit_factors()
	// so it is not possible. Moving the read() call causes an error.
	parameters::apply_units();
}

void PrintUsage(char* execname)
{
	logging::print_master(LOG_ERROR "Usage : %s [-abcdeimnptvz] [-(0-9)] [-s number] [-f scaling] parameters file\n", execname);
	logging::print_master(LOG_ERROR "\n-a : Monitor mass and angular momentum at each timestep\n");
	logging::print_master(LOG_ERROR "-b : Adjust azimuthal velocity to impose strict centrifugal balance at t=0\n");
	logging::print_master(LOG_ERROR "-c : Sloppy CFL condition (checked at each DT, not at each timestep)\n");
	logging::print_master(LOG_ERROR "-d : Print some debugging information on 'stdout' at each timestep\n");
	logging::print_master(LOG_ERROR "-e : Activate EU test problem torque file output\n");
	logging::print_master(LOG_ERROR "-f : Scale density array by 'scaling'. Useful to increase/decrease\n");
	logging::print_master(LOG_ERROR "     disk surface density after a restart, for instance.            \n");
	logging::print_master(LOG_ERROR "-i : tabulate Sigma profile as given by restart files\n");
	logging::print_master(LOG_ERROR "-n : Disable simulation. The program just reads parameters file\n");
	logging::print_master(LOG_ERROR "-o : Overrides output directory of input file.\n");
	logging::print_master(LOG_ERROR "-p : Give profiling information at each time step\n");
	logging::print_master(LOG_ERROR "-s : Restart simulation, taking #'number' files as initial conditions\n");
	logging::print_master(LOG_ERROR "-v : Verbose mode. Tells everything about parameters file\n");
	logging::print_master(LOG_ERROR "-z : fake sequential built when evaluating sums on HD meshes\n");
	logging::print_master(LOG_ERROR "-(0-9) : only write initial (or restart) HD meshes,\n");
	logging::print_master(LOG_ERROR "     proceed to the next nth output and exit\n");
	logging::print_master(LOG_ERROR "     This option must stand alone on one switch (-va -4 is legal, -v4a is not)\n");
	PersonalExit(1);
}

double TellNbOrbits(double time)
{
	return time/2.0/PI*sqrt(constants::G*1.0/1.0/1.0/1.0);
}

double TellNbOutputs(double time)
{
	return (time/DT/NINTERM);
}

void TellEverything()
{
	double temp;

	if (!CPU_Master)
		return;

	logging::print_master(LOG_VERBOSE "Disc properties:\n");
	logging::print_master(LOG_VERBOSE "----------------\n");
	logging::print_master(LOG_VERBOSE "Inner Radius          : %g\n", RMIN);
	logging::print_master(LOG_VERBOSE "Outer Radius          : %g\n", RMAX);
	logging::print_master(LOG_VERBOSE "Aspect Ratio          : %g\n", ASPECTRATIO);
	logging::print_master(LOG_VERBOSE "VKep at inner edge    : %.3g\n", sqrt(constants::G*1.0*(1.-0.0)/RMIN));
	logging::print_master(LOG_VERBOSE "VKep at outer edge    : %.3g\n", sqrt(constants::G*1.0/RMAX));
	/*
	logging::print_master(LOG_VERBOSE "boundary_inner        : %i\n", parameters::boundary_inner);
	logging::print_master(LOG_VERBOSE "boundary_outer        : %i\n", parameters::boundary_outer);
	*/
	//temp=2.0*PI*parameters::sigma0/(2.0-SIGMASLOPE)*(pow(RMAX,2.0-SIGMASLOPE) - pow(RMIN,2.0-SIGMASLOPE));	/* correct this and what follows... */
	//logging::print_master(LOG_VERBOSE "Initial Disk Mass             : %g\n", temp);
	//temp=2.0*PI*parameters::sigma0/(2.0-SIGMASLOPE)*(1.0 - pow(RMIN,2.0-SIGMASLOPE));
	//logging::print_master(LOG_VERBOSE "Initial Mass inner to r=1.0  : %g \n", temp);
	//temp=2.0*PI*parameters::sigma0/(2.0-SIGMASLOPE)*(pow(RMAX,2.0-SIGMASLOPE) - 1.0);
	//logging::print_master(LOG_VERBOSE "Initial Mass outer to r=1.0  : %g \n", temp);
	logging::print_master(LOG_VERBOSE "Travelling time for acoustic density waves :\n");
	temp = 2.0/3.0/ASPECTRATIO*(pow(RMAX,1.5)-pow(RMIN,1.5));
	logging::print_master(LOG_VERBOSE " * From Rmin to Rmax  : %.2g = %.2f orbits ~ %.1f outputs\n", temp, TellNbOrbits(temp), TellNbOutputs(temp));
	temp = 2.0/3.0/ASPECTRATIO*(pow(RMAX,1.5)-pow(1.0,1.5));
	logging::print_master(LOG_VERBOSE " * From r=1.0 to Rmax: %.2g = %.2f orbits ~ %.1f outputs\n", temp, TellNbOrbits(temp), TellNbOutputs(temp));
	temp = 2.0/3.0/ASPECTRATIO*(pow(1.0,1.5)-pow(RMIN,1.5));
	logging::print_master(LOG_VERBOSE " * From r=1.0 to Rmin: %.2g = %.2f orbits ~ %.1f outputs\n", temp, TellNbOrbits(temp), TellNbOutputs(temp));
	temp = 2.0*PI*sqrt(RMIN*RMIN*RMIN/constants::G/1.0);
	logging::print_master(LOG_VERBOSE "Orbital time at Rmin  : %.3g ~ %.2f outputs\n", temp, TellNbOutputs(temp));
	temp = 2.0*PI*sqrt(RMAX*RMAX*RMAX/constants::G/1.0);
	logging::print_master(LOG_VERBOSE "Orbital time at Rmax  : %.3g ~ %.2f outputs\n", temp, TellNbOutputs(temp));
	logging::print_master(LOG_VERBOSE "Sound speed :\n");
	logging::print_master(LOG_VERBOSE " * At unit radius     : %.3g\n", ASPECTRATIO*sqrt(constants::G*1.0));
	logging::print_master(LOG_VERBOSE " * At outer edge      : %.3g\n", ASPECTRATIO*sqrt(constants::G*1.0/RMAX));
	logging::print_master(LOG_VERBOSE " * At inner edge      : %.3g\n", ASPECTRATIO*sqrt(constants::G*1.0/RMIN));
	logging::print_master(LOG_VERBOSE "Grid properties:\n");
	logging::print_master(LOG_VERBOSE "----------------\n");
	logging::print_master(LOG_VERBOSE "Number of (local) rings  : %d\n", NRadial);
	logging::print_master(LOG_VERBOSE "Number of (global) rings : %d\n", GlobalNRadial);
	logging::print_master(LOG_VERBOSE "Number of sectors        : %d\n", NAzimuthal);
	logging::print_master(LOG_VERBOSE "Total (local) cells      : %d\n", NRadial*NAzimuthal);
	logging::print_master(LOG_VERBOSE "Total (gobal) cells      : %d\n", GlobalNRadial*NAzimuthal);
	logging::print_master(LOG_VERBOSE "Outputs properties:\n");
	logging::print_master(LOG_VERBOSE "-------------------\n");
	logging::print_master(LOG_VERBOSE "Time increment between outputs : %.3f = %.3f orbits\n", NINTERM*DT, TellNbOrbits(NINTERM*DT));
}
