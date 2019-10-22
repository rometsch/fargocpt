/**
	\file options.cpp
	\author Tobias Mueller <Tobias_Mueller@twam.info>

	Command line options are handled in this file.
*/

#include "options.h"
#include "LowTasks.h"
#include "global.h"
#include "logging.h"
#include "parameters.h"
#include <getopt.h>
#include <string.h>
#include <unistd.h>

extern bool parameters::Adiabatic;

namespace options
{

bool restart = false;
unsigned int restart_from;
char *parameter_file = NULL;
bool memory_usage = false;
bool disable = false;

static const char short_options[] = "-dr:s:vqbcenom";

static const struct option long_options[] = {
    {"debug", no_argument, NULL, 'd'},
    {"verbose", no_argument, NULL, 'v'},
    {"quiet", no_argument, NULL, 'q'},
    {"restart", required_argument, NULL, 'r'},
    {"seed", required_argument, NULL, 's'},
    {0, 0, 0, 0}};

void usage(int argc, char **argv)
{
    (void)argc;
    logging::print_master(
	LOG_ERROR
	"Usage: %s [options] configfile\n\n"
	"-d | --debug        Print some debugging information on 'stdout' at each timestep\n"
	"-r | --restart <n>  Restart calculation at timestep <n>\n"
	"-v | --verbose      Verbose mode. Tells everything about parameters file\n"
	"-q | --quiet        Only print errors and warnings\n"
	"-b |                Adjust azimuthal velocity to impose strict centrifugal balance at t=0\n"
	"-c |                Sloppy CFL condition (checked at each DT, not at each timestep)\n"
	"-e |                Activate EU test problem torque file output\n"
	"-n |                Disable simulation. The program just reads parameters file\n"
	"-o |                Overrides output directory of input file.\n"
	"-m |                estimate memory usage and print out\n"
	"",
	argv[0]);
}

void parse(int argc, char **argv)
{
    for (;;) {
	int index, c = 0;

	c = getopt_long(argc, argv, short_options, long_options, &index);

	if (c == EOF)
	    break;

	switch (c) {
	case 0: // getopt_long() flag
	    break;

	case 1: // no option
	    parameter_file = new char[strlen(optarg) + 1];
	    strcpy(parameter_file, optarg);
	    break;

	case 'd':
#ifdef NDEBUG
	    die("option -d for debugging is only valid if not compiled with NDEBUG");
#else
	    logging::print_level = 5;
	    debug = YES;
#endif
	    break;
	case 'r':
	    restart = true;
	    restart_from = atoi(optarg);
	    logging::print_master(LOG_INFO
				  "Restarting simulation at timestep %u\n",
				  restart_from);
	    break;
	case 'v':
	    logging::print_level = 4;
	    break;

	case 'q':
	    logging::print_level = 2;
	    break;

	case 'b':
	    CentrifugalBalance = YES;
	    break;

	case 'c':
	    SloppyCFL = YES;
	    break;

	case 'e':
	    Stockholm = YES;
	    break;

	case 'n':
	    disable = true;
	    break;

	case 'm':
	    memory_usage = true;
	    break;

	default:
	    usage(argc, argv);
	    PersonalExit(EXIT_FAILURE);
	}
    }

    // parameter file MUST be specified
    if (parameter_file == NULL) {
	usage(argc, argv);
	PersonalExit(EXIT_FAILURE);
    }
}

} // namespace options
