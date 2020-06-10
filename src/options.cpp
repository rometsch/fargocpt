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
#include "start_mode.h"
#include "util.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <string>
#include <unistd.h>

// Taken from
// https://stackoverflow.com/questions/5820810/case-insensitive-string-comp-in-c
static int strcicmp(char const *a, char const *b)
{
    for (;; a++, b++) {
	int d = tolower((unsigned char)*a) - tolower((unsigned char)*b);
	if (d != 0 || !*a)
	    return d;
    }
}

namespace options
{
char *parameter_file;
bool memory_usage = false;
bool disable = false;

bool is_number(std::string s)
{
    for (std::uint32_t i = 0; i < s.length(); i++)
	if (isdigit(s[i]) == false)
	    return false;

    return true;
}

static const char short_options[] = "-dvqbcnmh";

static const struct option long_options[] = {
    {"help", no_argument, NULL, 'h'},
    {"debug", no_argument, NULL, 'd'},
    {"verbose", no_argument, NULL, 'v'},
    {"quiet", no_argument, NULL, 'q'},
    {"start", no_argument, NULL, 's'},
    {"restart", optional_argument, NULL, 'r'},
    {"auto", no_argument, NULL, 'a'},
    {0, 0, 0, 0}};

void usage(int argc, char **argv)
{
    (void)argc;
    logging::print_master(
	LOG_ERROR,
	"Usage: %s [options] start|restart <N>|auto configfile\n\n"
	"start                  Start a new simulation from scratch\n"
	"restart <N>            Restart from an old simulation output, latest if no N specified\n"
	"auto                   Same as restart if output files are present, otherwise same as start\n"
	"-d | --debug           Print some debugging information on 'stdout' at each timestep\n"
	"-v | --verbose         Verbose mode. Tells everything about parameters file\n"
	"-q | --quiet           Only print errors and warnings\n"
	"-b |                   Adjust azimuthal velocity to impose strict centrifugal balance at t=0\n"
	"-c |                   Sloppy CFL condition (checked at each DT, not at each timestep)\n"
	"-n |                   Disable simulation. The program just reads parameters file\n"
	"-m |                   estimate memory usage and3df59ec2b9a5f5a5c00028320b2a31224a3686ecprint out\n"
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
	    if (parameter_file) {
		usage(argc, argv);
		die("Input Error: Parameter file must be the last argument");
	    }
	    if (start_mode::mode == start_mode::mode_none) {
		if (!strcicmp(optarg, "start")) {
		    start_mode::mode = start_mode::mode_start;
		} else if (!strcicmp(optarg, "restart")) {
		    start_mode::mode = start_mode::mode_restart;
		} else if (!strcicmp(optarg, "auto")) {
		    start_mode::mode = start_mode::mode_auto;
		} else {
		    usage(argc, argv);
			logging::print(LOG_ERROR, "Invalid start mode '%s'\n", optarg);
		    PersonalExit(0);
		}
	    } else {
		if (start_mode::mode == start_mode::mode_restart) {
		    if (is_number(optarg)) {
			start_mode::restart_from = std::atoi(optarg);

			if (start_mode::restart_from < 0) {
			    die("Input Error: Restart file can not be negative\n");
			}
			break;
		    }
		}

		parameter_file = new char[strlen(optarg) + 1];
		strcpy(parameter_file, optarg);
	    }
	    break;

	case 'd':
#ifdef NDEBUG
	    die("Option -d for debugging is only valid if not compiled with NDEBUG");
#else
	    logging::print_level = 5;
	    debug = YES;
#endif
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

	case 'n':
	    disable = true;
	    break;

	case 'm':
	    memory_usage = true;
	    break;

	case 'h':
	    usage(argc, argv);
	    PersonalExit(EXIT_SUCCESS);
	    break;
	default:
	    usage(argc, argv);
	    PersonalExit(EXIT_FAILURE);
	}
    }

    // parameter file MUST be specified
    if (parameter_file == NULL) {
	usage(argc, argv);
	die("Input error: no parameter file specified");
    }
}

} // namespace options
