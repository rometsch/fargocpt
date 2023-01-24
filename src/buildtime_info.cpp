#include "buildtime_info.h"
#include "logging.h"

void print_buildtimeinfo() {
	    // print some information about program
    logging::print_master(LOG_INFO "fargo: This file was compiled on %s, %s.\n",
			  __DATE__, __TIME__);
#ifdef GIT_COMMIT
    logging::print_master(LOG_INFO "fargo: Last git commit: %s\n", GIT_COMMIT);
#endif
#ifdef GIT_CHANGED
    logging::print_master(
	LOG_INFO "fargo: Files changed since git commit: %s\n", GIT_CHANGED);
#endif
#ifdef _GNU_SOURCE
    logging::print_master(LOG_INFO
			  "fargo: This version of FARGO used _GNU_SOURCE\n");
#endif
#ifdef NDEBUG
    logging::print_master(
	LOG_INFO
	"fargo: This version of FARGO used NDEBUG. So no assertion checks!\n");
#else
    logging::print_master(
	LOG_INFO
	"fargo: This version of FARGO does assertion checks! Compile with NDEBUG to speed up!\n");
#endif

}