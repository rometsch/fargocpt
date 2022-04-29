#include "backtrace.h"
#include "handle_signals.h"
#include "global.h"
#include "logging.h"
#include <csignal>
#include <functional>
#include "LowTasks.h"

static void heartbeat()
{

    logging::print_master(LOG_INFO
			  "\nInteractive status requested with SIGUSR1\n");
    logging::print_master(LOG_INFO "hydro dt = %g\n", hydro_dt);
    logging::print_master(LOG_INFO "output number = %d\n", N_output);
    logging::print_master(LOG_INFO "outer loop iteration = %d\n", N_outer_loop);
    logging::print_master(LOG_INFO "N hydro step = %d\n", N_hydro_iter);
    logging::print_master(LOG_INFO "PhysicalTime = %g\n", PhysicalTime);
};

static void handleSIGUSR1(__attribute__((unused)) int signum) { heartbeat(); }

static void handleSIGUSR2(__attribute__((unused)) int signum) { PrintTrace(); }

static void handleSIGTERM(__attribute__((unused)) int signum) {

	SIGTERM_RECEIVED = true;
}


void register_signal_handlers()
{
	SIGTERM_RECEIVED = false;
    signal(SIGUSR1, handleSIGUSR1);
    signal(SIGUSR2, handleSIGUSR2);

	signal(SIGTERM, handleSIGTERM);

}
