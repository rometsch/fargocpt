#include "backtrace.h"
#include "handle_signals.h"
#include "global.h"
#include "logging.h"
#include <csignal>
#include <functional>
#include "LowTasks.h"
#include <unistd.h> 

static void heartbeat()
{

    logging::print_master(LOG_INFO
			  "\nInteractive status requested with SIGUSR1\n");
    logging::print_master(LOG_INFO "hydro dt = %g\n", hydro_dt);
    logging::print_master(LOG_INFO "output number = %d\n", N_output);
    logging::print_master(LOG_INFO "outer loop iteration = %d\n", N_outer_loop);
    logging::print_master(LOG_INFO "N hydro step = %d\n", N_hydro_iter);
    logging::print_master(LOG_INFO "PhysicalTime = %g\n", PhysicalTime);
}

static void handleSIGUSR1(__attribute__((unused)) int signum) { heartbeat(); }

static void handleSIGUSR2(__attribute__((unused)) int signum) { PrintTrace(); }

static void handleSIGTERM(__attribute__((unused)) int signum) {
    write(0, "SIGTERM recieved\n", 16);
	SIGTERM_RECEIVED = 1;
}


static void registerSIGUSR1() {
    struct sigaction sa;

    sa.sa_handler = handleSIGUSR1;
    sa.sa_flags = 0; // SA_RESTART;
    sigemptyset(&sa.sa_mask);

    if (sigaction(SIGUSR1, &sa, NULL) == -1) {
        perror("sigaction SIGUSR1");
        exit(1);
    }
}

static void registerSIGUSR2() {
    struct sigaction sa;

    sa.sa_handler = handleSIGUSR2;
    sa.sa_flags = 0; // SA_RESTART;
    sigemptyset(&sa.sa_mask);

    if (sigaction(SIGUSR2, &sa, NULL) == -1) {
        perror("sigaction SIGUSR2");
        exit(1);
    }
}

static void registerSIGTERM() {
    struct sigaction sa;

    sa.sa_handler = handleSIGTERM;
    sa.sa_flags = 0; // SA_RESTART;
    sigemptyset(&sa.sa_mask);

    if (sigaction(SIGTERM, &sa, NULL) == -1) {
        perror("sigaction SIGTERM");
        exit(1);
    }
}


void register_signal_handlers()
{
    registerSIGUSR1();
    registerSIGUSR2();
	SIGTERM_RECEIVED = 0;
    registerSIGTERM();
}
