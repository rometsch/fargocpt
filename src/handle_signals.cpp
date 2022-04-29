#include <csignal>
#include <functional>
#include "backtrace.h"
#include "global.h"
#include "logging.h"

void heartbeat () {

	logging::print_master(LOG_INFO "\nInteractive status requested with SIGUSR1\n");
	logging::print_master(LOG_INFO "hydro dt = %g\n", hydro_dt);
	logging::print_master(LOG_INFO "output number = %d\n", N_output);
	logging::print_master(LOG_INFO "outer loop iteration = %d\n", N_outer_loop);
	logging::print_master(LOG_INFO "N hydro step = %d\n", N_hydro_iter);
	logging::print_master(LOG_INFO "PhysicalTime = %g\n", PhysicalTime);

};

void handleSIGUSR1( int signum ) {
	heartbeat();
}

void handleSIGUSR2( int signum ) {
	PrintTrace();
}

void register_signal_handlers () {
    signal(SIGUSR1, handleSIGUSR1);
    signal(SIGUSR2, handleSIGUSR2);
}
