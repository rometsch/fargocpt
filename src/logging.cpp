/**
	\file logging.cpp
	\author Tobias Mueller <Tobias_Mueller@twam.info>
*/

#include "logging.h"

#include "global.h"
#include "output.h"
#include "parameters.h"
#include "simulation.h"
#include <chrono>
#include <cmath>

namespace logging
{

/// print timestamps?
char time_format = 0;

/// messages with level <= print_level are printed (either to stderr or stdout,
/// depending on error_level)
char print_level = 3;

/// messages with level <= error_level are printed to stderr
char error_level = 0;

// time keeping
std::chrono::steady_clock::time_point realtime_start;
std::chrono::steady_clock::time_point realtime_last_log;
unsigned int n_last_log;

int vprint(const char *fmt, va_list args)
{
    char time_buf[80];

    // default log level for messages is info
    int current_level = 3;

    // get specific level for this message
    if (fmt[0] == '<') {
	unsigned char c = fmt[1];
	if (c && fmt[2] == '>') {
	    if (c >= '0' && c <= '5') {
		current_level = c - '0';
		fmt += 3;
	    }
	}
    }

    if (current_level <= print_level) {
	if (time_format) {
	    time_t ti;
	    struct tm *ts;

	    time(&ti);
	    switch (time_format) {
	    case 1: // print timestamp
		sprintf(time_buf, "%i", (int)ti);
		break;
	    case 2: // print UTC time
		ts = gmtime(&ti);
		strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", ts);
		break;
	    case 3: // print local time
		ts = localtime(&ti);
		strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S %Z",
			 ts);
		break;
	    }
	}

	char *buf;
	int res = vasprintf(&buf, fmt, args);
	if (!time_format) {
	    fprintf(current_level <= error_level ? stderr : stdout, "[%0*i] %s",
			(int)(std::log(CPU_Number) / std::log(10) + 1), CPU_Rank, buf);
	} else {
	    fprintf(current_level <= error_level ? stderr : stdout,
			"[%0*i %s] %s", (int)(std::log(CPU_Number) / std::log(10) + 1),
		    CPU_Rank, time_buf, buf);
	}
	free(buf);
	return res;
    }

    return 0;
}

int print(const char *fmt, ...)
{
    va_list args;
    int r;

    va_start(args, fmt);
    r = vprint(fmt, args);
    va_end(args);

    return r;
}

int print_master(const char *fmt, ...)
{
    if (!CPU_Master)
	return 0;

    va_list args;
    int r;

    va_start(args, fmt);
    r = vprint(fmt, args);
    va_end(args);

    return r;
}

void start_timer()
{
    realtime_start = std::chrono::steady_clock::now();
    realtime_last_log = realtime_start;
}


void print_runtime_final()
{
    std::chrono::steady_clock::time_point realtime_end =
	std::chrono::steady_clock::now();
    double realtime = std::chrono::duration_cast<std::chrono::microseconds>(
			  realtime_end - realtime_start)
			  .count();

    double time_per_step_ms = 0.0;
    if (sim::N_hydro_iter != 0) {
	time_per_step_ms = realtime / (1000.0 * sim::N_hydro_iter);
    }
    logging::print_master(
	LOG_INFO
	"-- Final: Total Hydrosteps %d, Physical Time %.2f, Realtime %.2f seconds, Time per Step: %.2f milliseconds\n",
	sim::N_hydro_iter, sim::PhysicalTime, realtime / 1000000.0, time_per_step_ms);
}

void print_runtime_info()
{
    // Print a line with information about the runtime: current hyrdro step,
    // average runtime, ... depending on whether enough real time or number of
    // hydro steps have passed since the last log

    std::chrono::steady_clock::time_point realtime_now;
    double realtime = 0.0;
    double realtime_since_last = 0.0;

    if (parameters::log_after_real_seconds > 0.0) {
	// need to get corrent time anyways
	realtime_now = std::chrono::steady_clock::now();
	realtime_since_last =
	    std::chrono::duration_cast<std::chrono::microseconds>(
		realtime_now - realtime_last_log)
		.count();
    }

    // Do we have to log because enough steps passed?
    bool log_bc_steps =
	parameters::log_after_steps > 0 &&
	(sim::N_hydro_iter - n_last_log) >= parameters::log_after_steps;
    // Do we have to log because enough real time passed?
    bool log_bc_time =
	parameters::log_after_real_seconds > 0 &&
	realtime_since_last / 1000000.0 > parameters::log_after_real_seconds;
    if (log_bc_steps || log_bc_time) {
	if (log_bc_steps) {
	    // get current time if not happend already
	    realtime_now = std::chrono::steady_clock::now();
	    realtime_since_last =
		std::chrono::duration_cast<std::chrono::microseconds>(
		    realtime_now - realtime_last_log)
		    .count();
	}
	realtime = std::chrono::duration_cast<std::chrono::microseconds>(
		       realtime_now - realtime_start)
		       .count();
	double time_per_step_ms = 0.0;
	if ((sim::N_hydro_iter - n_last_log) != 0) {
	    time_per_step_ms =
		realtime_since_last / (1000.0 * (sim::N_hydro_iter - n_last_log));
	}

	logging::print_master(
	    LOG_INFO
	    "Logging info: snapshot %d, monitor %d, hydrostep %d, time inside simulation %f, dt %.3e, realtime %.2f s, timeperstep %.2f ms\n",
	    sim::N_snapshot, sim::N_monitor, sim::N_hydro_iter, sim::PhysicalTime, sim::last_dt,
	    realtime / 1000000.0, time_per_step_ms);

	n_last_log = sim::N_hydro_iter;
	realtime_last_log = realtime_now;
    }
}

} // namespace logging
