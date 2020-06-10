/**
	\file logging.cpp
	\author Tobias Mueller <Tobias_Mueller@twam.info>
*/

#include "logging.h"

#include "global.h"
#include "parameters.h"
#include <chrono>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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



// int myvprint(const char *fmt, va_list args)
// {
//     char time_buf[80];

//     // default log level for messages is info
//     int current_level = 3;

//     // get specific level for this message
//     if (fmt[0] == '<') {
// 	unsigned char c = fmt[1];
// 	if (c && fmt[2] == '>') {
// 	    if (c >= '0' && c <= '5') {
// 		current_level = c - '0';
// 		fmt += 3;
// 	    }
// 	}
//     }

//     if (current_level <= print_level) {
// 	if (time_format) {
// 	    time_t ti;
// 	    struct tm *ts;

// 	    time(&ti);
// 	    switch (time_format) {
// 	    case 1: // print timestamp
// 		sprintf(time_buf, "%i", (int)ti);
// 		break;
// 	    case 2: // print UTC time
// 		ts = gmtime(&ti);
// 		strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", ts);
// 		break;
// 	    case 3: // print local time
// 		ts = localtime(&ti);
// 		strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S %Z",
// 			 ts);
// 		break;
// 	    }
// 	}

// 	char *buf;
// 	int res = vasprintf(&buf, fmt, args);
// 	if (!time_format) {
// 	    fprintf(current_level <= error_level ? stderr : stdout, "[%0*i] %s",
// 		    (int)(log(CPU_Number) / log(10) + 1), CPU_Rank, buf);
// 	} else {
// 	    fprintf(current_level <= error_level ? stderr : stdout,
// 		    "[%0*i %s] %s", (int)(log(CPU_Number) / log(10) + 1),
// 		    CPU_Rank, time_buf, buf);
// 	}
// 	free(buf);
// 	return res;
//     }

//     return 0;
// }

// int print_flagged(const unsigned int flag, const char *fmt, va_list args)
// {
//     char time_buf[80];

//     // log level from flag value
//     int current_level = flag;

//     if (current_level <= print_level) {
// 	if (time_format) {
// 	    time_t ti;
// 	    struct tm *ts;

// 	    time(&ti);
// 	    switch (time_format) {
// 	    case 1: // print timestamp
// 		sprintf(time_buf, "%i", (int)ti);
// 		break;
// 	    case 2: // print UTC time
// 		ts = gmtime(&ti);
// 		strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", ts);
// 		break;
// 	    case 3: // print local time
// 		ts = localtime(&ti);
// 		strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S %Z",
// 			 ts);
// 		break;
// 	    }
// 	}

// 	char *buf;
// 	int res = vasprintf(&buf, fmt, args);
// 	if (time_format) {
// 	    fprintf(current_level <= error_level ? stderr : stdout,
// 		    "[%0*i %s] %s", (int)(log(CPU_Number) / log(10) + 1),
// 		    CPU_Rank, time_buf, buf);
// 	} else {
// 	    fprintf(current_level <= error_level ? stderr : stdout, "[%0*i] %s",
// 		    (int)(log(CPU_Number) / log(10) + 1), CPU_Rank, buf);
// 	}
// 	free(buf);
// 	return res;
//     }

//     return 0;
// }

int error_master(const char *fmt, ...);
int warning_master(const char *fmt, ...);
int notice_master(const char *fmt, ...);
int info_master(const char *fmt, ...);
int verbose_master(const char *fmt, ...);
int debug_master(const char *fmt, ...);

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
    if (N_iter != 0) {
	time_per_step_ms = realtime / (1000.0 * N_iter);
    }
    logging::info_master("-- Final: Total Hydrosteps %d, Physical Time %.2f, Realtime %.2f seconds, Time per Step: %.2f milliseconds\n",
	N_iter, PhysicalTime, realtime / 1000000.0, time_per_step_ms);
}

void print_runtime_info(unsigned int output_number,
			unsigned int time_step_coarse, double dt)
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
    bool log_bc_steps = parameters::log_after_steps > 0 &&
			(N_iter - n_last_log) >= parameters::log_after_steps;
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
	if ((N_iter - n_last_log) != 0) {
	    time_per_step_ms =
		realtime_since_last / (1000.0 * (N_iter - n_last_log));
	}
	logging::info_master("output %d, timestep %d, hydrostep %d, time inside simulation %f, dt %.3e, realtime %.2f s, timeperstep %.2f ms\n",
	    output_number, time_step_coarse, N_iter, PhysicalTime, dt,
	    realtime / 1000000.0, time_per_step_ms);
	n_last_log = N_iter;
	realtime_last_log = realtime_now;
    }
}

} // namespace logging
