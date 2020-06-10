#ifndef LOGGING_H
#define LOGGING_H

#include "fmt/printf.h"
#include "global.h"

#include <stdarg.h>

#define LOG_ERROR 0   /* error conditions                     */
#define LOG_WARNING 1 /* warning conditions                   */
#define LOG_NOTICE 2  /* normal but significant condition     */
#define LOG_INFO 3    /* informational                        */
#define LOG_VERBOSE 4 /* verbose                              */
#define LOG_DEBUG 5   /* debug-level messages                 */

namespace logging
{

extern char print_level;
extern char error_level;

void v_fmt_print(const char *format, fmt::format_args args);

// template <typename... Args>
// void fmt_print(const char* format, const Args & ... args) {
//   vreport_error(format, fmt::make_format_args(args...));
// }

int print_flagged(const unsigned int flag, const char *fmt, va_list args);

template <typename... Args> void log(const char *format, const Args &... args)
{
    fmt::print("[{}] ", CPU_Rank);
    fmt::printf(format, args...);
}
template <typename... Args>
void log(const int log_level, const char *format, const Args &... args)
{
    if (log_level <= print_level) {
	log(format, args...);
    }
}

template <typename... Args>
void log_master(const char *format, const Args &... args)
{
    if (CPU_Master) {
	log(format, args...);
    }
}
template <typename... Args>
void log_master(const int log_level, const char *format, const Args &... args)
{
    if (log_level <= print_level) {
	log_master(format, args...);
    }
}

/***
 * Available logging functions:
 * 
 * logging::error
 * logging::warning
 * logging::notice
 * logging::info
 * logging::verbose
 * logging::debug
 * 
 * Each one is also available as a cpu master version (print only on master cpu):
 * 
 * logging::error_master
 * logging::warning_master
 * logging::notice_master
 * logging::info_master
 * logging::verbose_master
 * logging::debug_master
 * 
 ***/


template <typename... Args> void error(const char *format, const Args &... args)
{
    log(LOG_ERROR, format, args...);
}

template <typename... Args>
void warning(const char *format, const Args &... args)
{
    log(LOG_WARNING, format, args...);
}

template <typename... Args>
void notice(const char *format, const Args &... args)
{
    log(LOG_NOTICE, format, args...);
}

template <typename... Args> void info(const char *format, const Args &... args)
{
    log(LOG_INFO, format, args...);
}

template <typename... Args>
void verbose(const char *format, const Args &... args)
{
    log(LOG_VERBOSE, format, args...);
}

template <typename... Args> void debug(const char *format, const Args &... args)
{
    log(LOG_DEBUG, format, args...);
}

/*** 
 * Functions to log on master cpu only.
 ***/ 

template <typename... Args>
void error_master(const char *format, const Args &... args)
{
    log_master(LOG_ERROR, format, args...);
}

template <typename... Args>
void warning_master(const char *format, const Args &... args)
{
    log_master(LOG_WARNING, format, args...);
}

template <typename... Args>
void notice_master(const char *format, const Args &... args)
{
    log_master(LOG_NOTICE, format, args...);
}

template <typename... Args>
void info_master(const char *format, const Args &... args)
{
    log_master(LOG_INFO, format, args...);
}

template <typename... Args>
void verbose_master(const char *format, const Args &... args)
{
    log_master(LOG_VERBOSE, format, args...);
}

template <typename... Args>
void debug_master(const char *format, const Args &... args)
{
    log_master(LOG_DEBUG, format, args...);
}

void print_runtime_info(unsigned int output_number,
			unsigned int time_step_coarse, double dt);
void print_runtime_final();
void start_timer();

} // namespace logging

#endif // LOGGING_H
