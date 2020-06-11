#ifndef LOGGING_H
#define LOGGING_H

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

int log(const char *fmt, ...);
int log_master(const char *fmt, ...);


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

 int error(const char* format, ...);
 int warning(const char* format, ...);
 int notice(const char* format, ...);
 int info(const char* format, ...);
 int verbose(const char* format, ...);
 int debug(const char* format, ...);

 int error_master(const char* format, ...);
 int warning_master(const char* format, ...);
 int notice_master(const char* format, ...);
 int info_master(const char* format, ...);
 int verbose_master(const char* format, ...);
 int debug_master(const char* format, ...);

void print_runtime_info(unsigned int output_number,
			unsigned int time_step_coarse, double dt);
void print_runtime_final();
void start_timer();

} // namespace logging

#endif // LOGGING_H
