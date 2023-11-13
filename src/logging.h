#pragma once

#include <stdarg.h>
#include <chrono>
#include <string>


#define LOG_ERROR "<0>"	  /* error conditions                     */
#define LOG_WARNING "<1>" /* warning conditions                   */
#define LOG_NOTICE "<2>"  /* normal but significant condition     */
#define LOG_INFO "<3>"	  /* informational                        */
#define LOG_VERBOSE "<4>" /* verbose                              */
#define LOG_DEBUG "<5>"	  /* debug-level messages                 */

class t_data;

namespace logging
{

extern char print_level;
extern char error_level;

extern std::chrono::steady_clock::time_point realtime_start;

int vprint(const char *fmt, va_list args);
int print(const char *fmt, ...);
int print_master(const char *fmt, ...);

void print_runtime_info();
void print_runtime_final();
void start_timer();

void init_logfiles(const std::string outdir);
void finalize();

} // namespace logging
