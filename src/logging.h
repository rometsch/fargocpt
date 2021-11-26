#ifndef LOGGING_H
#define LOGGING_H

#include <stdarg.h>

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

int vprint(const char *fmt, va_list args);
int print(const char *fmt, ...);
int print_master(const char *fmt, ...);

void print_runtime_info(t_data& data, unsigned int output_number,
			unsigned int time_step_coarse, double dt);
void print_runtime_final();
void start_timer();

} // namespace logging

#endif // LOGGING_H
