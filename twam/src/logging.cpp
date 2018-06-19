/**
	\file logging.cpp
	\author Tobias Mueller <Tobias_Mueller@twam.info>
*/

#include "logging.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "global.h"

namespace logging {

/// print timestamps?
char time_format = 0;

/// messages with level <= print_level are printed (either to stderr or stdout, depending on error_level)
char print_level = 3;

/// messages with level <= error_level are printed to stderr
char error_level = 0;

int vprint(const char* fmt, va_list args)
{
	char time_buf[80];

	// default log level for messages is info
	int current_level = 3;

	// get specific level for this message
	if (fmt[0] == '<') {
		unsigned char c = fmt[1];
		if (c && fmt[2] == '>') {
			if (c >= '0' && c<='5') {
				current_level = c - '0';
				fmt += 3;
			}
		}
	}

	if (current_level <= print_level) {
		if (time_format) {
			time_t ti;
			struct tm  *ts;

			time(&ti);
			switch (time_format) {
				case 1: // print timestamp
					sprintf(time_buf, "%i",(int)ti);
					break;
				case 2: // print UTC time
					ts = gmtime(&ti);
					strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", ts);
					break;
				case 3: // print local time
					ts = localtime(&ti);
					strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S %Z", ts);
					break;
			}
		}

		char *buf;
		int res = vasprintf(&buf, fmt, args);
		if (!time_format) {
			fprintf(current_level <= error_level ? stderr : stdout, "[%0*i] %s",(int)(log(CPU_Number)/log(10)+1), CPU_Rank, buf);
		} else {
			fprintf(current_level <= error_level ? stderr : stdout, "[%0*i %s] %s",(int)(log(CPU_Number)/log(10)+1), CPU_Rank, time_buf, buf);			
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
	if (!CPU_Master) return 0;

	va_list args;
	int r;

	va_start(args, fmt);
	r = vprint(fmt, args);
	va_end(args);

	return r;
}

}
