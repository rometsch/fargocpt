#ifndef _GNU_SOURCE
#pragma once

#include <stdarg.h>

#ifndef HAVE_ASPRINTF
int vasprintf(char **ret, const char *format, va_list ap);
int asprintf(char **ret, const char *format, ...);
#endif

#endif
