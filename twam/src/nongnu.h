#ifndef _GNU_SOURCE
#ifndef _NONGNU_H_
#define _NONGNU_H_

#include <stdarg.h>

#ifndef HAVE_ASPRINTF
int vasprintf(char **ret, const char *format, va_list ap);
int asprintf(char **ret, const char *format, ...);
#endif

#endif
#endif
