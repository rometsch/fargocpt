#ifndef _GNU_SOURCE
#ifndef NONGNU_H
#define NONGNU_H

#include <stdarg.h>

#ifndef HAVE_ASPRINTF
int vasprintf(char **ret, const char *format, va_list ap);
int asprintf(char **ret, const char *format, ...);
#endif

#endif // NONGNU_H
#endif
