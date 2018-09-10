#ifndef _GNU_SOURCE

#include "nongnu.h"
#include <stdlib.h>

#ifndef HAVE_ASPRINTF

/* found on http://unixpapa.com/incnote/stdio.html */
int vasprintf(char **ret, const char *format, va_list ap)
{
	va_list ap2;
	
	int len= 100;        /* First guess at the size */
	
	if ((*ret= (char *)malloc(len)) == NULL)
		return -1;
	
	while (1) {
		int nchar;
		va_copy(ap2, ap);
		nchar = vsnprintf(*ret, len, format, ap2);

		if (nchar > -1 && nchar < len)
			return nchar;
		
		if (nchar > len) {
			len= nchar+1;
		} else {
			len*= 2;
		}
		
		if ((*ret= (char *)realloc(*ret, len)) == NULL) {
			free(*ret);
			return -1;
		}
	}
}

int asprintf(char **ret, const char *format, ...)
{
	va_list ap;
	int nc;
	va_start (ap, format);
	nc = vasprintf(ret, format, ap);
	va_end(ap);
	return nc;
}
#endif /*HAVE_ASPRINTF*/

#endif
