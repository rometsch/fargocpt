#include <stdlib.h>
#include <stddef.h>
#include <execinfo.h>

#include "backtrace.h"
#include "logging.h"
#include "LowTasks.h"

/**
	prints backtrace for debugging
*/
void PrintTrace(void)
{
	void *array[BACKTRACE_MAXDEPTH];
	size_t size;
	char **strings;
	size_t i;

	size = backtrace(array, BACKTRACE_MAXDEPTH);
	strings = backtrace_symbols(array, size);

	logging::print(LOG_ERROR "Obtained %zd stack frames.\n", size);

	for (i = 2; i < size; i++)
		logging::print(LOG_ERROR "%s\n", strings[i]);

	free(strings);
}
