#include <cxxabi.h>
#include <execinfo.h>
#include <iostream>
#include <signal.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ucontext.h>
#include <unistd.h>

#include "LowTasks.h"
#include "backtrace.h"
#include "global.h"
#include "logging.h"

/**
	prints backtrace for debugging
	code from https://stackoverflow.com/a/1925461
*/
void PrintTrace()
{

	
	std::cerr << std::endl;
	std::cerr << std::endl;
	std::cerr <<"Backtrace for your convenience:" << std::endl;
	std::cerr << "----------------------------------------------------------------------" << std::endl;

    void *array[BACKTRACE_MAXDEPTH];
    int size = backtrace(array, BACKTRACE_MAXDEPTH);
    char **messages = backtrace_symbols(array, size);

    // skip first stack frame (points here)
    for (int i = 2; i < size && messages != NULL; ++i) {
	char *mangled_name = 0, *offset_begin = 0, *offset_end = 0;

	// find parantheses and +address offset surrounding mangled name
	for (char *p = messages[i]; *p; ++p) {
	    if (*p == '(') {
		mangled_name = p;
	    } else if (*p == '+') {
		offset_begin = p;
	    } else if (*p == ')') {
		offset_end = p;
		break;
	    }
	}

	// if the line could be processed, attempt to demangle the symbol
	if (mangled_name && offset_begin && offset_end &&
	    mangled_name < offset_begin) {
	    *mangled_name++ = '\0';
	    *offset_begin++ = '\0';
	    *offset_end++ = '\0';

	    int status;
	    char *real_name = abi::__cxa_demangle(mangled_name, 0, 0, &status);

	    // if demangling is successful, output the demangled function name
	    if (status == 0) {
		std::cerr << "[" << CPU_Rank << "] [bt]: (" << i << ") "
			  << messages[i] << " : " << real_name << "+"
			  << offset_begin << offset_end << std::endl;

	    }
	    // otherwise, output the mangled function name
	    else {
		std::cerr << "[" << CPU_Rank << "] [bt]: (" << i << ") "
			  << messages[i] << " : " << mangled_name << "+"
			  << offset_begin << offset_end << std::endl;
	    }
	    free(real_name);
	}
	// otherwise, print the whole line
	else {
	    std::cerr << "[" << CPU_Rank << "] [bt]: (" << i << ") "
		      << messages[i] << std::endl;
	}
    }
    free(messages);
}
