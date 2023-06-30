#include <cxxabi.h>
#include <execinfo.h>
#include <iostream>
#include <signal.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ucontext.h>
#include <unistd.h>

#include "backtrace.h"
#include "global.h"
#include "output.h"

void Backtrace(std::ostream &);
void Backtrace(std::ostream &out) {
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
		out << "[" << CPU_Rank << "] [bt]: (" << i << ") "
			  << messages[i] << " : " << real_name << "+"
			  << offset_begin << offset_end << std::endl;

	    }
	    // otherwise, output the mangled function name
	    else {
		out << "[" << CPU_Rank << "] [bt]: (" << i << ") "
			  << messages[i] << " : " << mangled_name << "+"
			  << offset_begin << offset_end << std::endl;
	    }
	    free(real_name);
	}
	// otherwise, print the whole line
	else {
	    out << "[" << CPU_Rank << "] [bt]: (" << i << ") "
		      << messages[i] << std::endl;
	}
    }
    free(messages);
}

/**
	prints backtrace for debugging
	code from https://stackoverflow.com/a/1925461
*/
void PrintTrace(bool tofile)
{
	if (tofile) {
		const std::string filename = output::outdir + "backtrace.txt";
		std::ofstream out(filename);
		Backtrace(out);
		std::cerr << "[" << CPU_Rank << "] [bt]: backtrace written to " << filename << std::endl;
	} else {
		Backtrace(std::cerr);
	}
}
