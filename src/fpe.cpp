/**
	\file fpe.c

	Floating point exeception handling
*/

#include <fenv.h>
#include <vector>

#include "LowTasks.h"
#include "backtrace.h"
#include "fpe.h"
#include "global.h"
#include "logging.h"
#include "output.h"

/**
	Handles a floating point exception.
*/
void fpe_sighandler(int sig, siginfo_t *info, void *secret)
{
    const char *type;
    (void)sig;
    (void)secret;

    switch (info->si_code) {
    case FPE_INTDIV:
	type = "integer divide by zero";
	break;

    case FPE_INTOVF:
	type = "integer overflow";
	break;

    case FPE_FLTDIV:
	type = "floating-point divide by zero";
	break;

    case FPE_FLTOVF:
	type = "floating-point overflow";
	break;

    case FPE_FLTUND:
	type = "floating-point underflow";
	break;

    case FPE_FLTRES:
	type = "floating-point inexact result";
	break;

    case FPE_FLTINV:
	type = "floating-point invalid operation";
	break;

    case FPE_FLTSUB:
	type = "subscript out of range";
	break;

    default:
	type = "undefined";
	break;
    }

    logging::print(
	LOG_ERROR
	"CPU #%d has been signaled with a floating point exception of type '%s'.\n",
	CPU_Rank, type);

    PrintTrace();

    // cause abnormal process termination
    abort();
}

/**
	Initializes floating point exeception handling.
*/
void setfpe()
{
#ifdef _TRAP_FPE
    /* Can also be any of the value below:
	    FE_INEXACT           inexact result
	    FE_DIVBYZERO         division by zero
	    FE_UNDERFLOW         result not representable due to underflow
	    FE_OVERFLOW          result not representable due to overflow
	    FE_INVALID           invalid operation
	    FE_ALL_EXCEPT        bitwise OR of all supported exceptions
    */
    feenableexcept(FE_DIVBYZERO | FE_INVALID);

    // install our signal handler
    struct sigaction sa;

    sa.sa_sigaction = &fpe_sighandler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART | SA_SIGINFO;

    sigaction(SIGFPE, &sa, NULL);
#endif
}
