/**
	\file fpe.c

	Floating point exeception handling
*/

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifndef __USE_GNU
#define __USE_GNU
#endif

#include <fenv.h>
#include <vector>

#include "LowTasks.h"
#include "backtrace.h"
#include "fpe.h"
#include "global.h"
#include "logging.h"
#include "output.h"


//void register_exception_handler()
// void setfpe()
// {
//  struct sigaction sigact;

//  sigact.sa_sigaction = crit_err_hdlr;
//  sigact.sa_flags = SA_RESTART | SA_SIGINFO;

//  if (sigaction(SIGSEGV, &sigact, (struct sigaction *)NULL) != 0)
//  {
//   fprintf(stderr, "error setting signal handler for %d (%s)\n",
//     SIGSEGV, strsignal(SIGSEGV));

//   exit(EXIT_FAILURE);
//  }

//  if (sigaction(SIGFPE, &sigact, (struct sigaction *)NULL) != 0)
//  {
//   fprintf(stderr, "error setting signal handler for %d (%s)\n",
//     SIGFPE, strsignal(SIGFPE));

//   exit(EXIT_FAILURE);
//  }

// }

/***
	Use GNU compiler extension to trap divide by zero
	as floating point exception
 */
void enable_trap_fpe_gnu() {
    /* Can also be any of the value below:
	    FE_INEXACT           inexact result
	    FE_DIVBYZERO         division by zero
	    FE_UNDERFLOW         result not representable due to underflow
	    FE_OVERFLOW          result not representable due to overflow
	    FE_INVALID           invalid operation
	    FE_ALL_EXCEPT        bitwise OR of all supported exceptions
    */
    feenableexcept(FE_DIVBYZERO | FE_INVALID);
}

/***
	Disable usage of 
    GNU compiler extension to trap divide by zero
	as floating point exception
 */
void disable_trap_fpe_gnu() {
    /* Can also be any of the value below:
	    FE_INEXACT           inexact result
	    FE_DIVBYZERO         division by zero
	    FE_UNDERFLOW         result not representable due to underflow
	    FE_OVERFLOW          result not representable due to overflow
	    FE_INVALID           invalid operation
	    FE_ALL_EXCEPT        bitwise OR of all supported exceptions
    */
    fedisableexcept(FE_DIVBYZERO | FE_INVALID);
}


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
	enable_trap_fpe_gnu();

    // install our signal handler
    struct sigaction sa;

    sa.sa_sigaction = &fpe_sighandler;
	//sa.sa_sigaction = &crit_err_hdlr;
	//sa.sa_sigaction = &bt_sighandler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART | SA_SIGINFO;

    sigaction(SIGFPE, &sa, NULL);
#endif
}

