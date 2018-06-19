#ifndef _FPE_H_
#define _FPE_H_

#include <signal.h>

void fpe_sighandler(int sig, siginfo_t *info, void *secret);
void setfpe();

#endif
