#ifndef FPE_H
#define FPE_H

#include <signal.h>

void fpe_sighandler(int sig, siginfo_t *info, void *secret);
void setfpe();
void enable_trap_fpe_gnu();
void disable_trap_fpe_gnu();

#endif // FPE_H
