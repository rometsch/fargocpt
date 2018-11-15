#ifndef FPE_H
#define FPE_H

#include <signal.h>

void fpe_sighandler(int sig, siginfo_t *info, void *secret);
void setfpe();

#endif // FPE_H
