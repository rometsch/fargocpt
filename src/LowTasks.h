#ifndef _LOWTASKS_H_
#define _LOWTASKS_H_

#include <stdlib.h>
#include "types.h"

double GetGlobalIFrac(double r);

void PersonalExit(int returncode);

t_polargrid* CreatePolarGrid(unsigned int Nr, unsigned int Ns, const char* name);

void MultiplyPolarGridbyConstant(t_polargrid* arraysrc, double constant);


void die(const char *err, ...);
void die_errno(const char *fmt, ...);

#endif
