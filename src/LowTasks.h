#ifndef LOWTASKS_H
#define LOWTASKS_H

#include "types.h"
#include <stdlib.h>

void PersonalExit(int returncode);

t_polargrid *CreatePolarGrid(unsigned int Nr, unsigned int Ns,
			     const char *name);

void MultiplyPolarGridbyConstant(t_polargrid *arraysrc, double constant);

void die(const char *err, ...);
void die_errno(const char *fmt, ...);

#endif // LOWTASKS_H
