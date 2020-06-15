#ifndef LOWTASKS_H
#define LOWTASKS_H

#include "types.h"
#include <stdlib.h>
#include <string>

void PersonalExit(int returncode);

t_polargrid *CreatePolarGrid(unsigned int Nr, unsigned int Ns,
			     const char *name);

void MultiplyPolarGridbyConstant(t_polargrid *arraysrc, double constant);

void die(const char *err, ...);
void die_errno(const char *fmt, ...);

std::string lowercase(const std::string& s);

#endif // LOWTASKS_H
