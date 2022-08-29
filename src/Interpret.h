#ifndef INTERPRET_H
#define INTERPRET_H

#include "data.h"

// void var(char* name, void* ptr, int type, int necessary, char* deflt);
void ReadVariables(const std::string &filename, t_data &data, int argc, char **argv);
void PrintUsage(char *execname);
double TellNbOrbits(double time);
double TellNbOutputs(double time);
void TellEverything();

#endif // INTERPRET_H
