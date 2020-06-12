#ifndef INTERPRET_H
#define INTERPRET_H

#include "data.h"

// void var(char* name, void* ptr, int type, int necessary, char* deflt);
void ReadVariables(char *filename, t_data &data, int argc, char **argv);

void interpret_disk_parameters();
void create_outputdir(const char* filename);
void copy_config_file(const char *filename, int argc, char **argv);
void interpret_transport_algo();
void interpret_output_timesteps();
void interpret_grid();
void interpret_coorindate_system();
void warn_about_indirect_term_flag();
void interpret_equation_of_state(const char* filename);
void interpret_Nbody_interactions();
void interpret_viscosity();
void interpret_boundary();
void sanitize_output_dir_string();

void PrintUsage(char *execname);
double TellNbOrbits(double time);
double TellNbOutputs(double time);
void TellEverything();

#endif // INTERPRET_H
