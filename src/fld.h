#pragma once

#include "data.h"

namespace fld {

extern t_polargrid Ka, Kb, A, B, C, D, E, Trad, Told;

// Parameters
/// enable radiative diffusion
extern bool radiative_diffusion_enabled;
/// omega for SOR in radiative diffusion
extern double radiative_diffusion_omega;
/// enable automatic omega in SOR in radiative diffusion
extern bool radiative_diffusion_omega_auto_enabled;
/// maximum iterations in SOR in radiative diffusion
extern unsigned int radiative_diffusion_max_iterations;
/// tolerance for the radiative diffusion
extern double radiative_diffusion_tolerance;
// enable 2d test
extern bool radiative_diffusion_test_2d;
// constant density for 2d fld test
extern double radiative_diffusion_test_2d_density;
// constant K for 2d fld test
extern double radiative_diffusion_test_2d_K;
// number of diffusion steps to perform for the 2d test
extern unsigned int radiative_diffusion_test_2d_steps;
// enable 1d test
extern bool radiative_diffusion_test_1d;
// save internal grids of the fld module at each snapshot
extern bool radiative_diffusion_dump_data;
// check solution of linear system
extern bool radiative_diffusion_check_solution;

void init(const unsigned int Nrad, const unsigned int Naz);
void finalize();
void radiative_diffusion(t_data &data, const double current_time, const double dt);
void handle_output();
void config();
void write_logfile(const std::string filename);
}