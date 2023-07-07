#pragma once

#include "data.h"

namespace fld {

extern t_polargrid Ka, Kb, A, B, C, D, E, Trad, Told;

void init(const unsigned int Nrad, const unsigned int Naz);
void finalize();
void radiative_diffusion(t_data &data, const double current_time, const double dt);
void handle_output();

}