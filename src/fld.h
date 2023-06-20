#pragma once

#include "data.h"

namespace fld {

void init(const unsigned int Nrad, const unsigned int Naz);
void finalize();
void radiative_diffusion(t_data &data, const double current_time, const double dt);

}