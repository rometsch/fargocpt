#pragma once

#include "data.h"

namespace cfl {

void init(t_data &data);
double condition_cfl(t_data &data, const double dt_global_input = 0.0);
}