#ifndef CFL_H
#define CFL_H

#include "data.h"
#include "types.h"

namespace cfl {

void init(t_data &data);
double condition_cfl(t_data &data, const double dt_global_input = 0.0);
}

#endif