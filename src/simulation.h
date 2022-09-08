#ifndef SIMULATION_H
#define SIMULATION_H

#include "data.h"
#include "nbody/planet.h"
#include "types.h"
#include "global.h"

namespace sim {

extern double hydro_dt;
extern double last_dt;

void run_simulation(t_data &data);

}

#endif // SIMULATION_H
