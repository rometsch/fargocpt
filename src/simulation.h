#ifndef SIMULATION_H
#define SIMULATION_H

#include "data.h"
#include "nbody/planet.h"
#include "types.h"
#include "global.h"

namespace sim {

extern double hydro_dt;
extern double last_dt;
extern double dtemp;

void init(t_data &data);
void run(t_data &data);

}

#endif // SIMULATION_H
