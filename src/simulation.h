#pragma once

#include "data.h"

namespace sim {

extern double hydro_dt;
extern double last_dt;
extern double dtemp;

extern double time, timeInitial;
extern unsigned int N_snapshot;
extern unsigned int N_monitor;
extern unsigned long int N_hydro_iter;

void init(t_data &data);
void run(t_data &data);
void handle_outputs(t_data &data);
double CalculateTimeStep(t_data &data);

}
