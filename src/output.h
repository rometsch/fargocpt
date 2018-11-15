#ifndef OUTPUT_H
#define OUTPUT_H

#include "types.h"
#include "data.h"
#include <vector>

namespace output {

void check_free_space(t_data &data);

void write_grids(t_data &data, int index, int iter, double phystime);
void write_quantities(t_data &data);
void write_misc(unsigned int timestep);
void write_disk_quantities(t_data &data, unsigned int timestep, bool force_update);
void write_torques(t_data &data, unsigned int timestep, bool force_update);
void write_lightcurves(t_data &data, unsigned int timestep, bool force_update);
void write_coarse_time(unsigned int coarseOutputNumber, unsigned int fineOutputNumber);

double get_misc(unsigned int timestep, unsigned int column);

}

#endif // OUTPUT_H
