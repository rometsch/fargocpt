#ifndef _OUTPUT_H_
#define _OUTPUT_H_

#include "types.h"
#include "data.h"
#include <vector>

namespace output {

void check_free_space(t_data &data);

void write_grids(t_data &data, int index, int iter, double phystime);
void write_quantities(t_data &data);
void write_misc(t_data &data, unsigned int timestep);
void write_disk_quantities(t_data &data, unsigned int timestep, bool force_update);
void write_torques(t_data &data, unsigned int timestep, bool force_update);
void write_lightcurves(t_data &data, unsigned int timestep, bool force_update);
void write_coarse_time(unsigned int coarseOutputNumber, unsigned int fineOutputNumber, double physicalTime);

double get_misc(unsigned int timestep, unsigned int column);

}

#endif
