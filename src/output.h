#ifndef OUTPUT_H
#define OUTPUT_H

#include "data.h"
#include "types.h"
#include <map>
#include <string>
#include <vector>

namespace output
{

void check_free_space(t_data &data);

void write_grids(t_data &data, int index, int iter, double phystime);
void write_quantities(t_data &data, unsigned int timestep,
		      unsigned int nTimeStep, bool force_update);
void write_misc(unsigned int timestep);
void write_torques(t_data &data, unsigned int timestep, bool force_update);
void write_massflow_info(t_data &data);
void write_1D_info(t_data &data);
void write_massflow(t_data &data, unsigned int timestep);
void write_lightcurves(t_data &data, unsigned int timestep, bool force_update);
void write_coarse_time(unsigned int coarseOutputNumber,
		       unsigned int fineOutputNumber);

std::vector<double> reduce_disk_quantities(t_data &data, unsigned int timestep,
					   bool force_update);

double get_misc(unsigned int timestep, std::string variable);
double get_from_ascii_file(std::string filename, unsigned int timestep,
			   unsigned int column);
std::string get_version(std::string filename);
std::string text_file_variable_description(
    const std::map<const std::string, const int> &variables,
    const std::map<const std::string, const std::string> &units);
} // namespace output

#endif // OUTPUT_H
