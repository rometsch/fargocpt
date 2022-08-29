#ifndef OUTPUT_H
#define OUTPUT_H

#include "data.h"
#include "types.h"
#include <map>
#include <string>
#include <vector>

namespace output
{


extern std::string snapshot_dir;
extern std::string last_snapshot_dir;
extern std::string outdir;

struct misc_entry {
    int timestep;
    unsigned int nTimeStep;
    double PhysicalTime;
    double OmegaFrame;
    double FrameAngle;
    double dtemp;
    double last_dt;
    unsigned int N_iter;
};

void check_free_space(t_data &data);
void cleanup_autosave();

// void write_full_output(t_data &data, const std::string &snapshot_id);
void write_full_output(t_data &data, const std::string &snapshot_id,
		       const bool register_snapshot = true);
void write_output_version();
void write_grids(t_data &data, int index, int iter, double phystime);
void write_quantities(t_data &data, bool force_update);
void write_misc();
void write_torques(t_data &data, bool force_update);
void write_massflow_info(t_data &data);
void write_1D_info(t_data &data);
void write_massflow(t_data &data, unsigned int timestep);
void write_lightcurves(t_data &data, unsigned int timestep, bool force_update);
void write_coarse_time(unsigned int coarseOutputNumber,
		       unsigned int fineOutputNumber);

std::vector<double> reduce_disk_quantities(t_data &data, unsigned int timestep,
					   bool force_update,
					   const double quantitiy_radius);

int load_misc();
std::string get_version(std::string filename);
std::string text_file_variable_description(
    const std::map<const std::string, const int> &variables,
    const std::map<const std::string, const std::string> &units);

std::string get_last_snapshot_id();
std::int32_t get_latest_output_num(const std::string &snapshot_id);

void CheckAngularMomentumConservation(t_data &data);

} // namespace output

#endif // OUTPUT_H
