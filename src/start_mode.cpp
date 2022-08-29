#include "start_mode.h"
#include "LowTasks.h"
#include "global.h"
#include "logging.h"
#include "output.h"

#include "util.h"
#include <cstdint>
#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace start_mode
{

StartMode mode = mode_none;
std::int32_t restart_from = -1;
std::int32_t restart_debug = false;

static std::string &rtrim(std::string &str,
			  const std::string &chars = "\t\n\v\f\r ")
{
    str.erase(str.find_last_not_of(chars) + 1);
    return str;
}

void configure_start_mode()
{
    if (CPU_Master) {

	switch (mode) {
	case mode_start:
	    if (!std::experimental::filesystem::is_empty(OUTPUTDIR)) {
		std::string backup_path = OUTPUTDIR;
		rtrim(backup_path, "/");
		backup_path += "_bak";

		for (uint32_t i = 1;
		     std::experimental::filesystem::exists(backup_path); i++) {
		    backup_path = OUTPUTDIR;
		    rtrim(backup_path, "/");
		    backup_path += "_bak" + std::to_string(i);
		}
		logging::print_master(LOG_INFO
				      "%s is not empty, backing up as %s\n",
				      OUTPUTDIR.c_str(), backup_path.c_str());
		std::experimental::filesystem::rename(OUTPUTDIR, backup_path);
		std::experimental::filesystem::create_directory(OUTPUTDIR);
	    }
	    break;
	case mode_auto:
	    if (std::experimental::filesystem::is_empty(OUTPUTDIR)) {
		logging::print_master(
		    LOG_INFO "No output found, starting fresh simulation\n");
		mode = mode_start;
		break;
	    } else {
		mode = mode_restart;
		// continue with case mode_restart to get the last outputfile
		[[fallthrough]];
	    }
	case mode_restart:
	    if (restart_from < 0) {
		const std::string last_id = output::get_last_snapshot_id();
		snapshot_dir = OUTPUTDIR + "/snapshots/" + last_id;
		restart_from = output::get_latest_output_num(last_id);
	    } else {
		// const std::string last_id = std::to_string(restart_from);
		snapshot_dir = OUTPUTDIR + "/snapshots/" +
			       std::to_string(restart_from);
		// restart_from = output::get_latest_output_num(last_id);
	    }

	    if (restart_from < 0) {
		mode = mode_start;
	    }

	    break;
	case mode_debug:
	    if (restart_from < 0) {
		const std::string last_id = output::get_last_snapshot_id();
		restart_from = output::get_latest_output_num(last_id);
	    }

	    restart_debug = true;
	    mode = mode_restart;

	    if (restart_from < 0) {
		die("Can't restart, no valid output file debugmisc.bin found. Was the simulation run with \"DebugOutputs	YES?\"");
	    }
	    break;
	default:
	    die("Invalid start_mode");
	}
    }

    MPI_Bcast(&restart_debug, 1, MPI_INT32_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&restart_from, 1, MPI_INT32_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mode, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);

    int line_size = snapshot_dir.size();
    MPI_Bcast(&line_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!CPU_Master)
	snapshot_dir.resize(line_size);
    MPI_Bcast(const_cast<char *>(snapshot_dir.data()), line_size, MPI_CHAR, 0,
	      MPI_COMM_WORLD);
}

} // namespace start_mode
