#include "start_mode.h"
#include "LowTasks.h"
#include "global.h"
#include "logging.h"

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
				      OUTPUTDIR, backup_path.c_str());
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
		restart_from = get_latest_output_num();
	    }

	    if (restart_from < 0) {
		die("Can't restart, no valid output file found. Check misc.dat");
	    }
	    break;
	default:
	    die("Invalid start_mode");
	}
    }
    MPI_Bcast(&restart_from, 1, MPI_INT32_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mode, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
}

static std::string get_last_line(std::ifstream &in)
{
    std::string line;
    while (in >> std::ws && std::getline(in, line))
	;

    return line;
}

std::int32_t get_latest_output_num()
{
    std::ifstream misc_file;
    std::experimental::filesystem::path path;

    path = OUTPUTDIR;
    path /= "misc.dat";

    misc_file.open(path);

    if (!misc_file.is_open()) {
	return -1;
    }

    // find last line
    std::string last_line = get_last_line(misc_file);
    misc_file.close();

    std::string first_word =
	last_line.substr(0, last_line.find_first_of(" \t"));

    if (first_word.length() < 1 || !is_number(first_word)) {
	return -1;
    } else {
	return std::stoi(first_word);
    }
}

} // namespace start_mode
