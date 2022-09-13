#include "hydro_dt_logger.h"
#include<cmath>
#include "global.h"
#include "logging.h"
#include "LowTasks.h"
#include <fstream>
#include <iostream>
#include "start_mode.h"

hydro_dt_logger::hydro_dt_logger()
{
	m_N_hydro_iter_DT = 0;
	m_sum_hydro_dt = 0.0;
	m_sum_hydro_dt_sq = 0.0;
	m_min_hydro_dt = 1.0e300; // Infinity
	m_max_hydro_dt = 0.0;
}


void hydro_dt_logger::update(const double hydro_dt)
{
	m_N_hydro_iter_DT++;
	m_sum_hydro_dt += hydro_dt;
	m_sum_hydro_dt_sq += std::pow(hydro_dt, 2);

	m_min_hydro_dt = std::min(m_min_hydro_dt, hydro_dt);
	m_max_hydro_dt = std::max(m_max_hydro_dt, hydro_dt);
}

void hydro_dt_logger::reset()
{
	m_N_hydro_iter_DT = 0;
	m_sum_hydro_dt = 0.0;
	m_sum_hydro_dt_sq = 0.0;
	m_min_hydro_dt = 1.0e300;
	m_max_hydro_dt = 0.0;
}


void hydro_dt_logger::write(const unsigned int coarseOutputNumber,
							const unsigned int fineOutputNumber)
{
	FILE *fd = 0;
	char *fd_filename;
	static bool fd_created = false;

	if (CPU_Master) {

		if (asprintf(&fd_filename, "%s%s", OUTPUTDIR, "timestepLogging.dat") == -1) {
			logging::print_master(LOG_ERROR
								  "Not enough memory for string buffer.\n");
			PersonalExit(1);
		}
		// check if file exists and we restarted
		if ((start_mode::mode == start_mode::mode_restart) && !(fd_created)) {
			fd = fopen(fd_filename, "r");
			if (fd) {
				fd_created = true;
				fclose(fd);
			}
		}

		// open logfile
		if (!fd_created) {
			fd = fopen(fd_filename, "w");
		} else {
			fd = fopen(fd_filename, "a");
		}
		if (fd == NULL) {
			logging::print_master(
						LOG_ERROR "Can't write 'timestepLogging.dat' file. Aborting.\n");
			PersonalExit(1);
		}

		free(fd_filename);

		if (!fd_created) {
			// print header
			fprintf(
						fd,
						"# Time log for the hydro timestep size. Each entry averaged over one DT\n"
						"# One DT is %.18g (code) and %.18g (cgs). Time unit is: %.18g\n"
						"# Syntax: coarse output step <tab> fine output step <tab> NumHydrosteps in last DT <tab> mean dt <tab> min dt <tab> max dt <tab> std dev\n",
						DT, DT * units::time.get_cgs_factor(), units::time.get_cgs_factor());
			fd_created = true;
		}

		double mean_dt;
		double std_dev;
		if(m_N_hydro_iter_DT > 0){
		mean_dt = m_sum_hydro_dt / m_N_hydro_iter_DT;

		const double mean_sq_dt = m_sum_hydro_dt_sq / m_N_hydro_iter_DT;
		std_dev = std::sqrt(mean_sq_dt - std::pow(mean_dt, 2));
		} else {
			mean_dt = 0.0;
			std_dev = 0.0;
		}

		fprintf(fd, "%u\t%u\t%u\t%#.16e\t%#.16e\t%#.16e\t%#.16e\n", coarseOutputNumber, fineOutputNumber,
				m_N_hydro_iter_DT, mean_dt, m_min_hydro_dt, m_max_hydro_dt, std_dev);
		fclose(fd);

		reset();
	}
}

void hydro_dt_logger::dump()
{
	hydro_dt_logger_variables vars;
	vars.N_hydro_iter_DT = m_N_hydro_iter_DT;
	vars.sum_hydro_dt = m_sum_hydro_dt;
	vars.sum_hydro_dt_sq = m_sum_hydro_dt_sq;
	vars.min_hydro_dt = m_min_hydro_dt;
	vars.max_hydro_dt = m_max_hydro_dt;

	std::ofstream wf;

	char *filename = 0;
	if (asprintf(&filename, "%sdtLoggerDump.bin", OUTPUTDIR) == -1) {
		logging::print(LOG_ERROR "Not enough memory!\n");
		PersonalExit(1);
	}

	wf = std::ofstream(filename,
			   std::ios::out |
				   std::ios::binary);

	if (!wf) {
	logging::print(LOG_ERROR "Can't write %s file. Aborting.\n", filename);
	die("End\n");
	}

	wf.write((char *)(&vars), sizeof(hydro_dt_logger_variables));
	wf.close();
	free(filename);

}


void hydro_dt_logger::read()
{
	hydro_dt_logger_variables vars;

	char *filename = 0;
	if (asprintf(&filename, "%sdtLoggerDump.bin", OUTPUTDIR) == -1) {
		logging::print(LOG_ERROR "Not enough memory!\n");
		PersonalExit(1);
	}

	std::ifstream rf(filename, std::ofstream::binary | std::ios::in);

	if (!rf) {
	logging::print(LOG_ERROR "Can't read %s file. Aborting.\n", filename);
	die("End\n");
	}

	rf.read((char *)&vars, sizeof(hydro_dt_logger_variables));

	rf.close();
	free(filename);

	m_N_hydro_iter_DT = vars.N_hydro_iter_DT;
	m_sum_hydro_dt = vars.sum_hydro_dt;
	m_sum_hydro_dt_sq = vars.sum_hydro_dt_sq;
	m_min_hydro_dt = vars.min_hydro_dt;
	m_max_hydro_dt = vars.max_hydro_dt;

}
