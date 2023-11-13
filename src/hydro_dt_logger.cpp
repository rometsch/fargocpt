#include "hydro_dt_logger.h"
#include <cmath>
#include "global.h"
#include "logging.h"
#include "LowTasks.h"
#include "output.h"
#include "start_mode.h"
#include "parameters.h"
#include "simulation.h"
#include <fstream>
#include <iostream>
#include <chrono>


hydro_dt_logger::hydro_dt_logger()
{
	m_N_hydro_iter_DT = 0;
	m_N_hydro_in_last_interval = 0;
	m_sum_hydro_dt = 0.0;
	m_sum_hydro_dt_sq = 0.0;
	m_min_hydro_dt = 1.0e300; // Infinity
	m_max_hydro_dt = 0.0;
	m_N_hydro_last = 0;
	m_realtime_last = std::chrono::steady_clock::now();
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
	m_N_hydro_in_last_interval= m_N_hydro_iter_DT;
	m_N_hydro_iter_DT = 0;
	m_sum_hydro_dt = 0.0;
	m_sum_hydro_dt_sq = 0.0;
	m_min_hydro_dt = 1.0e300;
	m_max_hydro_dt = 0.0;
	m_N_hydro_last = sim::N_hydro_iter;
	m_realtime_last = std::chrono::steady_clock::now();
}


void hydro_dt_logger::write(const unsigned int coarseOutputNumber,
							const unsigned int fineOutputNumber)
{
	FILE *fd = 0;
	static bool fd_created = false;

	if (CPU_Master) {

		const std::string filename = output::outdir + "monitor/timestepLogging.dat";

		// check if file exists and we restarted
		if ((start_mode::mode == start_mode::mode_restart) && !(fd_created)) {
			fd = fopen(filename.c_str(), "r");
			if (fd) {
				fd_created = true;
				fclose(fd);
			}
		}

		// open logfile
		if (!fd_created) {
			fd = fopen(filename.c_str(), "w");
		} else {
			fd = fopen(filename.c_str(), "a");
		}
		if (fd == NULL) {
			logging::print_master(
						LOG_ERROR "Can't write 'timestepLogging.dat' file. Aborting.\n");
			PersonalExit(1);
		}

		if (!fd_created) {
			// print header
			fprintf(
						fd,
						"#version: 2\n"
						"#FargoCPT Time log for the hydro timestep size.\n"
						"#Each entry averaged over one monitor_timestep\n"
						"#One monitor_timestep is %.18g (code) and %.18g (cgs). Time unit is: %.18g\n"
						"#version: 1.1\n"
						"#variable: 0 | snapshot number | 1\n"
						"#variable: 1 | monitor number | 1\n"
						"#variable: 2 | hydrostep number | 1\n"
						"#variable: 3 | Number of Hydrosteps in last monitor_timestep | 1\n"
						"#variable: 4 | time | %s\n"
						"#variable: 5 | walltime | s\n"
						"#variable: 6 | walltime per hydrostep | ms\n"
						"#variable: 7 | mean dt | %s\n"
						"#variable: 8 | min dt | %s\n"
						"#variable: 9 | std dev dt | %s\n",
						parameters::monitor_timestep, parameters::monitor_timestep * units::time.get_code_to_cgs_factor(), units::time.get_code_to_cgs_factor(),
						units::time.get_cgs_factor_symbol().c_str(),units::time.get_cgs_factor_symbol().c_str(),units::time.get_cgs_factor_symbol().c_str(),units::time.get_cgs_factor_symbol().c_str());
			fd_created = true;
		}

		double realtime_sec = 0.0;

		const std::chrono::steady_clock::time_point realtime_now = std::chrono::steady_clock::now();
		realtime_sec = 
			std::chrono::duration_cast<std::chrono::microseconds>(realtime_now - logging::realtime_start)
			.count();
		realtime_sec /= 1000000.0; // to seconds

		const double realtime_since_last_ms =
		std::chrono::duration_cast<std::chrono::microseconds>(
		    realtime_now - m_realtime_last)
		    .count() / 1000.0; // to ms

		double time_per_step_ms = 0.0;
		const long int Nsteps = sim::N_hydro_iter - m_N_hydro_last;
		if ((Nsteps) != 0) {
			time_per_step_ms =
			realtime_since_last_ms / Nsteps;
		}

		double mean_dt;
		double std_dev;
		if(m_N_hydro_iter_DT > 0){
			mean_dt = m_sum_hydro_dt / m_N_hydro_iter_DT;

			const double mean_sq_dt = m_sum_hydro_dt_sq / m_N_hydro_iter_DT;
			std_dev = std::sqrt(std::abs(mean_sq_dt - std::pow(mean_dt, 2)));
		} else {
			mean_dt = 0.0;
			std_dev = 0.0;
		}
		fprintf(fd, "%u\t%u\t%lu\t%u\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\n", coarseOutputNumber, fineOutputNumber,
				sim::N_hydro_iter, m_N_hydro_iter_DT, sim::time, realtime_sec, time_per_step_ms, mean_dt, m_min_hydro_dt, m_max_hydro_dt, std_dev);
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

	const std::string filename  = output::outdir + "monitor/dtLoggerDump.bin";

	wf = std::ofstream(filename, std::ios::out | std::ios::binary);

	if (!wf) {
	logging::print(LOG_ERROR "Can't write %s file. Aborting.\n", filename.c_str());
	die("End\n");
	}

	wf.write((char *)(&vars), sizeof(hydro_dt_logger_variables));
	wf.close();

}


void hydro_dt_logger::read()
{
	hydro_dt_logger_variables vars;

	const std::string filename  = output::outdir + "monitor/dtLoggerDump.bin";

	std::ifstream rf(filename, std::ofstream::binary | std::ios::in);

	if (!rf) {
	logging::print(LOG_ERROR "Can't read %s file. Aborting.\n", filename.c_str());
	die("End\n");
	}

	rf.read((char *)&vars, sizeof(hydro_dt_logger_variables));

	rf.close();

	m_N_hydro_iter_DT = vars.N_hydro_iter_DT;
	m_sum_hydro_dt = vars.sum_hydro_dt;
	m_sum_hydro_dt_sq = vars.sum_hydro_dt_sq;
	m_min_hydro_dt = vars.min_hydro_dt;
	m_max_hydro_dt = vars.max_hydro_dt;

}


unsigned int hydro_dt_logger::get_N_hydro_in_last_interval() const {
	unsigned int rv = std::max(m_N_hydro_iter_DT, m_N_hydro_in_last_interval);
	return rv;
}

