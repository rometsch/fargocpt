#pragma once

#include <chrono>

struct hydro_dt_logger_variables {
	unsigned int N_hydro_iter_DT;
	double sum_hydro_dt;
	double sum_hydro_dt_sq;
	double min_hydro_dt;
	double max_hydro_dt;
};

class hydro_dt_logger
{

	unsigned int m_N_hydro_iter_DT;
	double m_sum_hydro_dt;
	double m_sum_hydro_dt_sq;
	double m_min_hydro_dt;
	double m_max_hydro_dt;
	unsigned int m_N_hydro_last;
	unsigned int m_N_hydro_in_last_interval;
	std::chrono::steady_clock::time_point m_realtime_last;
	void reset();

public:
	hydro_dt_logger();
	void update(const double hydro_dt);
	void write(const unsigned int coarseOutputNumber,
				const unsigned int fineOutputNumber);
	unsigned int get_N_hydro_in_last_interval() const;
	void dump();
	void read();
};
