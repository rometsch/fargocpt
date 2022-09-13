#ifndef HYDRO_DT_LOGGER_H
#define HYDRO_DT_LOGGER_H

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
	void reset();

public:
	hydro_dt_logger();
	void update(const double hydro_dt);
	void write(const unsigned int coarseOutputNumber,
				const unsigned int fineOutputNumber);
	void dump();
	void read();
};

#endif // HYDRO_DT_LOGGER_H
