#ifndef TIMESTEP_RINGBUFFER_H
#define TIMESTEP_RINGBUFFER_H

#include "planetary_system.h"

class timestep_ringbuffer
{
public:
	timestep_ringbuffer();
	~timestep_ringbuffer();
	void init(const int len, const double factor, const double dt_start);
	void print_state();
	void reinit(const int len, const double factor, const double dt_start);
	void update(const double average_time, const double dt);
	double get_mean_dt();

private:
	int m_length;
	int m_state;
	int *m_counts;
	double *m_total_times;
	double m_dt_factor;
};

#endif // TIMESTEP_RINGBUFFER_H
