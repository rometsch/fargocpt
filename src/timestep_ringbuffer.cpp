#include<cstdio>
#include "timestep_ringbuffer.h"
#include "logging.h"

timestep_ringbuffer::timestep_ringbuffer()
{
}

void timestep_ringbuffer::init(const int len, const double factor, const double dt_start){

	m_state = 0;
	m_length = len;
	if(m_length == 0)
	{
		return;
	}
	m_counts = new int[len];
	m_total_times = new double[len];
	m_dt_factor = factor;

	for(int i = 0; i < len; ++i){
		m_counts[i] = 0;
		m_total_times[i] = 0.0;
	}
	m_counts[m_state] = 1;
	m_total_times[m_state] = dt_start*10.0;
}

void timestep_ringbuffer::reinit(const int len, const double factor, const double dt_start){

	delete[] m_counts;
	delete[] m_total_times;

	m_state = 0;
	m_length = len;

	if(m_length == 0)
	{
		return;
	}

	m_counts = new int[len];
	m_total_times = new double[len];
	m_dt_factor = factor;

	for(int i = 0; i < len; ++i){
		m_counts[i] = 0;
		m_total_times[i] = 0.0;
	}
	m_counts[m_state] = 1;
	m_total_times[m_state] = dt_start*10.0;
}

timestep_ringbuffer::~timestep_ringbuffer(){
	delete m_counts;
	delete m_total_times;
}

void timestep_ringbuffer::print_state() const{
	logging::print(LOG_INFO "Buffer length = %d	Buffer state = %d\n", m_length, m_state);

	for(int i = 0; i < m_length; ++i){
		logging::print(LOG_INFO "Buffer %d count = %d time = %5e\n", i, m_counts[i], m_total_times[i]);
	}

}

void timestep_ringbuffer::update(const double average_time, const double dt)
{
	m_total_times[m_state] += dt;
	m_counts[m_state]++;

	if(m_total_times[m_state] > average_time){
		m_state = (m_state + 1) % m_length;
		m_counts[m_state] = 0;
		m_total_times[m_state] = 0.0;
	}

}

double timestep_ringbuffer::get_mean_dt(){

	print_state();

	if(m_length == 0 || m_dt_factor == 0.0){
		return 1.0e100; // infinity
	}

	int count = 0;
	int time = 0.0;
	for(int i = 0; i < m_length; ++i){
		count += m_counts[i];
		time += m_total_times[i];

		printf("loop (%d %d) 	time = %.3e %.3e	count = %d %d	factor = %.3e\n",
					   i, m_length,
					   time,  m_total_times[i], count, m_counts[i], m_dt_factor);
	}

	if(count == 0){
		return 1.0e100; // infinity
	}

	print_state();

	printf("Mean dt = %.3e	time = %.3e	count = %.3e	factor = %.3e\n",
				   time / (double)count * m_dt_factor,
				   time, count, m_dt_factor);

	return time / ((double)count) * m_dt_factor;
}
