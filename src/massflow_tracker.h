#ifndef MASSFLOW_TRACKER_H
#define MASSFLOW_TRACKER_H


#include "nbody/planetary_system.h"

class t_massflow_tracker
{
public:
	t_massflow_tracker();
	~t_massflow_tracker();
	void init(t_planetary_system &nbody_sys);
	void update_mass(const double delta_mass);
	void write_to_file();
	void read_from_file();
	void update_mass_accretion(const double dt);
	double get_mdot();

private:
	double m_delta_mass;
	double m_averaging_time;
	double m_mdot;
};

#endif // MASSFLOW_TRACKER_H
