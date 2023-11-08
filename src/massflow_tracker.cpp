#include<cstdio>
#include <fstream>

#include "massflow_tracker.h"
#include "boundary_conditions/boundary_conditions.h"
#include "global.h"
#include "output.h"

t_massflow_tracker::t_massflow_tracker()
{
}

t_massflow_tracker::~t_massflow_tracker(){
}

void t_massflow_tracker::init(t_planetary_system& nbody_sys){

	double averaging_time = 1.0e-12; // smaller than any dt, such that alpha = 1
	if(nbody_sys.get_number_of_planets() > 1){
	const double P = nbody_sys.get_planet(boundary_conditions::rof_planet).get_orbital_period();
	averaging_time = P * boundary_conditions::rof_averaging_time;
	}
	m_delta_mass = 0.0;
	m_averaging_time = averaging_time;
	m_mdot = 0.0;
}

void t_massflow_tracker::write_to_file(){

	if(CPU_Master && boundary_conditions::rochlobe_overflow){
	std::string filename = output::snapshot_dir + "/massflow_tracker" + ".bin";
	std::ofstream fs(filename, std::ios::out | std::ios::binary);
	fs.write((char *) &m_delta_mass, sizeof(m_delta_mass));
	fs.write((char *) &m_averaging_time, sizeof(m_averaging_time));
	fs.write((char *) &m_mdot, sizeof(m_mdot));
	}
}

void t_massflow_tracker::read_from_file(){

	if(CPU_Master && boundary_conditions::rochlobe_overflow){
	std::string filename = output::snapshot_dir + "/massflow_tracker" + ".bin";
	std::ifstream fs(filename, std::ios::in | std::ios::binary);
	fs.read((char *) &m_delta_mass, sizeof(m_delta_mass));
	fs.read((char *) &m_averaging_time, sizeof(m_averaging_time));
	fs.read((char *) &m_mdot, sizeof(m_mdot));
	}
}

void t_massflow_tracker::update_mass(const double delta_mass)
{
	if(CPU_Master && boundary_conditions::rochlobe_overflow){
	m_delta_mass += delta_mass; // sum over inner ring
	}
}

void t_massflow_tracker::update_mass_accretion(const double dt){

	if(CPU_Master && boundary_conditions::rochlobe_overflow){
	if(dt > 0.0 ){
	const double alpha = std::min(dt/m_averaging_time, 1.0);
	const double mdot_last = m_mdot;
	const double mdot = m_delta_mass / dt;
	m_mdot = (1.0 - alpha) * mdot_last + alpha * mdot;
	} else {
		m_mdot = 0.0;
	}
	m_delta_mass = 0.0; // reset delta_mass
	MPI_Send(&m_mdot, 1, MPI_DOUBLE, CPU_Highest, 0, MPI_COMM_WORLD);
	}

	if (CPU_Rank == CPU_Highest && boundary_conditions::rochlobe_overflow) {
		MPI_Recv(&m_mdot, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
				 &global_MPI_Status);
	}
	return;
}

double t_massflow_tracker::get_mdot(){
	return m_mdot;
}
