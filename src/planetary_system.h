#ifndef PLANETARY_SYSTEM_H
#define PLANETARY_SYSTEM_H

#include "planet.h"
#include "rebound/rebound.h"
#include <vector>

class t_planetary_system
{
  private:
    // list of all planets
    std::vector<t_planet *> m_planets;
	bool planet_restart_legacy;

  public:
    struct reb_simulation *m_rebound;
    t_planetary_system();
    ~t_planetary_system();

    inline void add_planet(t_planet *planet)
    {
	int file_id_corrector = 0;
	if(planet_restart_legacy){
		file_id_corrector = 1;
	}
	m_planets.push_back(planet);
	planet->set_planet_number(get_number_of_planets()-file_id_corrector);
    }
    inline unsigned int get_number_of_planets(void) const
    {
	return m_planets.size();
    }
    inline t_planet &get_planet(unsigned int number) const
    {
	return *(m_planets.at(number));
    }
    // inline void delete_planet(unsigned int number) { delete
    // m_planets.at(number); m_planets.erase(m_planets.begin()+number); }

    void initialize_planet_legacy(t_planet *planet, double mass,
				  double semi_major_axis, double eccentricity,
				  double phi);
    void initialize_planet_jacobi(t_planet *planet, double mass,
				  double semi_major_axis, double eccentricity,
				  double omega, double true_anomaly);
    void initialize_planet_jacobi_adjust_first_two(
	t_planet *planet, double mass, double semi_major_axis,
	double eccentricity, double omega, double true_anomaly);

    Pair get_hydro_frame_center_position();
    Pair get_hydro_frame_center_velocity();
    double get_hydro_frame_center_mass();
    void move_to_hydro_frame_center();

    void update_global_hydro_frame_center_mass();
    void calculate_orbital_elements();

    double get_mass();
    double get_mass(unsigned int n);
    Pair get_center_of_mass();
    Pair get_center_of_mass(unsigned int n);
    Pair get_center_of_mass_velocity(unsigned int n);

    void initialize_default_star();
    void init_rebound();
	void read_from_file(char *filename);
    void list_planets();
    void rotate(double angle);
    void restart(unsigned int timestep);

    void create_planet_files();
    void write_planets(unsigned int timestep, bool big_file);

	void integrate(double time, double dt);
	void copy_data_to_rebound();
	void copy_data_from_rebound();
	void correct_velocity_for_disk_accel();
};

#endif // PLANETARY_SYSTEM_H
