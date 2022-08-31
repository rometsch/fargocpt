#ifndef PLANETARY_SYSTEM_H
#define PLANETARY_SYSTEM_H

#include "planet.h"
#include "config.h"
#include "rebound/rebound.h"
#include <vector>

class t_planetary_system
{
  private:
    // list of all planets
    std::vector<t_planet *> m_planets;

  public:
    struct reb_simulation *m_rebound;
    t_planetary_system();
    ~t_planetary_system();

    inline void add_planet(t_planet *planet)
    {
	m_planets.push_back(planet);
	planet->set_planet_number(get_number_of_planets());
    }
    inline unsigned int get_number_of_planets(void) const
    {
	return m_planets.size();
    }
    inline t_planet &get_planet(unsigned int number) const
    {
	return *(m_planets.at(number));
    }
  
    void initialize_planet_jacobi(t_planet *planet, double mass,
				  double semi_major_axis, double eccentricity,
				  double omega, double true_anomaly);
    void initialize_planet_jacobi_adjust_first_two(
	t_planet *planet, double mass, double semi_major_axis,
	double eccentricity, double omega, double true_anomaly);

    Pair get_hydro_frame_center_position() const;
    Pair get_hydro_frame_center_velocity() const;
    double get_hydro_frame_center_mass() const;
    void move_to_hydro_frame_center();

    void update_global_hydro_frame_center_mass();
    void calculate_orbital_elements();

    double get_mass() const;
    double get_mass(unsigned int n) const;
    Pair get_center_of_mass() const;
    Pair get_center_of_mass(unsigned int n) const;
    Pair get_center_of_mass_velocity(unsigned int n) const;
    Pair get_center_of_mass_velocity() const;

    void init_rebound();
    void init_system(const std::string &filename);
    void config_consistency_checks();
    void init_corotation_body();
    void init_hydro_frame_center();
    void derive_config();

    void init_planet(config::Config& cfg);
    void list_planets();
    void rotate(double angle);
    void restart();

    void create_planet_files();
    void write_planets(int file_type);

    void integrate(double time, double dt);
    void copy_data_to_rebound();
    void copy_data_from_rebound();
    void correct_velocity_for_disk_accel();
    void correct_planet_accretion();
    void compute_dist_to_primary();
    void init_roche_radii();
    void update_roche_radii();
};

#endif // PLANETARY_SYSTEM_H
