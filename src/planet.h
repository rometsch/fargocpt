#ifndef PLANET_H
#define PLANET_H

#include "types.h"
#include <string>

class t_planet
{

    friend class t_planetary_system;

  private:
    /// mass
    double m_mass;
    /// x position
    double m_x;
    /// y position
    double m_y;
    /// x velocity
    double m_vx;
    /// y velocity
    double m_vy;
    /// accretion times^-1
    double m_acc;
    /// name
    char *m_name;
    /// planet number
    unsigned int m_planet_number;
    /// temperature
    double m_temperature;
    /// radius
    double m_radius;
    /// irradiate
    bool m_irradiate;
    /// rampup time
    double m_rampuptime;
    /// accelerations onto planet
    Pair m_disk_on_planet_acceleration;
    Pair m_nbody_on_planet_acceleration;
    /// orbital elements
    double m_semi_major_axis;
    double m_eccentricity;
    double m_mean_anomaly;
    double m_true_anomaly;
    double m_eccentric_anomaly;
    double m_pericenter_angle;

  public:
    // setter
    inline void set_mass(double value) { m_mass = value; }
    inline void set_x(double value) { m_x = value; }
    inline void set_y(double value) { m_y = value; }
    inline void set_vx(double value) { m_vx = value; }
    inline void set_vy(double value) { m_vy = value; }
    inline void set_acc(double value) { m_acc = value; }
    void set_name(const char *value);
    inline void set_planet_number(unsigned int value)
    {
	m_planet_number = value;
    }
    inline void set_temperature(double value) { m_temperature = value; }
    inline void set_radius(double value) { m_radius = value; }
    inline void set_irradiate(bool value) { m_irradiate = value; }
    inline void set_rampuptime(double value) { m_rampuptime = value; }
    inline void set_disk_on_planet_acceleration(Pair value)
    {
	m_disk_on_planet_acceleration = value;
    }
    inline void set_disk_on_planet_acceleration_x(double value)
    {
	m_disk_on_planet_acceleration.x = value;
    }
    inline void set_disk_on_planet_acceleration_y(double value)
    {
	m_disk_on_planet_acceleration.y = value;
    }
    inline void set_nbody_on_planet_acceleration(Pair value)
    {
	m_nbody_on_planet_acceleration = value;
    }
    inline void set_nbody_on_planet_acceleration_x(double value)
    {
	m_nbody_on_planet_acceleration.x = value;
    }
    inline void set_nbody_on_planet_acceleration_y(double value)
    {
	m_nbody_on_planet_acceleration.y = value;
    }

    // getter
    inline double get_mass(void) const { return m_mass; }
    double get_rampup_mass();
    inline double get_x(void) const { return m_x; }
    inline double get_y(void) const { return m_y; }
    inline double get_vx(void) const { return m_vx; }
    inline double get_vy(void) const { return m_vy; }
    inline double get_acc(void) const { return m_acc; }
    inline const char *get_name(void) const { return m_name; }
    inline unsigned int get_planet_number(void) const
    {
	return m_planet_number;
    }
    inline double get_temperature(void) const { return m_temperature; }
    inline double get_radius(void) const { return m_radius; }
    inline double get_irradiate(void) const { return m_irradiate; }
    inline double get_rampuptime(void) const { return m_rampuptime; }
    inline Pair get_disk_on_planet_acceleration(void) const
    {
	return m_disk_on_planet_acceleration;
    }
    inline double get_disk_on_planet_acceleration_x(void) const
    {
	return m_disk_on_planet_acceleration.x;
    }
    inline double get_disk_on_planet_acceleration_y(void) const
    {
	return m_disk_on_planet_acceleration.y;
    }
    inline Pair get_nbody_on_planet_acceleration(void) const
    {
	return m_nbody_on_planet_acceleration;
    }
    inline double get_nbody_on_planet_acceleration_x(void) const
    {
	return m_nbody_on_planet_acceleration.x;
    }
    inline double get_nbody_on_planet_acceleration_y(void) const
    {
	return m_nbody_on_planet_acceleration.y;
    }

    inline double get_semi_major_axis() const { return m_semi_major_axis; }
    inline double get_eccentricity() const { return m_eccentricity; }
    inline double get_mean_anomaly() const { return m_mean_anomaly; }
    inline double get_true_anomaly() const { return m_true_anomaly; }
    inline double get_eccentric_anomaly() const { return m_eccentric_anomaly; }
    inline double get_pericenter_angle() const { return m_pericenter_angle; }

    double get_r(void) const;
    double get_phi(void) const;
    double get_angular_momentum();
    double get_period();
    double get_omega();
	double get_rhill();

    void calculate_orbital_elements(double x, double y, double vx, double vy,
				    double com_mass);
    void set_orbital_elements_zero();

    void create_planet_file();
    void write(unsigned int timestep, bool big_file);
    void restart(unsigned int timestep);
    double get_value_from_file(unsigned int timestep,
			       std::string variable_name);
	~t_planet();
};

#endif // PLANET_H
