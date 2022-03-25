#ifndef PLANET_H
#define PLANET_H

#include "types.h"
#include <string>

struct planet_member_variables {
    unsigned int timestep;
    double m_mass;
    double m_x;
    double m_y;
    double m_vx;
    double m_vy;
    /// accretion times^-1
    double m_acc;
    double m_accreted_mass;
    unsigned int m_planet_number;
    double m_temperature;
    double m_radius;
    bool m_irradiate;
    double m_rampuptime;
    Pair m_disk_on_planet_acceleration;
    Pair m_nbody_on_planet_acceleration;

    /// orbital elements
    double m_semi_major_axis;
    double m_eccentricity;
    double m_mean_anomaly;
    double m_true_anomaly;
    double m_eccentric_anomaly;
    double m_pericenter_angle;

    double m_torque;
};

class t_planet
{

    friend class t_planetary_system;

  private:
    double m_mass;
    double m_x;
    double m_y;
    double m_vx;
    double m_vy;
    /// accretion times^-1
    double m_acc;
    double m_accreted_mass;
    std::string m_name;
    unsigned int m_planet_number;
    double m_temperature;
    double m_radius;
    bool m_irradiate;
    double m_rampuptime;
    Pair m_disk_on_planet_acceleration;
    Pair m_nbody_on_planet_acceleration;

    /// orbital elements
    double m_semi_major_axis;
    double m_eccentricity;
    double m_mean_anomaly;
    double m_true_anomaly;
    double m_eccentric_anomaly;
    double m_pericenter_angle;

    double m_torque;

  public:
	void print();
    inline void add_accreted_mass(double value) { m_accreted_mass += value; }
    // setter
    inline void set_mass(double value) { m_mass = value; }
    inline void set_x(double value) { m_x = value; }
    inline void set_y(double value) { m_y = value; }
    inline void set_vx(double value) { m_vx = value; }
    inline void set_vy(double value) { m_vy = value; }
    inline void set_acc(double value) { m_acc = value; }
    inline void set_torque(const double value) { m_torque = value; }
    void set_name(const std::string value);
    inline void set_planet_number(unsigned int value)
    {
	m_planet_number = value;
    }
    inline void set_temperature(double value) { m_temperature = value; }
	inline void set_planet_radial_extend(double value) { m_radius = value; }
    inline void set_irradiate(bool value) { m_irradiate = value; }
    inline void set_rampuptime(double value) { m_rampuptime = value; }
    inline void set_disk_on_planet_acceleration(Pair value)
    {
	m_disk_on_planet_acceleration = value;
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
	double get_rampup_mass() const;
    inline double get_x(void) const { return m_x; }
    inline double get_y(void) const { return m_y; }
    inline double get_vx(void) const { return m_vx; }
    inline double get_vy(void) const { return m_vy; }
    inline double get_acc(void) const { return m_acc; }
    inline const std::string &get_name(void) const { return m_name; }
    inline unsigned int get_planet_number(void) const
    {
	return m_planet_number;
    }
    inline double get_temperature(void) const { return m_temperature; }
	inline double get_planet_radial_extend(void) const { return m_radius; }
    inline double get_irradiate(void) const { return m_irradiate; }
    inline double get_rampuptime(void) const { return m_rampuptime; }
    inline const Pair get_disk_on_planet_acceleration(void) const
    {
	return m_disk_on_planet_acceleration;
    }
    inline const Pair get_nbody_on_planet_acceleration(void) const
    {
	return m_nbody_on_planet_acceleration;
    }

    inline double get_semi_major_axis() const { return m_semi_major_axis; }
    inline double get_eccentricity() const { return m_eccentricity; }
    inline double get_mean_anomaly() const { return m_mean_anomaly; }
    inline double get_true_anomaly() const { return m_true_anomaly; }
    inline double get_eccentric_anomaly() const { return m_eccentric_anomaly; }
    inline double get_pericenter_angle() const { return m_pericenter_angle; }
    inline double get_torque() const { return m_torque; }
    inline double get_accreted_mass() const { return m_accreted_mass; }

    double get_r(void) const;
    double get_phi(void) const;
    double get_angular_momentum() const;
	double get_orbital_period() const;
    double get_omega() const;
	double get_rhill() const;

    void calculate_orbital_elements(double x, double y, double vx, double vy,
				    double com_mass);
    void set_orbital_elements_zero();

    void copy(const planet_member_variables &other);
    void create_planet_file(bool debug_output);
    void write(const unsigned int timestep, const unsigned int file_type);
    void write_ascii(const char *filename, const unsigned int timestep) const;
	void write_binary(const char *filename, const unsigned int timestep, const bool debug) const;
    void restart(unsigned int timestep, bool debug);
	double get_value_from_file(unsigned int timestep, std::string variable_name);
    ~t_planet();
    t_planet();
};

#endif // PLANET_H
