#ifndef PLANET_H
#define PLANET_H

#include <string>

class t_planet
{
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
	/// feels the disk (ie migrates)
	bool m_feeldisk;
	/// feels other planets gravity
	bool m_feelother;
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

      public:
	// setter
	inline void set_mass(double value) { m_mass = value; }
	inline void set_x(double value) { m_x = value; }
	inline void set_y(double value) { m_y = value; }
	inline void set_vx(double value) { m_vx = value; }
	inline void set_vy(double value) { m_vy = value; }
	inline void set_acc(double value) { m_acc = value; }
	void set_name(const char *value);
	inline void set_feeldisk(bool value) { m_feeldisk = value; }
	inline void set_feelother(bool value) { m_feelother = value; }
	inline void set_planet_number(unsigned int value)
	{
		m_planet_number = value;
	}
	inline void set_temperature(double value) { m_temperature = value; }
	inline void set_radius(double value) { m_radius = value; }
	inline void set_irradiate(bool value) { m_irradiate = value; }
	inline void set_rampuptime(double value) { m_rampuptime = value; }

	// getter
	inline double get_mass(void) const { return m_mass; }
	inline double get_x(void) const { return m_x; }
	inline double get_y(void) const { return m_y; }
	inline double get_vx(void) const { return m_vx; }
	inline double get_vy(void) const { return m_vy; }
	inline double get_acc(void) const { return m_acc; }
	inline const char *get_name(void) const { return m_name; }
	inline bool get_feeldisk(void) const { return m_feeldisk; }
	inline bool get_feelother(void) const { return m_feelother; }
	inline unsigned int get_planet_number(void) const
	{
		return m_planet_number;
	}
	inline double get_temperature(void) const { return m_temperature; }
	inline double get_radius(void) const { return m_radius; }
	inline double get_irradiate(void) const { return m_irradiate; }
	inline double get_rampuptime(void) const { return m_rampuptime; }

	double get_distance();
	double get_r(void) const;
	double get_phi(void) const;
	double get_semi_major_axis();
	double get_angular_momentum();
	double get_period();
	double get_omega();
	double get_eccentricity();

	void create_planet_file();
	void write(unsigned int timestep, bool big_file);
	void restart(unsigned int timestep);
	double get_value_from_file(unsigned int timestep,
				   std::string variable_name);
};

#endif // PLANET_H
