#ifndef PLANETARY_SYSTEM_H
#define PLANETARY_SYSTEM_H

#include "planet.h"
#include <vector>

class t_planetary_system
{
  private:
    // list of all planets
    std::vector<t_planet *> m_planets;

  public:
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
    // inline void delete_planet(unsigned int number) { delete
    // m_planets.at(number); m_planets.erase(m_planets.begin()+number); }

    void read_from_file(char *filename);
    void list_planets();
    void rotate(double angle);
    void restart(unsigned int timestep);

    void create_planet_files();
    void write_planets(unsigned int timestep, bool big_file);
};

#endif // PLANETARY_SYSTEM_H
