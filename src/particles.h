#ifndef _PARTICLES_H_
#define _PARTICLES_H_

#include "particle.h"
#include "data.h"

namespace particles {

extern t_particle* particles;
extern unsigned int particles_size;
extern unsigned int local_number_of_particles;

void init(void);
void calculate_accelerations_from_star_and_planets(t_data& data);
void update_velocities_from_gas_drag(t_data &data, double dt);
void update_velocities_from_disk_gravity(t_data &data, double dt);
void integrate(t_data &data, double dt);
void write(unsigned int timestep);
void move(void);

};

#endif