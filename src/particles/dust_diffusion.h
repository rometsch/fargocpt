#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "../data.h"
#include "../particle.h"

namespace dust_diffusion
{
void init();
void diffuse_dust(t_data &data, std::vector<t_particle> &particles,
		  const double dt, const unsigned int N_particles);
void kick_particle(t_particle &particle, t_data &data, const double dt);
void compute_gas_diffusion_coefficient(t_data &data);
void compute_gas_density_radial_derivative(t_data &data);
} // namespace dust_diffusion

#endif // DIFFUSION_H
