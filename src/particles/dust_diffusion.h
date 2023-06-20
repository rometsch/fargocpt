#pragma once

#include "../data.h"
#include "particle.h"

namespace dust_diffusion
{
void init(t_data &data);
void diffuse_dust(t_data &data, std::vector<t_particle> &particles,
		  const double dt, const unsigned int N_particles);
double kick_length(t_particle &particle, t_data &data, const double dt);
void compute_gas_diffusion_coefficient(t_data &data);
void compute_gas_density_radial_derivative(t_data &data);
} // namespace dust_diffusion
