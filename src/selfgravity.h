#pragma once

#include "polargrid.h"

namespace selfgravity
{

extern double *g_radial;
extern double *g_azimuthal;

void mpi_init(void);
void mpi_finalize(void);

void init(t_data &data);
void update_sg_constants();

void compute(t_data &data, double dt, bool update);
void compute_FFT_density(t_polargrid &density);
void compute_acceleration(t_polargrid &density);
void compute_FFT_kernel();

void update_velocities(t_polargrid &v_radial, t_polargrid &v_azimuthal,
		       double dt);

void init_azimuthal_velocity(t_polargrid &v_azimuthal);
} // namespace selfgravity
