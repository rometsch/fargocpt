#ifndef SELFGRAVITY_H
#define SELFGRAVITY_H

#include "types.h"
#include "polargrid.h"

namespace selfgravity {

extern double* g_radial;
extern double* g_azimuthal;

void mpi_init(void);
void mpi_finalize(void);

void init();
void update_tobi_constants(t_data& data);

void compute(t_polargrid &density, t_polargrid &v_radial, t_polargrid &v_azimuthal, double dt, bool update);
void compute_FFT_density(t_polargrid &density);
void compute_acceleration(t_polargrid &density);
void compute_FFT_kernel();

void update_velocities(t_polargrid &v_radial, t_polargrid &v_azimuthal, double dt);

void init_azimuthal_velocity(t_polargrid &v_azimuthal);
void init_planetary_system(t_data &data);

}

#endif // SELFGRAVITY_H
