#pragma once

#include "data.h"

void init_radialarrays(void);
void resize_radialarrays(unsigned int size);

void init_physics(t_data &data);

void init_shock_tube_test(t_data &data);
void init_PVTE_shock_tube_test(t_data &data);
void init_spreading_ring_test(t_data &data);
void init_gas_density(t_data &data);
void add_gaussian_density_ring(t_data & data);
void init_eos_arrays(t_data &data);
void init_gas_energy(t_data &data);
void add_gaussian_energy_ring(t_data &data);
void init_gas_velocities(t_data &data);
void init_secondary_disk_densities(t_data &data);
void init_secondary_disk_energies(t_data &data);
void init_secondary_disk_velocities(t_data &data);
void renormalize_sigma_and_report(t_data &data);
