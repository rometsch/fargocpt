#pragma once

#include <cstdio>

void init_cell_finder(const double cell_growth_factor,
		      const double first_cell_size);
int get_rmed_id(const double r);
int get_rinf_id(const double r);

int get_inf_azimuthal_id(const double phi);
int get_med_azimuthal_id(const double phi);

unsigned int clamp_r_id_to_rmed_grid(const int cell_id, const bool is_vector);
unsigned int clamp_r_id_to_radii_grid(int cell_id, const bool is_vector);

unsigned int clamp_phi_id_to_grid(int cell_id);
