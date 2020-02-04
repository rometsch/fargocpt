#ifndef FIND_CELL_ID_H
#define FIND_CELL_ID_H

#include "parameters.h"
#include <cstdio>


void init_cell_finder(const double cell_growth_factor, const double first_cell_size);
int get_rmed_id(const double r);
int get_rinf_id(const double r);

unsigned int get_next_azimuthal_id(const unsigned int id_low);
int get_inf_azimuthal_id(const double phi);
int get_med_azimuthal_id(const double phi);

unsigned int clamp_r_id_to_rmed_grid(int cell_id);
unsigned int clamp_r_id_to_radii_grid(int cell_id);

unsigned int clamp_phi_id_to_grid(int cell_id);

#endif // FIND_CELL_ID_H
