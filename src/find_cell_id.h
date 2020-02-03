#ifndef FIND_CELL_ID_H
#define FIND_CELL_ID_H

#include "parameters.h"
#include <cstdio>


void init_cell_finder(parameters::t_radial_grid GRID, const double cell_growth_factor, const double first_cell_size);
unsigned int get_rmed_id(parameters::t_radial_grid GRID, const double r);
unsigned int get_rinf_id(parameters::t_radial_grid GRID, const double r);

unsigned int get_next_azimuthal_id(const unsigned int id_low);
unsigned int get_inf_azimuthal_id(const double phi);
unsigned int get_med_azimuthal_id(const double phi);

int clamp_id_to_grid(int cell_id);

/*
template <>
void init_cell_finder<parameters::t_radial_grid::logarithmic_spacing>(const double cell_growth_factor, const double first_cell_size);

template <>
void init_cell_finder<parameters::arithmetic_spacing>(const double cell_growth_factor, const double first_cell_size);

template <>
void init_cell_finder<parameters::exponential_spacing>(const double cell_growth_factor, const double first_cell_size);

template <>
void init_cell_finder<parameters::custom_spacing>(const double cell_growth_factor, const double first_cell_size);
*/


/*
template <parameters::t_radial_grid GRID>
int get_rmed_id(const double r);

template <parameters::t_radial_grid GRID>
int get_rinf_id(const double r);

template <parameters::t_radial_grid GRID>
int get_az_id(const double r);
*/

#endif // FIND_CELL_ID_H
