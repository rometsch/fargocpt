#pragma once

#include <string>

bool is_number(std::string s);

unsigned int get_next_azimuthal_id(const unsigned int id);
unsigned int get_prev_azimuthal_id(const unsigned int id);
unsigned int get_cell_id(const int nRadial, const int nAzimuthal);

void sum_without_ghost_cells(double &accumulator, const double &addend,
			     const unsigned int &n_radial);

double cutoff_outer(double point, double width, double x);
double cutoff_inner(double point, double width, double x);
bool is_big_endian(void);
bool is_distance_zero(double x);
