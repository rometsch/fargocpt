#ifndef SOURCEEULER_H
#define SOURCEEULER_H

#include "data.h"
#include "nbody/planet.h"
#include "types.h"

bool assure_minimum_value(t_polargrid &dst, double minimum_value);
bool assure_temperature_range(t_data &data);

void move_polargrid(t_polargrid &dst, t_polargrid &src);
void copy_polargrid(t_polargrid &dst, const t_polargrid& src);

void recalculate_derived_disk_quantities(t_data &data, const double current_time);
void recalculate_viscosity(t_data &data, const double current_time);
void init_euler(t_data &data, const double current_time);
void FreeEuler();

void update_with_sourceterms(t_data &data, const double dt);
void SubStep3(t_data &data, const double current_time, const double dt);
void radiative_diffusion(t_data &data, const double current_time, const double dt);
void compression_heating(t_data &data, const double dt);

void calculate_qplus(t_data &data);
void calculate_qminus(t_data &data, const double current_time);


void compute_sound_speed(t_data &data, const double current_time);
void compute_scale_height(t_data &data, const double current_time);
void compute_scale_height_old(t_data &data);
void compute_scale_height_nbody(t_data &data, const double current_time);
void compute_scale_height_center_of_mass(t_data &data);
void compute_pressure(t_data &data);
void compute_temperature(t_data &data);

void SetTemperatureFloorCeilValues(t_data &data, std::string filename,
				   int line);

void compute_heating_cooling_for_CFL(t_data &data, const double current_time);

#endif // SOURCEEULER_H
