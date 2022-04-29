#ifndef SOURCEEULER_H
#define SOURCEEULER_H

#include "data.h"
#include "types.h"

bool assure_minimum_value(t_polargrid &dst, double minimum_value);
bool assure_temperature_range(t_data &data);

void move_polargrid(t_polargrid &dst, t_polargrid &src);
void copy_polargrid(t_polargrid &dst, t_polargrid &src);
void SwitchPolarGrid(t_polargrid *dst, t_polargrid *src);

void recalculate_derived_disk_quantities(t_data &data, bool force_update);
void init_euler(t_data &data);
void FreeEuler();

void AlgoGas(t_data &data);

void update_with_sourceterms(t_data &data, double dt);
void update_with_artificial_viscosity(t_data &data, double dt);
void SubStep3(t_data &data, double dt);
void radiative_diffusion(t_data &data, double dt);

void calculate_qplus(t_data &data);
void calculate_qminus(t_data &data);

double condition_cfl(t_data &data, t_polargrid &v_radial,
			 t_polargrid &v_azimuthal, t_polargrid &soundspeed,
			 const double deltaT);

void compute_sound_speed(t_data &data, bool force_update);
void compute_scale_height(t_data &data, const bool force_update);
void compute_scale_height_old(t_data &data, const bool force_update);
void compute_scale_height_nbody(t_data &data, const bool force_update);
void compute_scale_height_center_of_mass(t_data &data, const bool force_update);
void compute_pressure(t_data &data, bool force_update);
void compute_temperature(t_data &data, bool force_update);
void compute_rho(t_data &data, bool force_update);

void ComputeViscousStressTensor(t_data &data);
void ComputeCircumPlanetaryMasses(t_data &data);
void SetTemperatureFloorCeilValues(t_data &data, std::string filename,
				   int line);

void compute_heating_cooling_for_CFL(t_data &data);

#endif // SOURCEEULER_H
