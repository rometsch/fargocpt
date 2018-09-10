#ifndef _SOURCEEULER_H_
#define _SOURCEEULER_H_

#include "data.h"
#include "types.h"

int DetectCrash(t_polargrid* array);
bool assure_minimum_value(t_polargrid &dst, double minimum_value);
bool assure_minimum_temperature(t_polargrid &energy, t_polargrid &density, double minimum_value);
bool assure_maximum_temperature(t_polargrid &energy, t_polargrid &density, double maximum_value);

void copy_polargrid(t_polargrid &dst, t_polargrid &src);
void SwitchPolarGrid(t_polargrid *dst, t_polargrid *src);

void init_euler(t_data &data);
void FreeEuler();

void AlgoGas(unsigned int nTimeStep, Force* force, t_data &data);

void SubStep1(t_data &data, double dt);
void SubStep2(t_data &data, double dt);
void SubStep3(t_data &data, double dt);
void radiative_diffusion(t_data &data, double dt);

void calculate_qplus(t_data &data, double dt);
void calculate_qminus(t_data &data, double dt);

double condition_cfl(t_data &data, t_polargrid &v_radial, t_polargrid &v_azimuthal, t_polargrid &soundspeed, double deltaT);

void compute_sound_speed(t_data &data, bool force_update);
void compute_aspect_ratio(t_data &data, bool force_update);
void compute_pressure(t_data &data, bool force_update);
void compute_temperature(t_data &data, bool force_update);
void compute_rho(t_data &data, bool force_update);

double CircumPlanetaryMass(t_data &data);

#endif
