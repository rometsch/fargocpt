#ifndef GAS_TORQUES_H
#define GAS_TORQUES_H

#include "data.h"

namespace gas_torques
{

void calculate_advection_torque(t_data &data, t_polargrid& t_adv, const double dt);
void calculate_viscous_torque(t_data &data, t_polargrid& t_visc, const double dt);
void calculate_gravitational_torque(t_data &data,  t_polargrid &t_grav, const double dt);

} // namespace gas_torques
#endif // GAS_TORQUES_H
