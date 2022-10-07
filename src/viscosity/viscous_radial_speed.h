#ifndef VISCOUS_RADIAL_SPEED_H
#define VISCOUS_RADIAL_SPEED_H

#include "../types.h"
#include "../polargrid.h"
#include "../radialarray.h"
#include "../data.h"

namespace viscous_speed
{

void init_vr_table_outer_boundary(t_data &data);
double lookup_initial_vr(const double r);
double get_vr_with_numerical_viscous_speed(const double r, const double mass);

} // viscous_speed

#endif // VISCOUS_RADIAL_SPEED_H
