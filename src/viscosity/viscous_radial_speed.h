#ifndef VISCOUS_RADIAL_SPEED_H
#define VISCOUS_RADIAL_SPEED_H

#include "../types.h"
#include "../polargrid.h"
#include "../radialarray.h"

namespace viscous_speed
{

double get_vr_with_numerical_viscous_speed(const double r, const double mass);

} // viscous_speed

#endif // VISCOUS_RADIAL_SPEED_H
