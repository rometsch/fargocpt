#ifndef VISCOUS_RADIAL_SPEED_H
#define VISCOUS_RADIAL_SPEED_H

#include "../types.h"
#include "../polargrid.h"
#include "../radialarray.h"

namespace viscous_speed
{

double get_vr_with_numerical_viscous_speed(const double r, const double mass);
double get_vr_with_numerical_viscous_speed_wrapper(const t_radialarray &Rcenter, const t_radialarray &Rinterface, const Pair center_pos, const double mass, const int nr, const int np);

} // viscous_speed

#endif // VISCOUS_RADIAL_SPEED_H
