#pragma once

#include "../types.h"
#include "../polargrid.h"
#include "../radialarray.h"
#include "../data.h"

namespace viscous_speed
{

void init_vr_table_boundary(t_data &data);
double lookup_initial_vr_outer(const double r);
double lookup_initial_vr_inner(const double r);
double get_vr_with_numerical_viscous_speed(const double r, const double mass);
double get_vr_outer_viscous_speed_correction_factor(const double r, const double mass);
} // viscous_speed
