#ifndef _BOUNDARY_CONDITIONS_H_
#define _BOUNDARY_CONDITIONS_H_

#include "data.h"
#include "types.h"

namespace boundary_conditions {

void apply_boundary_condition_temperature(t_data &data);
void apply_boundary_condition(t_data &data, double dt, bool final);
void open_boundary_inner(t_data &data);
void open_boundary_outer(t_data &data);
void reflecting_boundary_inner(t_data &data);
void reflecting_boundary_outer(t_data &data);
void viscous_outflow_boundary_inner(t_data &data);
void damping(t_data &data, double dt);
void mass_overflow(t_data &data);
void boundary_layer_inner_boundary(t_data &data);
void boundary_layer_outer_boundary(t_data &data);
void keplerian2d_boundary_inner(t_data &data);
void keplerian2d_boundary_outer(t_data &data);

}

#endif
