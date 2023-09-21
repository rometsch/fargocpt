#pragma once

#include "../data.h"
#include "../types.h"

namespace boundary_conditions
{

/****************************************
/// Basic
****************************************/

void initial_boundary_inner(t_data &data);
void initial_boundary_outer(t_data &data);

void open_boundary_inner(t_data &data);
void open_boundary_outer(t_data &data);

void zero_gradient_boundary_inner(t_data &data);
void zero_gradient_boundary_outer(t_data &data);

void reflecting_boundary_inner(t_data &data);
void reflecting_boundary_outer(t_data &data);

void keplerian2d_boundary_inner(t_data &data);
void keplerian2d_boundary_outer(t_data &data);

void viscous_outflow_boundary_inner(t_data &data);
} // namespace boundary_conditions