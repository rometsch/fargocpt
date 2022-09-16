#ifndef VISCOSITY_H
#define VISCOSITY_H

#include "../data.h"

namespace viscosity
{

void update_viscosity(t_data &data);
void compute_viscous_terms(t_data &data, bool include_artifical_viscosity);
void update_velocities_with_viscosity(t_data &data, const double dt);

} // namespace viscosity

#endif // VISCOSITY_H