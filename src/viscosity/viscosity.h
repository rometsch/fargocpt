#ifndef VISCOSITY_H
#define VISCOSITY_H

#include "../data.h"

namespace viscosity
{

double get_alpha(int nr, int naz, t_data &data);
void update_viscosity(t_data &data);
void compute_viscous_stress_tensor(t_data &data);
void update_velocities_with_viscosity(t_data &data, const double dt);

} // namespace viscosity

#endif // VISCOSITY_H
