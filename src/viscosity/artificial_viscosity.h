#pragma once

#include "../data.h"

namespace art_visc{

void update_with_artificial_viscosity(t_data &data, const double dt);
void update_with_artificial_viscosity_TW(t_data &data, const double dt);
void update_with_artificial_viscosity_SN(t_data &data, const double dt);

}
