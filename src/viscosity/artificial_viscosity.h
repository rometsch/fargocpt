#ifndef ARTIFICIAL_VISCOSITY_H
#define ARTIFICIAL_VISCOSITY_H

#include "../data.h"

namespace art_visc{

void update_with_artificial_viscosity(t_data &data, const double time, const double dt);
void update_with_artificial_viscosity_TW(t_data &data, const double dt);
void update_with_artificial_viscosity_TW_old(t_data &data, const double dt);
void update_with_artificial_viscosity_SN(t_data &data, const double dt);

}

#endif // ARTIFICIAL_VISCOSITY_H
