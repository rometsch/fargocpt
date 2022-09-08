#ifndef CFL_H
#define CFL_H

#include "data.h"
#include "types.h"

namespace cfl {

double condition_cfl(t_data &data, t_polargrid &v_radial,
		     t_polargrid &v_azimuthal, t_polargrid &soundspeed,
		     const double deltaT);

}

#endif