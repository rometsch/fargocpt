#ifndef _STRESS_H_
#define _STRESS_H_

#include "data.h"

namespace stress {

void calculate_Reynolds_stress(t_data &data);
void calculate_gravitational_stress(t_data &data);

}

#endif
