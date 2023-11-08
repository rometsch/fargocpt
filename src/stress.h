#pragma once
#include "data.h"

namespace stress
{

void calculate_Reynolds_stress(t_data &data);
void calculate_gravitational_stress(t_data &data);

} // namespace stress
