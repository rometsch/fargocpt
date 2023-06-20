#pragma once

#include "jsf.hpp"
#include "randutils.hpp"
#include "ziggurat.hpp"

namespace fargo_random
{

void init();
void seed(unsigned int);

double std_normal();
double uniform();

} // namespace fargo_random
