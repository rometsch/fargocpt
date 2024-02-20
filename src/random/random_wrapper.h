#pragma once

#include "jsf.hpp"
#include "randutils.hpp"
#include "ziggurat.hpp"

namespace fargo_random
{

void init();
void seed(unsigned int);

double get_std_normal();
double get_uniform_one();
double get_uniform_two_pi();
double get_uniform_dust_eccentricity();
int get_uniform256();

} // namespace fargo_random
