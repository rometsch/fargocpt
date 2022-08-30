#ifndef RANDOM_H
#define RANDOM_H

#include "jsf.hpp"
#include "randutils.hpp"
#include "ziggurat.hpp"

namespace fargo_random
{

void init();
void seed();

double std_normal();
double uniform();

} // namespace fargo_random

#endif // RANDOM_H
