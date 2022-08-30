#include "random.h"

#include <random>

namespace fargo_random
{

cxx::ziggurat_normal_distribution<double> std_normal_dist{};
std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

std::unique_ptr<jsf64> gen_jsf;
uint64_t a_seed[4];

void init()
{
    seed();
    gen_jsf.reset(new jsf64(a_seed[0]));
}

void seed()
{
    uint32_t seeds[8];
    std::array<uint32_t, 11> entropy = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    randutils::auto_seed_256 seeder(entropy);

    seeder.generate(seeds, seeds + 8);

    a_seed[0] = (uint64_t)seeds[0] << 32 | seeds[1];
    a_seed[1] = (uint64_t)seeds[2] << 32 | seeds[3];
    a_seed[2] = (uint64_t)seeds[4] << 32 | seeds[5];
    a_seed[3] = (uint64_t)seeds[6] << 32 | seeds[6];
}

double std_normal() { return std_normal_dist(*gen_jsf); }
double uniform() { return uniform_dist(*gen_jsf); }

} // namespace fargo_random