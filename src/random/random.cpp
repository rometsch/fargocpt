#include "random.h"
#include "../logging.h"

#ifdef _OPENMP
#include <omp.h>
#else
#pragma message "Your compiler does not support OpenMP, at least with the flags you're using."
#endif
#include <random>
#include <vector>



namespace fargo_random
{

std::vector<cxx::ziggurat_normal_distribution<double>> std_normal_dists;
std::vector<std::uniform_real_distribution<double>> uniform_dists;

std::vector<jsf64> gen_jsfs;
uint64_t a_seed[4];

static unsigned int get_number_required_rngs() {
    /* Return the number of required rngs per mpi process.
    Each openmp threads needs its own rngs because the used ones are not thread safe.
    Thus, the return value of this function corresponds to the number of threads used.
    */
    #ifdef _OPENMP
    unsigned int N = omp_get_max_threads();
    #else
    unsigned int N = 1;
    #endif
    return N;
}

void init() {

    const unsigned int Nrngs = get_number_required_rngs();
    logging::print_master(LOG_INFO "Initializing %d RNGs per MPI process.", Nrngs);

    for (unsigned int n=0; n<Nrngs; n++) {
        seed(n);
        gen_jsfs.emplace_back(jsf64(a_seed[n]));
        uniform_dists.emplace_back(std::uniform_real_distribution<double>(0.0, 1.0));
        std_normal_dists.emplace_back(cxx::ziggurat_normal_distribution<double>());
    }

}

void seed(unsigned int n)
{
    uint32_t seeds[8];
    std::array<uint32_t, 11> entropy = {1*(n+1), 2*(n+1), 3*(n+1), 4*(n+1), 5*(n+1), 6*(n+1), 7*(n+1), 8*(n+1), 9*(n+1), 10*(n+1), 11*(n+1)};
    randutils::auto_seed_256 seeder(entropy);

    seeder.generate(seeds, seeds + 8);

    a_seed[0] = (uint64_t)seeds[0] << 32 | seeds[1];
    a_seed[1] = (uint64_t)seeds[2] << 32 | seeds[3];
    a_seed[2] = (uint64_t)seeds[4] << 32 | seeds[5];
    a_seed[3] = (uint64_t)seeds[6] << 32 | seeds[6];
}

static unsigned int thread_num() {
    #ifdef _OPENMP
        const unsigned int n  = omp_get_thread_num();
    #else
        const unsigned int n = 0;
    #endif
    return n;
}

double std_normal() { 
    const unsigned int n = thread_num();
    return std_normal_dists[n](gen_jsfs[n]); 
}
double uniform() {
    const unsigned int n = thread_num();
    return uniform_dists[n](gen_jsfs[n]); 
}

} // namespace fargo_random