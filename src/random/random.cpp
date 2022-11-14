#include "random.h"

#ifdef _OPENMP
#include <omp.h>
#else
#pragma message "Your compiler does not support OpenMP, at least with the flags you're using."
#endif
#include <random>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>


namespace fargo_random
{

std::vector<cxx::ziggurat_normal_distribution<double>> std_normal_dists;
std::vector<std::uniform_real_distribution<double>> uniform_dists;

std::vector<jsf64> gen_jsfs;
uint64_t a_seed[4];

static void dump_numbers() {
    // std::vector<std::vector<double>> nums;
    const unsigned int N = 10000;
    // #ifdef _OPENMP
    //     const unsigned int Nthreads = omp_get_num_threads();
    // #else
    //     const unsigned int Nthreads = 1;
    // #endif
    // nums.resize(Nthreads);
    // for (unsigned int n=0; n<Nthreads; n++) {
    //     nums[n].resize(N);
    // }

    #pragma omp parallel for
    for (unsigned int n=0; n<6; n++) {
        #ifdef _OPENMP
            const unsigned int omp_id  = omp_get_thread_num();
        #else
            const unsigned int omp_id = 0;
        #endif

        const std::string filename = "rand/" + std::to_string(omp_id);
        std::ofstream out(filename);

        for (unsigned int i=0; i<N; i++) {
            out << std_normal() << "\t" << uniform() << std::endl;
        }

        out.close();
    }
    // std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!! Writing numbers !!!!!!!!!!!!" << std::endl;
    // std::cout << "Nthreads = " << Nthreads << " , N = " << N << std::endl;
    // for (unsigned int n=0; n<Nthreads; n++) {
    //     for (unsigned int i=0; i<N; i++) {
    //         std::cout << n << "\t" << nums[n][i] << std::endl;
    //     }
    // }

}

void init() {

    // #ifdef _OPENMP
    //     unsigned int omp_num = omp_get_num_threads();
    // #else
    //     unsigned int omp_num = 1;
    // #endif
    int omp_num;
    #ifdef _OPENMP
    #pragma omp parallel
    {
        // int omp_id  = omp_get_thread_num();
        omp_num = omp_get_num_threads();
    }
    #else
    // int omp_id  = 0;
    int omp_num = 1;
    #endif

    std::cout << "number of threads = " << omp_num << std::endl;

    for (int n=0; n<omp_num; n++) {
        seed(n);
        gen_jsfs.emplace_back(jsf64(a_seed[n]));
        uniform_dists.emplace_back(std::uniform_real_distribution<double>(0.0, 1.0));
        std_normal_dists.emplace_back(cxx::ziggurat_normal_distribution<double>());
    }

    dump_numbers();


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

double std_normal() { 
    #ifdef _OPENMP
        const unsigned int omp_id  = omp_get_thread_num();
    #else
        const unsigned int omp_id = 0;
    #endif
    return std_normal_dists[omp_id](gen_jsfs[omp_id]); 
}
double uniform() {
    #ifdef _OPENMP
        const unsigned int omp_id  = omp_get_thread_num();
    #else
        const unsigned int omp_id = 0;
    #endif
    return uniform_dists[omp_id](gen_jsfs[omp_id]); 
}

} // namespace fargo_random