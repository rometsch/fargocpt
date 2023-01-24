#pragma once

#ifdef _OPENMP
#include <omp.h>
#else
#pragma message "Your compiler does not support OpenMP, at least with the flags you're using."
#endif

#include <mpi.h>

void init_parallel(int argc, char *argv[]);
void finalize_parallel();