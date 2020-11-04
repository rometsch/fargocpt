#ifndef MPI_UTILS_H
#define MPI_UTILS_H
#include <string>

void mpi_error_check(const int error, const std::string msg);
void mpi_error_check(const int error);

#endif // MPI_UTILS_H
