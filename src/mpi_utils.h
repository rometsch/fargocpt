#ifndef MPI_UTILS_H
#define MPI_UTILS_H
#include <string>

void mpi_error_check(const int error, const std::string msg);
void mpi_error_check(const int error);

void mpi_error_check_file_write(const int error, const std::string filename);
void mpi_error_check_file_read(const int error, const std::string filename);

#endif // MPI_UTILS_H
