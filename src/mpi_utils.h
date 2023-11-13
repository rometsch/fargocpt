#pragma once

#include <string>

void mpi_error_check(const int error, const std::string msg);
void mpi_error_check(const int error);

void mpi_error_check_file_write(const int error, const std::string filename);
void mpi_error_check_file_read(const int error, const std::string filename);
