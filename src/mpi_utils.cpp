#include "mpi_utils.h"
#include "logging.h"
#include <mpi.h>
#include <stdlib.h>

void mpi_error_check(const int error, const std::string msg)
{
    const std::string err = LOG_ERROR;
    const std::string err_msg = err + msg;

    if (error != MPI_SUCCESS) {
	logging::print_master(err_msg.c_str());
    }
}

void mpi_error_check(const int error)
{
    int error_class, error_length;
    char error_string[MPI_MAX_ERROR_STRING + 1];

    if (error != MPI_SUCCESS) {

	// error class
	MPI_Error_class(error, &error_class);
	MPI_Error_string(error_class, error_string, &error_length);
	error_string[error_length] = 0;
	logging::print_master(LOG_ERROR "MPI error class: %s\n", error_string);

	// error code
	MPI_Error_string(error, error_string, &error_length);
	error_string[error_length] = 0;
	logging::print_master(LOG_ERROR "MPI error code: %s\n", error_string);

	MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
}

void mpi_error_check_file_write(const int error, const std::string filename)
{
    const std::string err_msg =
	"Error while writing to file" + filename +
	". Check file permissions and IO support of MPI library\n";
    mpi_error_check(error, err_msg);
}

void mpi_error_check_file_read(const int error, const std::string filename)
{
    const std::string err_msg =
	"Error while reading to file" + filename +
	". Check file permissions and IO support of MPI library\n";
    mpi_error_check(error, err_msg);
}
