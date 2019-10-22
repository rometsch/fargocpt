/**
	\file polargrid.cpp
	\author Tobias Mueller <Tobias_Mueller@twam.info>
*/

#include "polargrid.h"
#include "LowTasks.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "util.h"
#include <float.h>
#include <gsl/gsl_spline.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

t_polargrid::t_polargrid()
{
    Nrad = 0;
    Nsec = 0;
    m_name = NULL;
    m_unit = NULL;
    Field = NULL;
    m_scalar = true;
    m_write_1D = false;
    m_write_2D = false;
    m_calculate_on_write = false;
    m_write_max_max_1D = true;
    m_do_before_write = NULL;
}

t_polargrid::~t_polargrid()
{
    delete[] m_name;
    delete[] Field;
}

/**
	set size of polargrid. grid is cleared automatically
*/
void t_polargrid::set_size(ptrdiff_t size_radial, ptrdiff_t size_azimuthal)
{
    // delete old field
    delete[] Field;

    Nrad = size_radial;
    Nsec = size_azimuthal;

    // vector fields need one more cell in radial direction
    Field = new double[(m_scalar ? size_radial : size_radial + 1) *
		       (size_azimuthal)];

    clear();
}

/**
	set name of polargrid
*/
void t_polargrid::set_name(const char *name)
{
    // delete old name
    delete[] m_name;

    // aquire space for new name
    m_name = new char[strlen(name) + 1];

    strcpy(m_name, name);
}

/**
	set unit of polargrid
*/
void t_polargrid::set_unit(units::t_unit &unit) { m_unit = &unit; }

void t_polargrid::set_vector(bool value) { m_scalar = !value; }

void t_polargrid::set_scalar(bool value) { m_scalar = value; }

/**
	set all entries to 0
*/
void t_polargrid::clear()
{
    for (unsigned int n_radial = 0; n_radial <= get_max_radial(); ++n_radial) {
	for (unsigned int n_azimuthal = 0; n_azimuthal <= get_max_azimuthal();
	     ++n_azimuthal) {
	    operator()(n_radial, n_azimuthal) = 0.0;
	}
    }
}

void t_polargrid::write(unsigned int number, t_data &data) const
{
    if (get_write_1D() || get_write_2D() || m_calculate_on_write) {
	if (m_do_before_write != NULL) {
	    (*m_do_before_write)(data, number, false);
	}
    }

    if (get_write_1D())
	write1D(number);

    if (get_write_2D())
	write2D(number);
}

/**
	write polargrid to file

	\param number file number
*/
void t_polargrid::write2D(unsigned int number) const
{
    MPI_File fh;
    MPI_Status status;
    int error, error_class, error_length;
    unsigned int count;
    double *from;
    char *filename, error_string[MPI_MAX_ERROR_STRING + 1];

    if (asprintf(&filename, "%s/gas%s%i.dat", OUTPUTDIR, get_name(), number) <
	0) {
	die("Not enough memory!");
    }

    // try to open file
    error =
	MPI_File_open(MPI_COMM_WORLD, filename,
		      MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    if (error != MPI_SUCCESS) {
	logging::print_master(
	    LOG_ERROR
	    "Error while writing to file '%s'. Check file permissions and IO support of MPI library\n",
	    filename);

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
    free(filename);

    MPI_File_set_view(fh, 0, MPI_DOUBLE, MPI_DOUBLE,
		      const_cast<char *>("native"), MPI_INFO_NULL);
    MPI_File_seek(fh, (IMIN + Zero_or_active) * get_size_azimuthal(),
		  MPI_SEEK_SET);

    from = Field;
    count = get_size_radial();

    // only the hightest CPU should print the last cell of vector grids
    if (is_vector() && (CPU_Rank != CPU_Highest)) {
	count -= 1;
    }

    // We strip the first CPUOVERLAP rings if the current CPU is not the
    // innermost one
    if (CPU_Rank > 0) {
	from += CPUOVERLAP * get_size_azimuthal();
	count -= CPUOVERLAP;
    }

    // We strip the last CPUOVERLAP rings if the current CPU is not the
    // outermost one, equal to CPU_Highest in all cases
    if (CPU_Rank != CPU_Highest) {
	count -= CPUOVERLAP;
    }

    // if unit is set, copy values to temporary buffer and multiply them by
    // factor
    if (m_unit) {
	double *buffer;

	// copy data to temporary buffer
	buffer = new double[count * get_size_azimuthal()];

	// add unit factor
	memcpy(buffer, from, (count) * (get_size_azimuthal()) * sizeof(double));
	for (unsigned int i = 0; i < (count * get_size_azimuthal()); ++i) {
	    buffer[i] *= m_unit->get_cgs_factor();
	}

	// write data from temporary buffer
	MPI_File_write(fh, buffer, (count) * (get_size_azimuthal()), MPI_DOUBLE,
		       &status);

	// delete temporary buffer
	delete[] buffer;
    } else {
	// write data from buffer
	MPI_File_write(fh, from, (count) * (get_size_azimuthal()), MPI_DOUBLE,
		       &status);
    }

    // close file
    MPI_File_close(&fh);
}

/**
	write radial average values of polargrid to filename

	\param number file number
*/
void t_polargrid::write1D(unsigned int number) const
{
    MPI_File fh;
    MPI_Status status;
    int error, error_class, error_length;
    unsigned int count, from, number_of_values = 2;
    double *buffer;
    char *filename, error_string[MPI_MAX_ERROR_STRING + 1];

    // use Rmed or Rinf depending if this quantity is scalar or vector
    t_radialarray &radius = is_scalar() ? Rb : Ra;

    if (asprintf(&filename, "%s/gas%s1D%i.dat", OUTPUTDIR, get_name(), number) <
	0) {
	die("Not enough memory!");
    }

    // try to open file
    error =
	MPI_File_open(MPI_COMM_WORLD, filename,
		      MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    if (error != MPI_SUCCESS) {
	logging::print_master(
	    LOG_ERROR
	    "Error while writing to file '%s'. Check file permissions and IO support of MPI library\n",
	    filename);

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
    free(filename);

    if (m_write_max_max_1D) {
	// min/max need additional two values
	number_of_values += 2;
    }

    MPI_File_set_view(fh, 0, MPI_DOUBLE, MPI_DOUBLE,
		      const_cast<char *>("native"), MPI_INFO_NULL);
    MPI_File_seek(fh, (IMIN + Zero_or_active) * number_of_values, MPI_SEEK_SET);

    from = 0;
    count = get_size_radial();

    // only the hightest CPU should print the last cell of vector grids
    if (is_vector() && (CPU_Rank != CPU_Highest)) {
	count -= 1;
    }

    // We strip the first CPUOVERLAP rings if the current CPU is not the
    // innermost one
    if (CPU_Rank > 0) {
	from += CPUOVERLAP;
	count -= CPUOVERLAP;
    }

    // We strip the last CPUOVERLAP rings if the current CPU is not the
    // outermost one, equal to CPU_Highest in all cases
    if (CPU_Rank != CPU_Highest) {
	count -= CPUOVERLAP;
    }

    buffer = new double[number_of_values * count];

    for (unsigned int n_radial = 0; n_radial < count; ++n_radial) {
	buffer[number_of_values * n_radial] = radius[from + n_radial];
	buffer[number_of_values * n_radial + 1] = 0;

	for (unsigned int n_azimuthal = 0; n_azimuthal <= get_max_azimuthal();
	     ++n_azimuthal) {
	    buffer[number_of_values * n_radial + 1] += operator()(
		from + n_radial, n_azimuthal);
	}

	buffer[number_of_values * n_radial + 1] /= (double)get_size_azimuthal();
    }

    // if write min/max values, get them!
    if (m_write_max_max_1D) {
	for (unsigned int n_radial = 0; n_radial < count; ++n_radial) {
	    // min
	    buffer[number_of_values * n_radial + 2] = DBL_MAX;
	    // max
	    buffer[number_of_values * n_radial + 3] = -DBL_MAX;

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= get_max_azimuthal(); ++n_azimuthal) {
		buffer[number_of_values * n_radial + 2] =
		    min(buffer[number_of_values * n_radial + 2],
			operator()(from + n_radial, n_azimuthal));
		buffer[number_of_values * n_radial + 3] =
		    max(buffer[number_of_values * n_radial + 3],
			operator()(from + n_radial, n_azimuthal));
	    }
	}
    }

    // if unit is set multiply values by factor
    if (m_unit) {
	for (unsigned int n_radial = 0; n_radial < count; ++n_radial) {
	    buffer[number_of_values * n_radial + 1] *= m_unit->get_cgs_factor();

	    if (m_write_max_max_1D) {
		buffer[number_of_values * n_radial + 2] *=
		    m_unit->get_cgs_factor();
		buffer[number_of_values * n_radial + 3] *=
		    m_unit->get_cgs_factor();
	    }
	}
    }

    MPI_File_write_all(fh, buffer, count * number_of_values, MPI_DOUBLE,
		       &status);

    delete[] buffer;

    // close file
    MPI_File_close(&fh);
}

void t_polargrid::read2D(unsigned int number)
{
    char *filename;

    if (asprintf(&filename, "%s/gas%s%i.dat", OUTPUTDIR, get_name(), number) <
	0) {
	die("Not enough memory!");
    }

    read2D(filename);

    free(filename);
}

void t_polargrid::read2D(const char *_filename)
{
    MPI_File fh;
    MPI_Status status;
    MPI_Offset size;
    int error, error_class, error_length;
    unsigned int count;
    char *filename, error_string[MPI_MAX_ERROR_STRING + 1];

    filename = new char[strlen(_filename) + 1];
    strcpy(filename, _filename);

    error = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,
			  MPI_INFO_NULL, &fh);
    if (error != MPI_SUCCESS) {
	logging::print_master(
	    LOG_ERROR
	    "Error while reading from file '%s'. Check file permissions and IO support of MPI library\n",
	    filename);

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
    delete[] filename;

    // get file size
    MPI_File_get_size(fh, &size);

    count = (GlobalNRadial + (is_scalar() ? 0 : 1)) * Nsec * sizeof(double);

    if (count != size) {
	die("Filename '%s' has %u bytes but has to be %u bytes.", _filename,
	    size, count);
    }

    logging::print_master(LOG_INFO "Reading file '%s' with %u bytes.\n",
			  _filename, size);

    // allocate buffer and read file
    double *buffer_file = new double[count];

    MPI_File_set_view(fh, 0, MPI_DOUBLE, MPI_DOUBLE,
		      const_cast<char *>("native"), MPI_INFO_NULL);
    MPI_File_seek(fh, 0, MPI_SEEK_SET);
    MPI_File_read_all(fh, buffer_file, count, MPI_DOUBLE, &status);

    // close file
    MPI_File_close(&fh);

    // if unit is set multiply values by factor
    if (m_unit) {
	for (unsigned int n_radial = 0; n_radial < count; ++n_radial) {
	    buffer_file[n_radial] /= m_unit->get_cgs_factor();
	}
    }

    // copy data into polargrid
    for (unsigned int n_radial = 0; n_radial <= get_max_radial(); ++n_radial) {
	for (unsigned int n_azimuthal = 0; n_azimuthal < Nsec; ++n_azimuthal) {
	    operator()(n_radial, n_azimuthal) =
		buffer_file[(IMIN + n_radial) * Nsec + n_azimuthal];
	}
    }

    // delete file buffer
    delete[] buffer_file;
}

/**
	read radial average values of polargrid to filename

	\param number file number
*/
void t_polargrid::read1D(unsigned int number, bool skip_min_max)
{
    char *filename;

    if (asprintf(&filename, "%s/gas%s1D%i.dat", OUTPUTDIR, get_name(), number) <
	0) {
	die("Not enough memory!");
    }

    read1D(filename, skip_min_max);

    free(filename);
}

void t_polargrid::read1D(const char *_filename, bool skip_min_max)
{
    MPI_File fh;
    MPI_Status status;
    MPI_Offset size;
    int error, error_class, error_length;
    unsigned int count, number_of_values = 2;
    char *filename, error_string[MPI_MAX_ERROR_STRING + 1];

    // use Rmed or Rinf depending if this quantity is scalar or vector
    t_radialarray &radius = is_scalar() ? Rmed : Rinf;

    filename = new char[strlen(_filename) + 1];
    strcpy(filename, _filename);

    // try to open file
    error = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,
			  MPI_INFO_NULL, &fh);
    if (error != MPI_SUCCESS) {
	logging::print_master(
	    LOG_ERROR
	    "Error while reading from file '%s'. Check file permissions and IO support of MPI library\n",
	    filename);

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
    delete[] filename;

    if (skip_min_max) {
	// min/max need additional two values
	number_of_values += 2;
    }

    // get file size
    MPI_File_get_size(fh, &size);

    // calculate number of values in file
    count = size / number_of_values / sizeof(double);
    if (is_vector()) {
	count -= 1;
    }

    logging::print_master(
	LOG_INFO "Reading file '%s' with %u bytes. Reading %u values...\n",
	_filename, size, count);

    // allocate buffer and read file
    double *buffer_file = new double[number_of_values * count];

    MPI_File_set_view(fh, 0, MPI_DOUBLE, MPI_DOUBLE,
		      const_cast<char *>("native"), MPI_INFO_NULL);
    MPI_File_seek(fh, 0, MPI_SEEK_SET);
    MPI_File_read_all(fh, buffer_file, count * number_of_values, MPI_DOUBLE,
		      &status);

    // close file
    MPI_File_close(&fh);

    // if unit is set multiply values by factor
    if (m_unit) {
	for (unsigned int n_radial = 0; n_radial < count; ++n_radial) {
	    buffer_file[number_of_values * n_radial + 1] /=
		m_unit->get_cgs_factor();
	}
    }

    // allocate buffers for radius & value and copy values from buffer to it
    double *buffer_radius = new double[count];
    double *buffer_value = new double[count];

    for (unsigned int n_radial = 0; n_radial < count; ++n_radial) {
	buffer_radius[n_radial] = buffer_file[number_of_values * n_radial + 0];
	buffer_value[n_radial] = buffer_file[number_of_values * n_radial + 1];
    }

    // print warning if using spline values outside from provide range
    if ((RMIN < 2 * buffer_radius[0] - buffer_radius[1]) ||
	(RMAX > 2 * buffer_radius[count - 1] - buffer_radius[count - 2])) {
	logging::print_master(
	    LOG_WARNING
	    "Warning: '%s' covers radii from %.5lf to %.5lf, but you're using %.5lf to %.5lf! Spline interpolation isn't performing well here!\n",
	    _filename, buffer_radius[0], buffer_radius[count - 1], RMIN, RMAX);
    }

    // delete file buffer
    delete[] buffer_file;

    // allocate structures needed for spline interpolation
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, count);

    // init spline
    gsl_spline_init(spline, buffer_radius, buffer_value, count);

    // free temporay buffers for radius/value
    delete[] buffer_radius;
    delete[] buffer_value;

    // evaluate spline on all points needed and write into polargrid
    for (unsigned int n_radial = 0; n_radial <= get_max_radial(); ++n_radial) {
	for (unsigned int n_azimuthal = 0; n_azimuthal <= get_max_azimuthal();
	     ++n_azimuthal) {
	    operator()(n_radial, n_azimuthal) =
		gsl_spline_eval(spline, radius[n_radial], acc);
	}
    }

    // free spline interpolation structures
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
}

/**
	calculate how much bytes are need for one 1D output
*/
unsigned int t_polargrid::bytes_needed_1D()
{
    // factor 2 is because of radius and value
    unsigned int number_of_values = 2;

    if (m_write_max_max_1D)
	number_of_values += 2;

    return number_of_values * sizeof(double) * get_size_radial();
}

/**
	calculate how much bytes are needed for one 2D output
*/
unsigned int t_polargrid::bytes_needed_2D()
{
    // factor 2 is because of radius and value
    return 2.0 * sizeof(double) * get_size_azimuthal() * get_size_radial();
}

/**
	multiply each polargrid entry with a constant
*/
t_polargrid &t_polargrid::operator*=(double c)
{
    // TODO: Speed up by using only one loop and so skipping calculation of
    // cells. This should be done if backward-compability bug (size =
    // (Nrad+1)*(Nsec+1) ) is fixed.

    for (unsigned int n_radial = 0; n_radial < Nrad; ++n_radial) {
	for (unsigned int n_azimuthal = 0; n_azimuthal < Nsec; ++n_azimuthal) {
	    operator()(n_radial, n_azimuthal) *= c;
	}
    }

    return *this;
}

/**
	divide each polargrid entry with a constant
*/
t_polargrid &t_polargrid::operator/=(double c)
{
    // TODO: Speed up by using only one loop and so skipping calculation of
    // cells. This should be done if backward-compability bug (size =
    // (Nrad+1)*(Nsec+1) ) is fixed.

    for (unsigned int n_radial = 0; n_radial < Nrad; ++n_radial) {
	for (unsigned int n_azimuthal = 0; n_azimuthal < Nsec; ++n_azimuthal) {
	    operator()(n_radial, n_azimuthal) /= c;
	}
    }

    return *this;
}

unsigned int t_polargrid::get_memory_usage(ptrdiff_t size_radial,
					   ptrdiff_t size_azimuthal)
{
    return (m_scalar ? size_radial : size_radial + 1) * (size_azimuthal) *
	   sizeof(double);
}
