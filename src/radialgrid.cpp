#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <gsl/gsl_spline.h>
#include <mpi.h>
#include <sstream>

#include "LowTasks.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "mpi_utils.h"
#include "polargrid.h"
#include "radialarray.h"
#include "radialgrid.h"

t_radialgrid::t_radialgrid()
{
    m_size_radial = 0;
    m_name = NULL;
    m_unit = NULL;
    m_data = NULL;
    m_scalar = true;
    m_write_1D = false;
    m_calculate_on_write = false;
    m_do_before_write = NULL;
    m_clear_after_write = false;
}

t_radialgrid::~t_radialgrid()
{
    delete[] m_name;
    delete[] m_data;
}

/**
	set size of radial grid. grid is cleared automatically
*/
void t_radialgrid::set_size(ptrdiff_t size_radial)
{
    // delete old field
    delete[] m_data;

    m_size_radial = size_radial;

    // vector fields need one more cell in radial direction
    m_data = new double[(m_scalar ? size_radial : size_radial + 1)];

    clear();
}

/**
	set name of radial grid
*/
void t_radialgrid::set_name(const char *name)
{
    // delete old name
    delete[] m_name;

    // aquire space for new name
    m_name = new char[strlen(name) + 1];

    strcpy(m_name, name);
}

/**
	set unit of radial grid
*/
void t_radialgrid::set_unit(units::t_unit &unit) { m_unit = &unit; }

void t_radialgrid::set_vector(bool value) { m_scalar = !value; }

void t_radialgrid::set_scalar(bool value) { m_scalar = value; }

/**
	set all entries to 0
*/
void t_radialgrid::clear()
{
    std::memset(m_data, 0, sizeof(*m_data) * get_size_radial());
}

void t_radialgrid::write_radialgrid(unsigned int number, t_data &data)
{
    if (!get_write_1D()) {
	return;
    }

    if (m_do_before_write != NULL) {
	(*m_do_before_write)(data, number, false);
    }

    write1D(number);

    if (m_clear_after_write) {
	clear();
    }
}

/**
	 Write the data to binary file using the provided filename.

	 \param string filename : filename to write the data to
	 \param int number : output number
	 \param t_data data : data construct
	 \param bool one_file : write array without radii to one single file
*/
/* NO LONGER USED; writing repeatedly to single filed with mpi corrupted the
files after too many entries. void t_radialgrid::write(std::string filename,
unsigned int number, t_data &data, bool one_file, bool force_write)
{
    if (!force_write && !get_write_1D()) {
	return;
    }
    if (m_do_before_write != NULL) {
	(*m_do_before_write)(data, number, false);
    }

    write1D(filename, one_file);

    if (m_clear_after_write) {
	clear();
    }
}*/

/**
	write radial average values of polargrid to filename

	\param string filename
    \param bool one_file : write array without radii to one single file
*/
void t_radialgrid::write1D(unsigned int timestep) const
{

    /**
	    write radial average values of polargrid to filename

	    \param number file number
    */
    std::stringstream filename;
    filename << OUTPUTDIR << "/gas" << get_name() << timestep << ".dat";
    write1D(filename.str(), false);
}

void t_radialgrid::write1D(std::string filename, bool one_file) const
{
    MPI_File fh;
    MPI_Status status;
    int error;
    unsigned int count, from, number_of_values = 2;
    double *buffer;

    // use Rmed or Rinf depending if this quantity is scalar or vector
    t_radialarray &radius = is_scalar() ? Rmed : Rinf;

    // only write values and no radii for one_file
    if (one_file) {
	number_of_values = 1;
    }
    // try to open file
    if (one_file) {
	// append to existing file for consecutive output of arrays
	error = MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
			      MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL,
			      &fh);
	// if file doesn't exist yet, create it
	if (error != MPI_SUCCESS) {
	    error = MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
				  MPI_MODE_WRONLY | MPI_MODE_CREATE,
				  MPI_INFO_NULL, &fh);
	}

    } else {
	// make a new file for each output
	error = MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
			      MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,
			      &fh);
    }
    mpi_error_check_file_write(error, filename);

    // Move file pointer to correct position
    MPI_File_set_view(fh, 0, MPI_DOUBLE, MPI_DOUBLE,
		      const_cast<char *>("native"), MPI_INFO_NULL);
    if (one_file) {
	MPI_File_seek(fh, (IMIN + Zero_or_active) * number_of_values,
		      MPI_SEEK_END | MPI_SEEK_SET);
    } else {
	MPI_File_seek(fh, (IMIN + Zero_or_active) * number_of_values,
		      MPI_SEEK_SET);
    }
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

    if (one_file) {
	for (unsigned int n_radial = 0; n_radial < count; ++n_radial) {
	    buffer[number_of_values * n_radial] = operator()(from + n_radial);
	}
    } else {
	for (unsigned int n_radial = 0; n_radial < count; ++n_radial) {
	    buffer[number_of_values * n_radial] = radius[from + n_radial];
	    buffer[number_of_values * n_radial + 1] = operator()(from +
								 n_radial);
	}
    }

    // if unit is set multiply values by factor
    unsigned int value_offset = 1;
    if (one_file)
	value_offset = 0;
    if (m_unit) {
	for (unsigned int n_radial = 0; n_radial < count; ++n_radial) {
	    buffer[number_of_values * n_radial + value_offset] *=
		m_unit->get_cgs_factor();
	}
    }

    MPI_File_write_all(fh, buffer, count * number_of_values, MPI_DOUBLE,
		       &status);

    delete[] buffer;

    // close file
    MPI_File_close(&fh);
}

/**
	read radial average values of polargrid to filename

	\param number file number
*/
void t_radialgrid::read1D(unsigned int number)
{
    char *filename;

    if (asprintf(&filename, "%s/gas%s%i.dat", OUTPUTDIR, get_name(), number) <
	0) {
	die("Not enough memory!");
    }

    read1D(filename);

    free(filename);
}

void t_radialgrid::read1D(const char *_filename)
{
    MPI_File fh;
    MPI_Status status;
    MPI_Offset size;
    unsigned int count, number_of_values = 2;
    const std::string filename = std::string(_filename);

    // use Rmed or Rinf depending if this quantity is scalar or vector
    t_radialarray &radius = is_scalar() ? Rmed : Rinf;

    // try to open file
    mpi_error_check_file_read(MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
					    MPI_MODE_RDONLY, MPI_INFO_NULL,
					    &fh),
			      filename);

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
	operator()(n_radial) = gsl_spline_eval(spline, radius[n_radial], acc);
    }

    // free spline interpolation structures
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
}

/**
	calculate how much bytes are need for one 1D output
*/
unsigned int t_radialgrid::bytes_needed_1D()
{
    // factor 2 is because of radius and value
    unsigned int number_of_values = 2;

    return number_of_values * sizeof(double) * get_size_radial();
}

unsigned int t_radialgrid::get_memory_usage(ptrdiff_t size_radial)
{
    return (m_scalar ? size_radial : size_radial + 1) * sizeof(double);
}

t_radialgrid &t_radialgrid::operator=(const t_polargrid &polargrid)
{
    if ((polargrid.get_size_radial() == get_size_radial()) &&
	(polargrid.is_scalar() == is_scalar())) {
	for (unsigned int n_radial = 0; n_radial <= get_max_radial();
	     ++n_radial) {
	    operator()(n_radial) = 0.0;
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= polargrid.get_max_azimuthal(); ++n_azimuthal) {
		operator()(n_radial) += polargrid(n_radial, n_azimuthal);
	    }
	    operator()(n_radial) /= (double)polargrid.get_size_azimuthal();
	}

    } else {
	if (polargrid.get_size_radial() != get_size_radial()) {
	    die("Size of '%s' (%u) and '%s' (%u) does not match!\n",
		polargrid.get_name(), polargrid.get_size_radial(), get_name(),
		get_size_radial());
	} else {
	    die("Type of '%s' and '%s' does not match!\n", polargrid.get_name(),
		get_name());
	}
    }

    return *this;
}
