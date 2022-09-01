#include "output.h"
#include "Force.h"
#include "LowTasks.h"
#include "Pframeforce.h"
#include "SideEuler.h"
#include "Theo.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "nongnu.h"
#include "options.h"
#include "parameters.h"
#include "quantities.h"
#include "start_mode.h"
#include "stress.h"
#include "util.h"
#include "viscosity.h"

#include <dirent.h>

#include "unistd.h" // for access()
#include <cfloat>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <sys/stat.h>
#include <sys/statvfs.h>

namespace output
{

const std::map<const std::string, const int> quantities_file_column_v2 = {
    {"time step", 0},
    {"physical time", 1},
    {"mass", 2},
    {"angular momentum", 3},
    {"total energy", 4},
    {"internal energy", 5},
    {"kinematic energy", 6},
    {"potential energy", 7},
    {"radial kinetic energy", 8},
    {"azimuthal kinetic energy", 9},
    {"eccentricity", 10},
    {"periastron", 11},
    {"qplus", 12},
    {"qminus", 13},
    {"pdivv", 14},
    {"delta mass inner positive", 15},
    {"delta mass inner negative", 16},
    {"delta mass outer positive", 17},
    {"delta mass outer negative", 18},
    {"delta mass wave damping positive", 19},
    {"delta mass wave damping negative", 20},
    {"delta mass floor density positive", 21}};

const std::map<const std::string, const int> quantities_file_column_v2_1 = {
    {"time step", 0},
    {"analysis time step", 1},
    {"physical time", 2},
    {"mass", 3},
    {"angular momentum", 4},
    {"total energy", 5},
    {"internal energy", 6},
    {"kinematic energy", 7},
    {"potential energy", 8},
    {"radial kinetic energy", 9},
    {"azimuthal kinetic energy", 10},
    {"eccentricity", 11},
    {"periastron", 12},
    {"qplus", 13},
    {"qminus", 14},
    {"pdivv", 15},
    {"delta mass inner positive", 16},
    {"delta mass inner negative", 17},
    {"delta mass outer positive", 18},
    {"delta mass outer negative", 19},
    {"delta mass wave damping positive", 20},
    {"delta mass wave damping negative", 21},
    {"delta mass floor density positive", 22}};

const std::map<const std::string, const int> quantities_file_column_v2_2 = {
    {"time step", 0},
    {"analysis time step", 1},
    {"physical time", 2},
    {"mass", 3},
    {"radius", 4},
    {"angular momentum", 5},
    {"total energy", 6},
    {"internal energy", 7},
    {"kinematic energy", 8},
    {"potential energy", 9},
    {"radial kinetic energy", 10},
    {"azimuthal kinetic energy", 11},
    {"eccentricity", 12},
    {"periastron", 13},
    {"qplus", 14},
    {"qminus", 15},
    {"pdivv", 16},
    {"delta mass inner positive", 17},
    {"delta mass inner negative", 18},
    {"delta mass outer positive", 19},
    {"delta mass outer negative", 20},
    {"delta mass wave damping positive", 21},
    {"delta mass wave damping negative", 22},
    {"delta mass floor density positive", 23}};

const std::map<const std::string, const int> quantities_file_column_v2_3 = {
    {"time step", 0},
    {"analysis time step", 1},
    {"physical time", 2},
    {"mass", 3},
    {"radius", 4},
    {"angular momentum", 5},
    {"total energy", 6},
    {"internal energy", 7},
    {"kinematic energy", 8},
    {"potential energy", 9},
    {"radial kinetic energy", 10},
    {"azimuthal kinetic energy", 11},
    {"eccentricity", 12},
    {"periastron", 13},
    {"viscous dissipation", 14},
    {"luminosity", 15},
    {"pdivv", 16},
    {"delta mass inner positive", 17},
    {"delta mass inner negative", 18},
    {"delta mass outer positive", 19},
    {"delta mass outer negative", 20},
	{"delta mass inner wave damping positive", 21},
	{"delta mass inner wave damping negative", 22},
	{"delta mass outer wave damping positive", 23},
	{"delta mass outer wave damping negative", 24},
	{"delta mass floor density positive", 25},
	{"aspect ratio", 26}};

const std::map<const std::string, const std::string> quantities_file_variables =
    {{"physical time", "time"},
     {"mass", "mass"},
     {"radius", "length"},
     {"angular momentum", "angular_momentum"},
     {"total energy", "energy"},
     {"internal energy", "energy"},
     {"kinematic energy", "energy"},
     {"potential energy", "energy"},
     {"viscous dissipation", "power"},
     {"luminosity", "power"},
     {"pdivv", "pressure per time"},
     {"radial kinetic energy", "energy"},
     {"azimuthal kinetic energy", "energy"},
     {"delta mass inner positive", "mass"},
     {"delta mass inner negative", "mass"},
     {"delta mass outer positive", "mass"},
     {"delta mass outer negative", "mass"},
	 {"delta mass inner wave damping positive", "mass"},
	 {"delta mass inner wave damping negative", "mass"},
	 {"delta mass outer wave damping positive", "mass"},
	 {"delta mass outer wave damping negative", "mass"},
     {"delta mass floor density positive", "mass"},
     {"time step", "1"},
     {"analysis time step", "1"},
     {"omega frame", "frequency"},
     {"lost mass", "mass"},
     {"frame angle", "frequency"},
     {"eccentricity", "1"},
     {"periastron", "1"},
     {"aspect ratio", "1"}};

const auto quantities_file_column = quantities_file_column_v2_3;

void check_free_space(t_data &data)
{
    char *directory_name;
    DIR *directory_pointer;
    struct statvfs fiData;

    if (asprintf(&directory_name, "%s/", OUTPUTDIR) < 0) {
	die("Not enough memory.");
    }

    // Create output directory if it doesn't exist
    if (CPU_Master) {
	struct stat buffer;
	if (stat(OUTPUTDIR, &buffer)) {
	    mkdir(OUTPUTDIR, 0700);
	}
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // check if output directory exists
    if ((directory_pointer = opendir(directory_name)) == NULL) {
	logging::print_master(LOG_ERROR "Output directory %s doesn't exist!\n",
			      OUTPUTDIR);
	die("Not output directory!");
    }

    closedir(directory_pointer);

    unsigned long int space_needed = 0l;
    unsigned int number_of_files = 0l;

    // go thru all polargrids in output_polargrids and check if they are going
    // to be put out
    for (unsigned int i = 0; i < t_data::N_POLARGRID_TYPES; ++i) {
	if (data[(t_data::t_polargrid_type)i].get_write_1D()) {
	    space_needed += data[(t_data::t_polargrid_type)i].bytes_needed_1D();
	    number_of_files += 1;
	}

	if (data[(t_data::t_polargrid_type)i].get_write_2D()) {
	    space_needed += data[(t_data::t_polargrid_type)i].bytes_needed_2D();
	    number_of_files += 1;
	}
    }

    space_needed *= NTOT / NINTERM;
    number_of_files *= NTOT / NINTERM;

    logging::print_master(LOG_INFO "Output information:\n");
    logging::print_master(LOG_INFO "   Output directory: %s\n", OUTPUTDIR);
    logging::print_master(LOG_INFO "    Number of files: %u\n",
			  number_of_files);
    logging::print_master(LOG_INFO "  Total output size: %.2f GB\n",
			  (double)space_needed / 1024.0 / 1024.0 / 1024.0);

    if ((statvfs(directory_name, &fiData)) < 0) {
	logging::print_master(
	    LOG_WARNING
	    "Couldn't stat filesystem. You have to check for enough free space manually!\n");
    } else {
	// free space of device, more precisely number of space available to
	// non-priv processes (like us)
	unsigned long int free_space = fiData.f_bavail * fiData.f_frsize;

	logging::print_master(LOG_INFO "    Space Available: %.2f GB\n",
			      (double)free_space / 1024.0 / 1024.0 / 1024.0);

	if (space_needed > free_space) {
	    logging::print_master(
		LOG_WARNING
		"There is not enough space for all outputs! The program will fail at same point!\n");
	}
    }

    free(directory_name);
}

void write_grids(t_data &data, int index, int iter, double phystime,
		 bool debug = false)
{
    if (!debug) {
	logging::print_master(
	    LOG_INFO
	    "Writing output %d, Timestep Number %d, Physical Time %f.\n",
	    index, iter, phystime);
    }

    // go thru all grids and write them
    for (unsigned int i = 0; i < t_data::N_POLARGRID_TYPES; ++i) {
	data[(t_data::t_polargrid_type)i].write_polargrid(index, data, debug);
    }

    if (!debug) {
	// go thru all grids and write them
	for (unsigned int i = 0; i < t_data::N_RADIALGRID_TYPES; ++i) {
	    data[(t_data::t_radialgrid_type)i].write_radialgrid(index, data);
	}
    }
}

/**

*/
void write_quantities(t_data &data, bool force_update)
{
    FILE *fd = 0;
    char *fd_filename;
    static bool fd_created = false;

    if (CPU_Master) {

	if (asprintf(&fd_filename, "%s%s", OUTPUTDIR, "Quantities.dat") == -1) {
	    logging::print_master(LOG_ERROR
				  "Not enough memory for string buffer.\n");
	    PersonalExit(1);
	}
	// check if file exists and we restarted
	if ((start_mode::mode == start_mode::mode_restart) && !(fd_created)) {
	    if (access(fd_filename, W_OK) != -1) {
		fd_created = true;
	    }
	}

	// open logfile
	if (!fd_created) {
	    fd = fopen(fd_filename, "w");
	} else {
	    fd = fopen(fd_filename, "a");
	}
	if (fd == NULL) {
	    logging::print_master(
		LOG_ERROR "Can't write 'Quantities.dat' file. Aborting.\n");
	    PersonalExit(1);
	}

	free(fd_filename);

	if (!fd_created) {
	    // print header
	    fprintf(fd, "#FargoCPT quantities file\n");
	    fprintf(fd, "#version: 2.2\n");
	    fprintf(fd, "%s",
		    text_file_variable_description(quantities_file_column,
						   quantities_file_variables)
			.c_str());
	    fd_created = true;
	}
    }

    double quantities_limit_radius;
    if (quantities_radius_limit < 0.0) {
	auto &primary = data.get_planetary_system().get_planet(0);
	// distance to primary is distance to secondary for the primary
	quantities_limit_radius = primary.get_distance_to_primary() *
				  primary.get_dimensionless_roche_radius();
    } else {
	quantities_limit_radius = quantities_radius_limit;
    }

    const auto disk_quantities = reduce_disk_quantities(
	data, N_output, force_update, quantities_limit_radius);
    const double disk_eccentricity = disk_quantities[0];
    const double disk_periastron = disk_quantities[1];

    // computate absolute deviation from start values (this has to be done on
    // all nodes!)

    const double totalMass =
	quantities::gas_total_mass(data, quantities_limit_radius);

    if (totalMass <= 0.0) { // If roche lobe is smaller than RMIN
	quantities::gas_total_mass(data, RMAX);
    }

    const double diskRadius = quantities::gas_disk_radius(data, totalMass);
    const double totalAngularMomentum =
	quantities::gas_angular_momentum(data, quantities_limit_radius);
    const double internalEnergy =
	quantities::gas_internal_energy(data, quantities_limit_radius);
    const double qplus =
	quantities::gas_viscous_dissipation(data, quantities_limit_radius);
    const double qminus =
	quantities::gas_luminosity(data, quantities_limit_radius);
    const double kinematicEnergy =
	quantities::gas_kinematic_energy(data, quantities_limit_radius);
    const double radialKinematicEnergy =
	quantities::gas_radial_kinematic_energy(data, quantities_limit_radius);
    const double azimuthalKinematicEnergy =
	quantities::gas_azimuthal_kinematic_energy(data,
						   quantities_limit_radius);

    if (!parameters::body_force_from_potential) {
	CalculateNbodyPotential(data, PhysicalTime);
    }

	const double gravitationalEnergy =
	-quantities::gas_reduce_mass_average(data, data[t_data::POTENTIAL], quantities_limit_radius);

    const double totalEnergy =
	internalEnergy + kinematicEnergy + gravitationalEnergy;

	if(!(parameters::heating_star_enabled || parameters::self_gravity)){
	quantities::compute_aspectratio(data, N_output, force_update);
	}
	const double scale_height =
	quantities::gas_reduce_mass_average(data, data[t_data::ASPECTRATIO], quantities_limit_radius);

    double pdivv_total = 0.0;
    double InnerPositive = 0.0;
    double InnerNegative = 0.0;
    double OuterPositive = 0.0;
    double OuterNegative = 0.0;
	double InnerWaveDampingPositive = 0.0;
	double OuterWaveDampingPositive = 0.0;
	double InnerWaveDampingNegative = 0.0;
	double OuterWaveDampingNegative = 0.0;
    double FloorPositive = 0.0;

    MPI_Reduce(&data.pdivv_total, &pdivv_total, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);
    MPI_Reduce(&MassDelta.InnerPositive, &InnerPositive, 1, MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);
    MPI_Reduce(&MassDelta.InnerNegative, &InnerNegative, 1, MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);
    MPI_Reduce(&MassDelta.OuterPositive, &OuterPositive, 1, MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);
    MPI_Reduce(&MassDelta.OuterNegative, &OuterNegative, 1, MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);
	MPI_Reduce(&MassDelta.InnerWaveDampingPositive, &InnerWaveDampingPositive, 1,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&MassDelta.OuterWaveDampingPositive, &OuterWaveDampingPositive, 1,
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&MassDelta.InnerWaveDampingNegative, &InnerWaveDampingNegative, 1,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&MassDelta.OuterWaveDampingNegative, &OuterWaveDampingNegative, 1,
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&MassDelta.FloorPositive, &FloorPositive, 1, MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);

    if (CPU_Master) {
	// print to logfile
	fprintf(
	    fd,
		"%u\t%u\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\n",
	    N_output, N_outer_loop, PhysicalTime, totalMass, diskRadius,
	    totalAngularMomentum, totalEnergy, internalEnergy, kinematicEnergy,
	    gravitationalEnergy, radialKinematicEnergy,
	    azimuthalKinematicEnergy, disk_eccentricity, disk_periastron, qplus,
	    qminus, pdivv_total, InnerPositive, InnerNegative, OuterPositive,
		OuterNegative, InnerWaveDampingPositive, InnerWaveDampingNegative,
		OuterWaveDampingPositive, OuterWaveDampingNegative,
	    FloorPositive, scale_height);

	// close file
	fclose(fd);
    }
    // set mass delta to 0
    MassDelta.reset();
}

/**
	log misc. data
*/

void write_misc(const bool debug_file)
{
    if (!CPU_Master) {
	return;
    }

    std::ofstream wf;

    std::string filename;
    if (debug_file) {
	filename = std::string(OUTPUTDIR) + "debugmisc.bin";
    } else {
	filename = std::string(OUTPUTDIR) + "misc.bin";
    }

    static bool fd_created = false;

    // check if file exists and we restarted
    if ((start_mode::mode == start_mode::mode_restart) && (!fd_created)) {
	wf = std::ofstream(filename.c_str(), std::ios::in | std::ios::binary);
	if (wf.good()) {
	    fd_created = true;
	}
	wf.close();
    }

    // open logfile
    if (!fd_created || debug_file) {
	wf = std::ofstream(filename.c_str(), std::ios::out | std::ios::binary);
	if (!fd_created) {
	    fd_created = true;
	}
    } else {
	wf = std::ofstream(filename.c_str(),
			   std::ios::out | std::ios::binary | std::ios::app);
    }

    if (!wf.is_open()) {
	logging::print_master(
	    LOG_ERROR "Can't write '%s' file in \"write_misc\". Aborting.\n",
	    filename.c_str());
	PersonalExit(1);
    }

    misc_entry misc{N_output,	N_outer_loop, PhysicalTime, OmegaFrame,
		    FrameAngle, dtemp,	      last_dt,	    N_hydro_iter};

    wf.write((char *)&misc, sizeof(misc));

    wf.close();
}

std::string get_version(std::string filename)
{
    std::string version;
    std::ifstream infile(filename);

    if (infile.fail()) {
	logging::print_master(LOG_ERROR "Error: File %s cannot be opened!\n",
			      filename.c_str());
	PersonalExit(1);
    }

    std::string line_start;
    while (infile >> line_start) {
	// is it the version line
	if (line_start == "#version:") {
	    infile >> version;
	    return version;
	} else if (line_start.substr(0, 1) == "#") {
	    // were still in the header, jump to next line
	    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	} else {
	    // not a comment line anymore, so nothing was found
	    break;
	}
    }
    // nothing found up until here, assume dafault version 1.
    return "1";
}

static std::string unit_descriptor(double value, std::string unit)
{
    // produce a string containing the pair of value and unit as
    // a string such as '1.7823468234...e16 g'
    // i.e. the number with format #.16e
    std::stringstream us;
    us.precision(16);
    us << std::scientific << value << " " << unit;
    return us.str();
}

std::string text_file_variable_description(
    const std::map<const std::string, const int> &variables,
    const std::map<const std::string, const std::string> &units)
{
    // construct a header string describing each variable in
    // its own line including the column and its unit. e.g.
    // #variable: 1 | PhysicalTime | s

    std::map<int, std::string> vars_by_column;
    for (auto const &ent : variables) {
	std::string name = ent.first;
	int column = ent.second;
	vars_by_column[column] = name;
    }

    // build a map with all code units
    // need to do this runtime since units can be set in
    // parameter file
    // could also do this somewhere in initialization...
    std::map<std::string, std::string> unit_descriptors = {
	{"mass", units::mass.get_cgs_factor_symbol()},
	{"mass accretion rate",
	 units::mass_accretion_rate.get_cgs_factor_symbol()},
	{"angular_momentum", units::angular_momentum.get_cgs_factor_symbol()},
	{"time", units::time.get_cgs_factor_symbol()},
	{"energy", units::energy.get_cgs_factor_symbol()},
	{"frequency", unit_descriptor(1.0 / units::time, "1/s")},
	{"1", "1"},
	{"length", units::length.get_cgs_factor_symbol()},
	{"velocity", units::velocity.get_cgs_factor_symbol()},
	{"power", units::power.get_cgs_factor_symbol()},
	{"specific power", unit_descriptor(units::power.get_cgs_factor() /
					       units::length.get_cgs_factor() *
					       units::length.get_cgs_factor(),
					   "erg cm2/s/g")},
	{"potential", units::potential.get_cgs_factor_symbol()},
	{"pressure per time", unit_descriptor(units::pressure.get_cgs_factor() /
						  units::time.get_cgs_factor(),
					      "dyn/cm/s")},
	{"torque", unit_descriptor(units::torque.get_cgs_factor(),
				   units::torque.get_cgs_symbol())}};

    std::string var_descriptor;
    for (auto const &ent : vars_by_column) {
	std::string column = std::to_string(ent.first);
	std::string name = ent.second;
	std::string unit = units.at(name);
	var_descriptor += "#variable: " + column + " | " + name + " | " +
			  unit_descriptors[unit] + "\n";
    }
    return var_descriptor;
}

double get_from_ascii_file(std::string filename, unsigned int timestep,
			   unsigned int column, bool debug_restart)
{
    unsigned int line_timestep = 0;
    std::ifstream infile(filename);
    std::string line_start;

    if (!infile.is_open()) {
	die("Error: could not open %s\n", filename.c_str());
    }

    if (debug_restart) { // just jump to the last line of the file
	infile.seekg(-1, std::ios_base::end);
	if (infile.peek() == '\n') {
	    infile.seekg(-1, std::ios_base::cur);
	    for (int i = infile.tellg(); i > 0; i--) {
		if (infile.peek() == '\n') {
		    // Found
		    infile.get();
		    break;
		}
		infile.seekg(i, std::ios_base::beg);
	    }
	}
	infile >> line_start; // read fist element, same as non debug version.
			      // Otherwise column is not correct.
    } else {
	while (infile >> line_start) {
	    // search the file until the correct timestep is found
	    if (line_start.substr(0, 1) == "#") {
		// jump to next line
		infile.ignore(std::numeric_limits<std::streamsize>::max(),
			      '\n');
	    } else {
		// check the timestep
		line_timestep = std::stoul(line_start);
		if (line_timestep == timestep) {
		    break;
		} else {
		    // jump to next line
		    infile.ignore(std::numeric_limits<std::streamsize>::max(),
				  '\n');
		}
	    }
	}
    }
    double rv = std::nan("1");
    // read as many times as needed to reach the desired value
    for (unsigned int i = 0; i < column; i++) {
	infile >> rv;
    }

    return rv;
}

int get_misc(const int timestep, const bool debug)
{
    std::string filename;
    if (debug) {
	filename = std::string(OUTPUTDIR) + "debugmisc.bin";
    } else {
	filename = std::string(OUTPUTDIR) + "misc.bin";
    }

    std::ifstream rf(filename, std::ios::in | std::ios::binary);

    if (!rf.is_open()) {
	logging::print_master(
	    LOG_ERROR "Can't read '%s' file in \"get_misc\". Aborting.\n",
	    filename.c_str());
	PersonalExit(1);
    }

    misc_entry misc;

    rf.read((char *)&misc, sizeof(misc));
    if (!debug) {
	while (misc.timestep != timestep && !rf.eof()) {
	    if (rf.eof()) {
		logging::print(LOG_ERROR
			       "Can't read %s at timestep %d. Aborting.\n",
			       filename.c_str(), timestep);
		die("End\n");
	    }
	    rf.read((char *)&misc, sizeof(misc_entry));
	}
	if (timestep != misc.timestep) {
	    logging::print(LOG_ERROR
			   "Can't find timestep %d in %s. Aborting.\n",
			   timestep, filename.c_str());
	    die("End\n");
	}
    }

    N_output = misc.timestep;
    N_outer_loop = misc.nTimeStep;
    PhysicalTime = misc.PhysicalTime;
    OmegaFrame = misc.OmegaFrame;
    FrameAngle = misc.FrameAngle;
    dtemp = misc.dtemp;
    last_dt = misc.last_dt;
    N_hydro_iter = misc.N_iter;

    rf.close();
    return misc.timestep;
}

void write_torques(t_data &data, unsigned int timestep, bool force_update)
{
    const int ns = data[t_data::DENSITY].Nsec;
    const double *cell_center_x = CellCenterX->Field;
    const double *cell_center_y = CellCenterY->Field;

    // do everything for all planets/stars
    for (unsigned int n_planet = 0;
	 n_planet < data.get_planetary_system().get_number_of_planets();
	 ++n_planet) {
	t_planet &planet = data.get_planetary_system().get_planet(n_planet);

	const double x = planet.get_x();
	const double y = planet.get_y();
	const double mass = planet.get_mass();

	// calculate smoothing length only once if not dependend on radius

	for (unsigned int n_radial = Zero_or_active; n_radial < Max_or_active;
	     ++n_radial) {

	    data[t_data::TORQUE_1D](n_radial) = 0.0;
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
		 ++n_azimuthal) {
		// calculate smoothing length if dependend on radius
		// i.e. for thickness smoothing with scale height at cell
		// location

		const double smooth =
		    compute_smoothing(data, n_radial, n_azimuthal);
		const int cell_id = n_azimuthal + n_radial * ns;
		const double xc = cell_center_x[cell_id];
		const double yc = cell_center_y[cell_id];
		const double cellmass =
		    Surf[n_radial] *
		    data[t_data::DENSITY](n_radial, n_azimuthal);
		const double dx = xc - x;
		const double dy = yc - y;
		const double dist_sm_2 =
		    std::pow(dx, 2) + std::pow(dy, 2) + std::pow(smooth, 2);
		const double dist_sm_3 = dist_sm_2 * std::sqrt(dist_sm_2);
		const double inv_dist_sm_3 = 1.0 / dist_sm_3;

		const double Fx =
		    constants::G * cellmass * dx * inv_dist_sm_3 * mass;
		const double Fy =
		    constants::G * cellmass * dy * inv_dist_sm_3 * mass;

		const double Torque = x * Fy - y * Fx;

		data[t_data::TORQUE](n_radial, n_azimuthal) = Torque;
		data[t_data::TORQUE_1D](n_radial) += Torque;
	    }
	}

	char *name;
	if (asprintf(&name, "1D_torque_planet%i_", n_planet) < 0) {
	    die("Not enough memory!");
	}
	data[t_data::TORQUE_1D].set_name(name);
	free(name);

	if (force_update == false) {
	    data[t_data::TORQUE_1D].write1D(timestep);
	}
    }
}

void write_1D_info(t_data &data)
{
    for (int i = 0; i < t_data::N_POLARGRID_TYPES; ++i) {
	if (data[t_data::t_polargrid_type(i)].get_write_1D()) {
	    char *tmp;

	    if (asprintf(&tmp, "%s/gas%s1D.info", OUTPUTDIR,
			 data[t_data::t_polargrid_type(i)].get_name()) < 0) {
		die("Not enough memory!");
	    }

	    int Nr = GlobalNRadial;

	    if (data[t_data::t_polargrid_type(i)].is_vector())
		Nr += 1;

	    const std::string filename_info = std::string(tmp);
	    std::ofstream info_ofs(filename_info);
	    info_ofs << "# version 0.1" << std::endl;
	    info_ofs
		<< "# " << data[t_data::t_polargrid_type(i)].get_name()
		<< " 1d radial, in first line alternating: radii | quantity | minimum quantity | maximum quantity"
		<< std::endl;
	    info_ofs << "# values at time in timestepCoarse.dat" << std::endl;
	    info_ofs << "Nr = " << Nr << std::endl;
	    std::string unit;
	    if (data[t_data::t_polargrid_type(i)].get_unit() != NULL) {
		unit = std::string(data[t_data::t_polargrid_type(i)]
				       .get_unit()
				       ->get_cgs_symbol());
	    } else {
		unit = "1";
	    }
	    info_ofs << "unit = " << unit << std::endl;
	    if (data[t_data::t_polargrid_type(i)].get_unit() != NULL) {
		info_ofs << "code_units_to_cgs_factor = "
			 << data[t_data::t_polargrid_type(i)]
				.get_unit()
				->get_cgs_factor()
			 << std::endl;
	    } else {
		info_ofs << "code_units_to_cgs_factor = " << 1.0 << std::endl;
	    }
	    info_ofs << "bigendian = " << is_big_endian() << std::endl;
	    info_ofs.close();

	    free(tmp);
	}
    }
}

/**
	Calculates eccentricity, semi major axis and periastron averaged with
   the radial cells
*/
std::vector<double> reduce_disk_quantities(t_data &data, unsigned int timestep,
					   bool force_update,
					   const double quantitiy_radius)
{
    double local_eccentricity = 0.0;
    double disk_eccentricity = 0.0;
	double local_mass = 0.0;
	double global_mass = 0.0;
    // double semi_major_axis = 0.0;
    // double local_semi_major_axis = 0.0;
	double cell_mass = 0.0;
    double periastron = 0.0;
    double local_periastron = 0.0;

    // calculate eccentricity, semi_major_axis and periastron grid
	quantities::calculate_disk_ecc_peri(data, timestep, force_update);

    // Loop thru all cells excluding GHOSTCELLS & CPUOVERLAP cells (otherwise
    // they would be included twice!)
    for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size; ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
	     ++n_azimuthal) {
	    if (Rmed[n_radial] <= quantitiy_radius) {
		// eccentricity and semi major axis weighted with cellmass
		cell_mass = data[t_data::DENSITY](n_radial, n_azimuthal) *
			     Surf[n_radial];
		local_eccentricity +=
		    data[t_data::ECCENTRICITY](n_radial, n_azimuthal) *
			cell_mass;
		// local_semi_major_axis +=
		// data[t_data::SEMI_MAJOR_AXIS](n_radial, n_azimuthal) *
		// local_mass;
		local_periastron +=
		    data[t_data::PERIASTRON](n_radial, n_azimuthal) *
			cell_mass;
		local_mass += cell_mass;
	    }
	}
    }

    // synchronize threads
    MPI_Reduce(&local_eccentricity, &disk_eccentricity, 1, MPI_DOUBLE, MPI_SUM,
	       0, MPI_COMM_WORLD);
    // MPI_Allreduce(&local_semi_major_axis, &semi_major_axis, 1, MPI_DOUBLE,
    // MPI_SUM, MPI_COMM_WORLD);
    MPI_Reduce(&local_periastron, &periastron, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);

	MPI_Reduce(&local_mass, &global_mass, 1, MPI_DOUBLE, MPI_SUM, 0,
		  MPI_COMM_WORLD);

	if (global_mass > 0.0) {
	disk_eccentricity /= global_mass;
	periastron /= global_mass;
    } else {
	disk_eccentricity = 0.0;
	periastron = 0.0;
    }

    std::vector<double> rv = {disk_eccentricity, periastron};

    return rv;
}

void write_lightcurves(t_data &data, unsigned int timestep, bool force_update)
{
    // calculate luminosity
    quantities::calculate_radial_luminosity(data, timestep, force_update);

    // calculate dissipation
    quantities::calculate_radial_dissipation(data, timestep, force_update);

    double *luminosity_values =
	new double[parameters::lightcurves_radii.size()];
    double *dissipation_values =
	new double[parameters::lightcurves_radii.size()];

    for (unsigned int i = 0; i < parameters::lightcurves_radii.size(); ++i) {
	luminosity_values[i] = 0.0;
	dissipation_values[i] = 0.0;
    }

    // all cpu except the lowest one expects data
    if (CPU_Rank > 0) {
	MPI_Recv(luminosity_values, parameters::lightcurves_radii.size(),
		 MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, NULL);
	MPI_Recv(dissipation_values, parameters::lightcurves_radii.size(),
		 MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, NULL);
    }

    unsigned int current_lightcurves_bin = 0;
    for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size; ++n_radial) {
	while ((current_lightcurves_bin <
		parameters::lightcurves_radii.size() - 1) &&
	       (parameters::lightcurves_radii[current_lightcurves_bin] <
		Rmed[n_radial])) {
	    current_lightcurves_bin++;
	}
	luminosity_values[current_lightcurves_bin] +=
	    data[t_data::LUMINOSITY_1D](n_radial);
	dissipation_values[current_lightcurves_bin] +=
	    data[t_data::DISSIPATION_1D](n_radial);
    }

    // all cpu except the hightest one sends data
    if (CPU_Rank < CPU_Highest) {
	MPI_Send(luminosity_values, parameters::lightcurves_radii.size(),
		 MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD);
	MPI_Send(dissipation_values, parameters::lightcurves_radii.size(),
		 MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD);
    }

    // the last process can write the data
    if (CPU_Rank == CPU_Highest) {
	// write luminosities
	FILE *fd = 0;
	char *fd_filename;
	static bool fd_created_luminosity = false;

	if (asprintf(&fd_filename, "%s%s", OUTPUTDIR, "luminosity.dat") == -1) {
	    logging::print_master(LOG_ERROR
				  "Not enough memory for string buffer.\n");
	    PersonalExit(1);
	}

	// check if file exists and we restarted
	if ((start_mode::mode == start_mode::mode_restart) &&
	    !(fd_created_luminosity)) {
	    fd = fopen(fd_filename, "r");
	    if (fd) {
		fd_created_luminosity = true;
	    }
	    fclose(fd);
	}
	// open logfile
	if (!fd_created_luminosity) {
	    fd = fopen(fd_filename, "w");
	} else {
	    fd = fopen(fd_filename, "a");
	}
	if (fd == NULL) {
	    logging::print_master(
		LOG_ERROR "Can't write 'luminosity.dat' file. Aborting.\n");
	    PersonalExit(1);
	}

	free(fd_filename);

	if (!fd_created_luminosity) {
	    // print header
	    fprintf(fd, "# PhysicalTime\tluminosities\n");
	    fd_created_luminosity = true;
	}

	fprintf(fd, "%.20e\t", PhysicalTime);

	for (unsigned int i = 1; i < parameters::lightcurves_radii.size();
	     ++i) {
	    fprintf(fd, "%.20e\t", luminosity_values[i]);
	}
	fprintf(fd, "\n");

	// close file
	fclose(fd);

	// write dissipation
	static bool fd_created_dissipation = false;

	if (asprintf(&fd_filename, "%s%s", OUTPUTDIR, "dissipation.dat") ==
	    -1) {
	    logging::print_master(LOG_ERROR
				  "Not enough memory for string buffer.\n");
	    PersonalExit(1);
	}

	// check if file exists and we restarted
	if ((start_mode::mode == start_mode::mode_restart) &&
	    !(fd_created_dissipation)) {
	    fd = fopen(fd_filename, "r");
	    if (fd) {
		fd_created_dissipation = true;
		fclose(fd);
	    }
	}
	// open logfile
	if (!fd_created_dissipation) {
	    fd = fopen(fd_filename, "w");
	} else {
	    fd = fopen(fd_filename, "a");
	}
	if (fd == NULL) {
	    logging::print_master(
		LOG_ERROR "Can't write 'dissipation.dat' file. Aborting.\n");
	    PersonalExit(1);
	}

	free(fd_filename);

	if (!fd_created_dissipation) {
	    // print header
	    fprintf(fd, "# PhysicalTime\tdissipation\n");
	    fd_created_dissipation = true;
	}

	fprintf(fd, "%.20e\t", PhysicalTime);

	for (unsigned int i = 1; i < parameters::lightcurves_radii.size();
	     ++i) {
	    fprintf(fd, "%.20e\t", dissipation_values[i]);
	}
	fprintf(fd, "\n");

	// close file
	fclose(fd);
    }

    delete[] luminosity_values;
    delete[] dissipation_values;
}

/**
Write for each coarse output step the corresponding fine grained output number
and the simulation time in cgs units.
*/
void write_coarse_time(unsigned int coarseOutputNumber,
		       unsigned int fineOutputNumber)
{
    FILE *fd = 0;
    char *fd_filename;
    static bool fd_created = false;

    if (CPU_Master) {

	if (asprintf(&fd_filename, "%s%s", OUTPUTDIR, "timeCoarse.dat") == -1) {
	    logging::print_master(LOG_ERROR
				  "Not enough memory for string buffer.\n");
	    PersonalExit(1);
	}
	// check if file exists and we restarted
	if ((start_mode::mode == start_mode::mode_restart) && !(fd_created)) {
	    fd = fopen(fd_filename, "r");
	    if (fd) {
		fd_created = true;
		fclose(fd);
	    }
	}

	// open logfile
	if (!fd_created) {
	    fd = fopen(fd_filename, "w");
	} else {
	    fd = fopen(fd_filename, "a");
	}
	if (fd == NULL) {
	    logging::print_master(
		LOG_ERROR "Can't write 'timeCoarse.dat' file. Aborting.\n");
	    PersonalExit(1);
	}

	free(fd_filename);

	if (!fd_created) {
	    // print header
	    fprintf(
		fd,
		"# Time log for course output.\n"
		"# One DT is %.18g (code) and %.18g (cgs).\n"
		"# Syntax: coarse output step <tab> fine output step <tab> physical time (cgs)\n",
		DT, DT * units::time.get_cgs_factor());
	    fd_created = true;
	}
    }

    if (CPU_Master) {
	fprintf(fd, "%u\t%u\t%#.16e\n", coarseOutputNumber, fineOutputNumber,
		PhysicalTime * units::time);
	fclose(fd);
    }
}

void write_ecc_peri_changes(const unsigned int coarseOutputNumber, const unsigned fineOutputNumber)
{
	FILE *fd = 0;
	char *fd_filename;
	static bool fd_created = false;

	if (CPU_Master) {

	if (asprintf(&fd_filename, "%s%s", OUTPUTDIR, "eccentricity_change.dat") == -1) {
		logging::print_master(LOG_ERROR
				  "Not enough memory for string buffer.\n");
		PersonalExit(1);
	}
	// check if file exists and we restarted
	if ((start_mode::mode == start_mode::mode_restart) && !(fd_created)) {
		fd = fopen(fd_filename, "r");
		if (fd) {
		fd_created = true;
		fclose(fd);
		}
	}

	// open logfile
	if (!fd_created) {
		fd = fopen(fd_filename, "w");
	} else {
		fd = fopen(fd_filename, "a");
	}
	if (fd == NULL) {
		logging::print_master(
		LOG_ERROR "Can't write 'eccentricity_change.dat' file. Aborting.\n");
		PersonalExit(1);
	}

	free(fd_filename);

	if (!fd_created) {
		// print header
		fprintf(
		fd,
		"# Different torques and eccentricity changes by update steps.\n"
		"# Time unit is %.18g (cgs).\n"
		"# Torque unit is %.18g (cgs).\n"
		"# Eccentricity unit is 1.\n"
		"# Syntax:\n"
		"# 0 : coarse output step\n"
		"# 1 : fine output step\n"
		"# 2 : PhysicalTime\n"

		"# 3 : ecc change from source terms\n"
		"# 4 : ecc change from artificial viscosity\n"
		"# 5 : ecc change from viscosity\n"
		"# 6 : ecc change from transport\n"
		"# 7 : ecc change from damping\n"

		"# 8 : Periastron change from source terms\n"
		"# 9 : Periastron change from artificial viscosity\n"
		"# 10: Periastron change from from viscosity\n"
		"# 11: Periastron change from transport\n"
		"# 12: Periastron change from damping\n",
		units::time.get_cgs_factor(),
		units::torque.get_cgs_factor());
		fd_created = true;
	}
	}

	if (CPU_Master) {
	fprintf(fd, "%u\t%u\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\n",
		coarseOutputNumber,
		fineOutputNumber,
		PhysicalTime,

		delta_ecc_source,
		delta_ecc_art_visc,
		delta_ecc_visc,
		delta_ecc_transport,
		delta_ecc_damp,

		delta_peri_source,
		delta_peri_art_visc,
		delta_peri_visc,
		delta_peri_transport,
		delta_peri_damp);
	fclose(fd);
	}

	delta_ecc_source = 0.0;
	delta_ecc_art_visc = 0.0;
	delta_ecc_visc = 0.0;
	delta_ecc_transport = 0.0;
	delta_ecc_damp = 0.0;

	delta_peri_source = 0.0;
	delta_peri_art_visc = 0.0;
	delta_peri_visc = 0.0;
	delta_peri_transport = 0.0;
	delta_peri_damp = 0.0;
}

} // namespace output
