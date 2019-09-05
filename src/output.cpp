#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include <sys/statvfs.h>
#include "stress.h"
#include "util.h"
#include "global.h"
#include "LowTasks.h"
#include "constants.h"
#include "output.h"
#include "nongnu.h"
#include "logging.h"
#include "viscosity.h"
#include "parameters.h"
#include "SideEuler.h"
#include "quantities.h"
#include "options.h"

#include "unistd.h" // for access()
#include <sys/stat.h>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <limits>
#include <sstream>


namespace output {

// info on variables in misc file
const std::map<const std::string, const int> misc_file_column_v1 = {
	{ "time step", 0 }
	,{ "physical time", 1 }
	,{ "omega frame", 2 }
	,{ "lost mass", 3 }
	,{ "frame angle", 4 } };

const std::map<const std::string, const int> misc_file_column_v2 = {
	{ "time step", 0 }
	,{ "physical time", 1 }
	,{ "omega frame", 2 }
	,{ "frame angle", 3 } };

auto misc_file_columns = misc_file_column_v2;

const std::map<const std::string, const std::string> misc_file_variables = {
	{ "time step", "1" },
	{ "physical time", "time" },
	{ "omega frame", "frequency" },
	{ "lost mass", "mass" },
	{ "frame angle", "1" }
	};

const std::map<const std::string, const int> quantities_file_column_v2 = {
{ "time step", 0 },
{ "physical time", 1 },
{ "mass", 2 },
{ "angular momentum", 3 },
{ "total energy", 4 },
{ "internal energy", 5 },
{ "kinematic energy", 6 },
{ "potential energy", 7 },
{ "radial kinetic energy", 8 },
{ "azimuthal kinetic energy", 9 },
{ "eccentricity", 10 },
{ "periastron", 11},
{ "qplus", 12 },
{ "qminus", 13 },
{ "pdivv", 14 },
{ "delta mass inner positive", 15 },
{ "delta mass inner negative", 16 },
{ "delta mass outer positive", 17 },
{ "delta mass outer negative", 18 },
{ "delta mass wave damping positive", 19 },
{ "delta mass wave damping negative", 20 },
{ "delta mass floor density positive", 21 }
};

const std::map<const std::string, const int> quantities_file_column_v2_1 = {
{ "time step", 0 },
{ "analysis time step", 1},
{ "physical time", 2 },
{ "mass", 3 },
{ "angular momentum", 4 },
{ "total energy", 5 },
{ "internal energy", 6 },
{ "kinematic energy", 7 },
{ "potential energy", 8 },
{ "radial kinetic energy", 9 },
{ "azimuthal kinetic energy", 10 },
{ "eccentricity", 11 },
{ "periastron", 12},
{ "qplus", 13 },
{ "qminus", 14 },
{ "pdivv", 15 },
{ "delta mass inner positive", 16 },
{ "delta mass inner negative", 17 },
{ "delta mass outer positive", 18 },
{ "delta mass outer negative", 19 },
{ "delta mass wave damping positive", 20 },
{ "delta mass wave damping negative", 21 },
{ "delta mass floor density positive", 22 }
};


const std::map<const std::string, const int> quantities_file_column_v2_2 = {
{ "time step", 0 },
{ "analysis time step", 1},
{ "physical time", 2 },
{ "mass", 3 },
{ "radius", 4 },
{ "angular momentum", 5 },
{ "total energy", 6 },
{ "internal energy", 7 },
{ "kinematic energy", 8 },
{ "potential energy", 9 },
{ "radial kinetic energy", 10 },
{ "azimuthal kinetic energy", 11 },
{ "eccentricity", 12 },
{ "periastron", 13},
{ "qplus", 14 },
{ "qminus", 15 },
{ "pdivv", 16 },
{ "delta mass inner positive", 17 },
{ "delta mass inner negative", 18 },
{ "delta mass outer positive", 19 },
{ "delta mass outer negative", 20 },
{ "delta mass wave damping positive", 21 },
{ "delta mass wave damping negative", 22 },
{ "delta mass floor density positive", 23 }
};

const std::map<const std::string, const int> quantities_file_column_v2_3 = {
{ "time step", 0 },
{ "analysis time step", 1},
{ "physical time", 2 },
{ "mass", 3 },
{ "radius", 4 },
{ "angular momentum", 5 },
{ "total energy", 6 },
{ "internal energy", 7 },
{ "kinematic energy", 8 },
{ "potential energy", 9 },
{ "radial kinetic energy", 10 },
{ "azimuthal kinetic energy", 11 },
{ "eccentricity", 12 },
{ "periastron", 13},
{ "qplus", 14 },
{ "qminus", 15 },
{ "pdivv", 16 },
{ "delta mass inner positive", 17 },
{ "delta mass inner negative", 18 },
{ "delta mass outer positive", 19 },
{ "delta mass outer negative", 20 },
{ "delta mass wave damping positive", 21 },
{ "delta mass wave damping negative", 22 },
{ "delta mass floor density positive", 23 },
{ "aspect ratio", 24 }
};

const std::map<const std::string, const std::string> quantities_file_variables = {
{ "physical time", "time" },
{ "mass", "mass" },
{ "radius", "length" },
{ "angular momentum", "angular_momentum" },
{ "total energy", "energy" },
{ "internal energy", "energy" },
{ "kinematic energy", "energy" },
{ "potential energy", "energy" },
{ "qplus", "specific power" },
{ "qminus", "specific power" },
{ "pdivv", "pressure per time" },
{ "radial kinetic energy", "energy" },
{ "azimuthal kinetic energy", "energy" },
{ "delta mass inner positive", "mass" },
{ "delta mass inner negative", "mass" },
{ "delta mass outer positive", "mass" },
{ "delta mass outer negative", "mass" },
{ "delta mass wave damping positive", "mass" },
{ "delta mass wave damping negative", "mass" },
{ "delta mass floor density positive", "mass" },
{ "time step", "1" },
{ "analysis time step", "1" },
{ "omega frame", "frequency" },
{ "lost mass", "mass" },
{ "frame angle", "frequency" },
{ "eccentricity", "1" },
{ "periastron", "1" },
{ "aspect ratio", "1" }
};


const auto quantities_file_column = quantities_file_column_v2_3;


void check_free_space(t_data &data)
{
	char *directory_name;
	DIR *directory_pointer;
	struct dirent *directory_entry;
	struct statvfs fiData;

	if (asprintf(&directory_name, "%s/",OUTPUTDIR) <0) {
		die("Not enough memory.");
    }

    // Create output directory if it doesn't exist
    if(CPU_Master)
    {
        struct stat buffer;
        if(stat(OUTPUTDIR, &buffer))
        {
            mkdir(OUTPUTDIR, 0700);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);


	// check if output directory exists
	if ((directory_pointer = opendir(directory_name)) == NULL) {
        logging::print_master(LOG_ERROR "Output directory %s doesn't exist!\n", OUTPUTDIR);
        die("Not output directory!");
    }

	while ((directory_entry = readdir(directory_pointer))) {
		if ((strcmp("..",directory_entry->d_name) != 0) && (strcmp(".",directory_entry->d_name) != 0)) {
				logging::print_master(LOG_NOTICE "Output directory %s is not empty!\n", OUTPUTDIR);
			break;
		}
	}

	closedir(directory_pointer);

	unsigned long int space_needed = 0l;
	unsigned int number_of_files = 0l;

	// go thru all polargrids in output_polargrids and check if they are going to be put out
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

	space_needed *= NTOT/NINTERM;
	number_of_files *= NTOT/NINTERM;

	logging::print_master(LOG_INFO "Output information:\n");
	logging::print_master(LOG_INFO "   Output directory: %s\n", OUTPUTDIR);
	logging::print_master(LOG_INFO "    Number of files: %u\n", number_of_files);
	logging::print_master(LOG_INFO "  Total output size: %.2f GB\n", (double)space_needed/1024.0/1024.0/1024.0);

	if((statvfs(directory_name,&fiData)) < 0 ) {
		logging::print_master(LOG_WARNING "Couldn't stat filesystem. You have to check for enough free space manually!\n");
	} else {
		// free space of device, more precisely number of space available to non-priv processes (like us)
		unsigned long int free_space = fiData.f_bavail * fiData.f_frsize;

		logging::print_master(LOG_INFO "    Space Available: %.2f GB\n", (double)free_space/1024.0/1024.0/1024.0);

		if (space_needed > free_space) {
			logging::print_master(LOG_WARNING "There is not enough space for all outputs! The program will fail at same point!\n");
		}
	}

	free(directory_name);
}

void write_grids(t_data &data, int index, int iter, double phystime)
{
	logging::print_master(LOG_INFO "Writing output %d, Timestep Number %d, Physical Time %f.\n", index, iter, phystime);

	// go thru all grids and write them
	for (unsigned int i = 0; i < t_data::N_POLARGRID_TYPES; ++i) {
		data[(t_data::t_polargrid_type)i].write(index, data);
	}

	// go thru all grids and write them
	for (unsigned int i = 0; i < t_data::N_RADIALGRID_TYPES; ++i) {
		data[(t_data::t_radialgrid_type)i].write(index, data);
	}
}

/**

*/
void write_quantities(t_data &data, unsigned int timestep, unsigned int nTimeStep, bool force_update)
{
	FILE* fd = 0;
	char* fd_filename;
	static bool fd_created = false;

	if (CPU_Master) {

		if (asprintf(&fd_filename, "%s%s", OUTPUTDIR, "Quantities.dat") == -1) {
			logging::print_master(LOG_ERROR "Not enough memory for string buffer.\n");
			PersonalExit(1);
		}
		// check if file exists and we restarted
		if ((options::restart) && !(fd_created)) {
			if(access( fd_filename, W_OK ) != -1)
			{
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
			logging::print_master(LOG_ERROR "Can't write 'Quantities.dat' file. Aborting.\n");
			PersonalExit(1);
		}

		free(fd_filename);

		if (!fd_created) {
			// print header
			fprintf(fd,"#FargoCPT quantities file\n");
			fprintf(fd,"#version: 2.2\n");
			fprintf(fd,"%s", text_file_variable_description(quantities_file_column, quantities_file_variables).c_str() );
			fd_created = true;
		}
	}

	auto disk_quantities = reduce_disk_quantities(data, timestep, force_update);
	double disk_eccentricity = disk_quantities[0];
	double disk_periastron = disk_quantities[1];

	// computate absolute deviation from start values (this has to be done on all nodes!)
	double totalMass = quantities::gas_total_mass(data);
	double diskRadius = quantities::gas_disk_radius(data, totalMass);
	double totalAngularMomentum = quantities::gas_angular_momentum(data);
	double internalEnergy = quantities::gas_internal_energy(data);
	double qplus = quantities::gas_qplus(data);
	double qminus = quantities::gas_qminus(data);
	double kinematicEnergy = quantities::gas_kinematic_energy(data);
	double radialKinematicEnergy = quantities::gas_radial_kinematic_energy(data);
	double azimuthalKinematicEnergy = quantities::gas_azimuthal_kinematic_energy(data);
	double gravitationalEnergy = quantities::gas_gravitational_energy(data);
	double totalEnergy = internalEnergy + kinematicEnergy + gravitationalEnergy;
	double scale_height = quantities::gas_aspect_ratio(data);


	double pdivv_total = 0.0;
	double InnerPositive = 0.0;
	double InnerNegative = 0.0;
	double OuterPositive = 0.0;
	double OuterNegative = 0.0;
	double WaveDampingPositive = 0.0;
	double WaveDampingNegative = 0.0;
	double FloorPositive	= 0.0;

	MPI_Reduce(&data.pdivv_total,				&pdivv_total,			1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&MassDelta.InnerPositive,		&InnerPositive,			1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&MassDelta.InnerNegative,		&InnerNegative,			1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&MassDelta.OuterPositive,		&OuterPositive,			1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&MassDelta.OuterNegative,		&OuterNegative,			1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&MassDelta.WaveDampingPositive,	&WaveDampingPositive,	1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&MassDelta.WaveDampingNegative,	&WaveDampingNegative,	1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&MassDelta.FloorPositive,		&FloorPositive,			1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


	if (CPU_Master) {
		// print to logfile
		fprintf(fd, "%u\t%u\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\t%#.16e\n",
			timestep,
            nTimeStep,
			PhysicalTime,
			totalMass,
			diskRadius,
			totalAngularMomentum,
			totalEnergy,
			internalEnergy,
			kinematicEnergy,
			gravitationalEnergy,
			radialKinematicEnergy,
			azimuthalKinematicEnergy,
			disk_eccentricity,
			disk_periastron,
			qplus,
			qminus,
			pdivv_total,
			InnerPositive,
			InnerNegative,
			OuterPositive,
			OuterNegative,
			WaveDampingPositive,
			WaveDampingNegative,
			FloorPositive,
			scale_height);

		// close file
		fclose(fd);
	}
	// set mass delta to 0
	MassDelta.reset();
}

/**
	log misc. data
*/
void write_misc(unsigned int timestep)
{
	FILE *fd = 0;
	char* fd_filename;
	static bool fd_created = false;

	if (CPU_Master) {
		if (asprintf(&fd_filename, "%s%s", OUTPUTDIR, "misc.dat") == -1) {
			logging::print_master(LOG_ERROR "Not enough memory for string buffer.\n");
			PersonalExit(1);
		}
		// check if file exists and we restarted
		if ((options::restart) && !(fd_created)) {
			fd = fopen(fd_filename, "r");
			if (fd) {
				fd_created = true;
			}
			fclose(fd);
		}
		// open logfile
		if (!fd_created) {
			fd = fopen(fd_filename, "w");
		} else {
			fd = fopen(fd_filename, "a");
		}
		if (fd == NULL) {
			logging::print_master(LOG_ERROR "Can't write 'misc.dat' file. Aborting.\n");
			PersonalExit(1);
		}

		free(fd_filename);

		if (!fd_created) {
			// print header
			fprintf(fd,"#FargoCPT misc file\n");
			fprintf(fd,"#version: 2\n");
			fprintf(fd,"%s", text_file_variable_description(misc_file_columns, misc_file_variables).c_str());
			fd_created = true;
		}
	}

	if (CPU_Master) {
		// print to logfile
		fprintf(fd, "%u\t%#.18g\t%#.18g\t%#.18g\n", timestep, PhysicalTime, OmegaFrame, FrameAngle);

		// close file
		fclose(fd);
	}
}

std::string get_version(std::string filename) {
	std::string version;
	std::ifstream infile(filename);

	if(infile.fail()){
		logging::print_master(LOG_ERROR "Error: File %s cannot be opened!\n", filename.c_str());
		PersonalExit(1);
	}

	std::string line_start;
	while (infile >> line_start) {
		// is it the version line
		if (line_start == "#version:") {
			infile >> version;
			return version;
		} else if (line_start.substr(0,1) == "#") {
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

std::string unit_descriptor(double value, std::string unit) {
	// produce a string containing the pair of value and unit as
	// a string such as '1.7823468234...e16 g'
	// i.e. the number with format #.16e
	std::stringstream us;
	us.precision(16);
	us << std::scientific << value << " " << unit;
	return us.str();
}

std::string text_file_variable_description(const std::map<const std::string, const int> &variables, const std::map<const std::string, const std::string> &units) {
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
		{ "mass" , units::mass.get_cgs_factor_symbol() },
		{ "angular_momentum" , units::angular_momentum.get_cgs_factor_symbol() },
		{ "time" , units::time.get_cgs_factor_symbol() },
		{ "energy" , units::energy.get_cgs_factor_symbol() },
		{ "frequency" , unit_descriptor( 1.0/units::time, "1/s")  },
		{ "1" , "1" },
		{ "length" , units::length.get_cgs_factor_symbol() },
		{ "velocity" , units::velocity.get_cgs_factor_symbol() },
		{ "power", units::power.get_cgs_factor_symbol() },
		{ "specific power", unit_descriptor(
											units::power.get_cgs_factor()/units::length.get_cgs_factor()*units::length.get_cgs_factor()
											, "erg cm2/s/g" ) },
		{ "pressure per time",  unit_descriptor(
											units::pressure.get_cgs_factor()/units::time.get_cgs_factor()
											, "dyn/cm/s" ) }
	};

	std::string var_descriptor;
	for (auto const &ent : vars_by_column) {
		std::string column = std::to_string(ent.first);
		std::string name = ent.second;
		std::string unit = units.at(name);
		var_descriptor += "#variable: " + column + " | "
			+ name + " | " + unit_descriptors[unit] + "\n";
	}
	return var_descriptor;
}

double get_from_ascii_file(std::string filename, unsigned int timestep, unsigned int column) {
	unsigned int line_timestep = 0;
	std::ifstream infile(filename);
	std::string line_start;

	while (infile >> line_start) {
		// search the file until the correct timestep is found
		if (line_start.substr(0,1) == "#") {
			// jump to next line
			infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		} else {
			// check the timestep
			line_timestep = std::stoul(line_start);
			if (line_timestep == timestep) {
				break;
			} else {
				// jump to next line
				infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			}
		}
	}
	double rv;
	// read as many times as needed to reach the desired value
	for (unsigned int i=0; i<column; i++) {
		infile >> rv;
	}
	return rv;
}

double get_misc(unsigned int timestep, std::string variable)
{
	unsigned int column = 0;
	std::string filename = std::string(OUTPUTDIR) + "misc.dat";
	std::string version = get_version(filename);

	if (version == "2") {
		if (variable == "timestep") column = 0;
		else if (variable == "physical time") column = 1;
		else if (variable == "omega frame") column = 2;
		else if (variable == "frame angle") column = 3;
		else {
			printf("Don't know variable '%s' in misc.dat v2\n", variable.c_str());
			PersonalExit(1);
		}
	}
	else if (version == "1") {
		if (variable == "timestep") column = 0;
		else if (variable == "physical time") column = 1;
		else if (variable == "omega frame") column = 2;
		else if (variable == "lost mass") column = 3;
		else if (variable == "frame angle") column = 4;
		else {
			printf("Don't know variable '%s' in misc.dat v1\n", variable.c_str());
			PersonalExit(1);
		}
	}
    double rv = get_from_ascii_file(filename, timestep, column);
	return rv;
}

void write_torques(t_data &data, unsigned int timestep, bool force_update) {
	double *global_torques = (double*)malloc(sizeof(*global_torques)*(data.get_planetary_system().get_number_of_planets()+1));
	double *local_torques = (double*)malloc(sizeof(*local_torques)*(data.get_planetary_system().get_number_of_planets()+1));

	// central star
	local_torques[0] = 0.0;

	// do everything for all planets/stars
	for (unsigned int n_planet = 0; n_planet < data.get_planetary_system().get_number_of_planets(); ++n_planet) {
		local_torques[n_planet+1] = 0;
		t_planet& planet = data.get_planetary_system().get_planet(n_planet);
		double smooth = parameters::thickness_smoothing * ASPECTRATIO_REF * pow(planet.get_distance(), 1.0+FLARINGINDEX);

		// hill radius
		double r_hill = pow(planet.get_mass()/(3.0*(M+planet.get_mass())),1.0/3.0)*planet.get_semi_major_axis();
		// 0.8 = cut off parameter
		double r_taper = 0.8 * r_hill;

		for (unsigned int n_radial = 0; n_radial <= data[t_data::TORQUE].get_max_radial(); ++n_radial) {
			data[t_data::TORQUE_1D](n_radial) = 0.0;

			for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::TORQUE].get_max_azimuthal(); ++n_azimuthal) {
				double phi = (double)n_azimuthal/(double)data[t_data::TORQUE].get_size_azimuthal()*2.0*PI;
				double dx = planet.get_x() - Rmed[n_radial]*cos(phi);
				double dy = planet.get_y() - Rmed[n_radial]*sin(phi);
				double m = data[t_data::DENSITY](n_radial, n_azimuthal)*Surf[n_radial];

				double distance2 = pow2(dx) + pow2(dy) + pow2(smooth);

				double taper = 1.0;

				if (sqrt(pow2(dx) + pow2(dy)) < 4.0*r_taper) {
					taper = 1.0/(exp(-(sqrt(pow2(dx)+pow2(dy))-r_taper)/(0.1*r_taper))+1.0);
				}

				double F_x = - dx * m * planet.get_mass() * pow(distance2, -1.5) * taper;
				double F_y = - dy * m * planet.get_mass() * pow(distance2, -1.5) * taper;

				data[t_data::TORQUE](n_radial, n_azimuthal) = planet.get_x()*F_y - planet.get_y()*F_x;
				data[t_data::TORQUE_1D](n_radial) += data[t_data::TORQUE](n_radial, n_azimuthal);
				if ((n_radial >= ((CPU_Rank == 0) ? GHOSTCELLS_B : CPUOVERLAP)) && (n_radial <= data[t_data::TORQUE].get_max_radial() - ( CPU_Rank == CPU_Highest ? GHOSTCELLS_B : CPUOVERLAP ))) {
					local_torques[n_planet+1] += data[t_data::TORQUE](n_radial, n_azimuthal);
				}
			}
		}

		char* name;
		if (asprintf(&name, "1D_torque_planet%i_", n_planet)<0) {
			die("Not enough memory!");
		}
		data[t_data::TORQUE_1D].set_name(name);
		free(name);

		if (force_update == false) {
			data[t_data::TORQUE_1D].write1D(timestep);
		}
	}

	MPI_Allreduce(local_torques, global_torques, data.get_planetary_system().get_number_of_planets()+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	if (CPU_Master) {
		// output
		FILE* fd = 0;
		char* fd_filename;
		static bool fd_created = false;

		if (asprintf(&fd_filename, "%s%s", OUTPUTDIR, "torques.dat") == -1) {
			logging::print_master(LOG_ERROR "Not enough memory for string buffer.\n");
			PersonalExit(1);
		}

		// check if file exists and we restarted
		if ((options::restart) && !(fd_created)) {
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
			logging::print_master(LOG_ERROR "Can't write 'torques.dat' file. Aborting.\n");
			PersonalExit(1);
		}

		free(fd_filename);

		if (!fd_created) {
			// print header
			fprintf(fd,"# \n");
			fd_created = true;
		}

		fprintf(fd,"%.20e", PhysicalTime);
		for (unsigned int i = 0; i <= data.get_planetary_system().get_number_of_planets(); ++i) {
			fprintf(fd,"\t%.20e", global_torques[i]*units::torque);
		}
		fprintf(fd, "\n");

		// close file
		fclose(fd);
	}

	free(local_torques);
	free(global_torques);
}

void write_1D_info(t_data &data) {
	for(int i = 0; i < t_data::N_POLARGRID_TYPES; ++i)
	{
		if (data[t_data::t_polargrid_type(i)].get_write_1D())
        {
			char *tmp;

			if (asprintf(&tmp, "%s/gas%s1D.info",OUTPUTDIR,data[t_data::t_polargrid_type(i)].get_name())<0) {
				die("Not enough memory!");
			}

			int Nr = GlobalNRadial;

			if(data[t_data::t_polargrid_type(i)].is_vector())
				Nr += 1;

			const std::string filename_info = std::string(tmp);
			std::ofstream info_ofs(filename_info);
			info_ofs << "# version 0.1" << std::endl;
			info_ofs << "# " <<  data[t_data::t_polargrid_type(i)].get_name() << " 1d radial, in first line alternating: radii | quantity | minimum quantity | maximum quantity" << std::endl;
			info_ofs << "# values at time in timestepCoarse.dat" << std::endl;
			info_ofs << "Nr = " << Nr << std::endl;
			std::string unit;
			if (data[t_data::t_polargrid_type(i)].get_unit() != NULL) {
				unit = std::string( data[t_data::t_polargrid_type(i)].get_unit()->get_cgs_symbol() );
			} else {
				unit = "1";
			}
			info_ofs << "unit = " << unit  << std::endl;
			info_ofs << "bigendian = " << is_big_endian() << std::endl;
			info_ofs.close();
		}
	}
}

void write_massflow_info(t_data &data) {
	const std::string filename_info = std::string(OUTPUTDIR) + "/gasMassFlow1D.info";
	std::ofstream info_ofs(filename_info);
	info_ofs << "# Mass flow 1d radial, first line radii, from second line on, values at time in Quantities.dat" << std::endl;
	info_ofs << "Nr = " << GlobalNRadial+1 << std::endl;
	info_ofs << "unit = " << data[t_data::MASSFLOW_1D].get_unit()->get_cgs_symbol() << std::endl;
	info_ofs <<	 "bigendian = " << is_big_endian() << std::endl;
	const std::string filename = std::string(OUTPUTDIR) + "/gasMassFlow1D.dat";
	std::ofstream ofs(filename,  std::ios::binary);
	ofs.write( (char*)Radii.array, sizeof(*Radii.array)*(GlobalNRadial+1) );
}

void write_massflow(t_data &data, unsigned int timestep) {
    (void) timestep;
	const std::string filename = std::string(OUTPUTDIR) + "/gasMassFlow1D.dat";
	data[t_data::MASSFLOW_1D].write(filename, TimeStep, data, true, true);
}

/**
	Calculates eccentricity, semi major axis and periastron averaged with the radial cells
*/
std::vector<double> reduce_disk_quantities(t_data &data, unsigned int timestep, bool force_update)
{
	double local_eccentricity = 0.0;
	double gas_total_mass = quantities::gas_total_mass(data);
	double disk_eccentricity = 0.0;
	//double semi_major_axis = 0.0;
	//double local_semi_major_axis = 0.0;
	double local_mass= 0.0;
	double periastron = 0.0;
	double local_periastron = 0.0;

	// calculate eccentricity, semi_major_axis and periastron grid
	quantities::calculate_disk_quantities(data, timestep, force_update);

	// Loop thru all cells excluding GHOSTCELLS & CPUOVERLAP cells (otherwise they would be included twice!)
	for (unsigned int n_radial = (CPU_Rank == 0) ? GHOSTCELLS_B : CPUOVERLAP; n_radial <= data[t_data::DENSITY].get_max_radial() - ( CPU_Rank == CPU_Highest ? GHOSTCELLS_B : CPUOVERLAP ); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
			// eccentricity and semi major axis weighted with cellmass
			local_mass = data[t_data::DENSITY](n_radial, n_azimuthal) * Surf[n_radial];
			local_eccentricity += data[t_data::ECCENTRICITY](n_radial, n_azimuthal) * local_mass;
			//local_semi_major_axis += data[t_data::SEMI_MAJOR_AXIS](n_radial, n_azimuthal) * local_mass;
			local_periastron += data[t_data::PERIASTRON](n_radial, n_azimuthal) * local_mass;
		}
	}

	// synchronize threads
	MPI_Reduce(&local_eccentricity, &disk_eccentricity, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//MPI_Allreduce(&local_semi_major_axis, &semi_major_axis, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Reduce(&local_periastron, &periastron, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	disk_eccentricity /= gas_total_mass;
	//semi_major_axis /= gas_total_mass;
	periastron /= gas_total_mass;

	std::vector<double> rv = {disk_eccentricity, periastron};

	return rv;
}

void write_lightcurves(t_data &data, unsigned int timestep, bool force_update)
{
	// calculate luminosity
	quantities::calculate_radial_luminosity(data, timestep, force_update);

	// calculate dissipation
	quantities::calculate_radial_dissipation(data, timestep, force_update);

	double* luminosity_values = new double[parameters::lightcurves_radii.size()];
	double* dissipation_values = new double[parameters::lightcurves_radii.size()];

	for (unsigned int i = 0; i < parameters::lightcurves_radii.size(); ++i) {
		luminosity_values[i] = 0.0;
		dissipation_values[i] = 0.0;
	}

	// all cpu except the lowest one expects data
	if (CPU_Rank > 0) {
		MPI_Recv(luminosity_values, parameters::lightcurves_radii.size(), MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, NULL);
		MPI_Recv(dissipation_values, parameters::lightcurves_radii.size(), MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, NULL);
	}

	unsigned int current_lightcurves_bin = 0;
	for (unsigned int n_radial = (CPU_Rank == 0) ? GHOSTCELLS_B : CPUOVERLAP; n_radial <= data[t_data::LUMINOSITY_1D].get_max_radial() - ( CPU_Rank == CPU_Highest ? GHOSTCELLS_B : CPUOVERLAP ); ++n_radial) {
		while ((current_lightcurves_bin < parameters::lightcurves_radii.size()-1) && (parameters::lightcurves_radii[current_lightcurves_bin] < Rmed[n_radial])) {
			current_lightcurves_bin++;
		}
		luminosity_values[current_lightcurves_bin] += data[t_data::LUMINOSITY_1D](n_radial);
		dissipation_values[current_lightcurves_bin] += data[t_data::DISSIPATION_1D](n_radial);
	}

	// all cpu except the hightest one sends data
	if (CPU_Rank < CPU_Highest) {
		MPI_Send(luminosity_values, parameters::lightcurves_radii.size(), MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD);
		MPI_Send(dissipation_values, parameters::lightcurves_radii.size(), MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD);
	}

	// the last process can write the data
	if (CPU_Rank == CPU_Highest) {
		// write luminosities
		FILE* fd = 0;
		char* fd_filename;
		static bool fd_created_luminosity = false;

		if (asprintf(&fd_filename, "%s%s", OUTPUTDIR, "luminosity.dat") == -1) {
			logging::print_master(LOG_ERROR "Not enough memory for string buffer.\n");
			PersonalExit(1);
		}

		// check if file exists and we restarted
		if ((options::restart) && !(fd_created_luminosity)) {
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
			logging::print_master(LOG_ERROR "Can't write 'luminosity.dat' file. Aborting.\n");
			PersonalExit(1);
		}

		free(fd_filename);

		if (!fd_created_luminosity) {
			// print header
			fprintf(fd,"# PhysicalTime\tluminosities\n");
			fd_created_luminosity = true;
		}

		fprintf(fd,"%.20e\t",PhysicalTime);

		for (unsigned int i = 1; i < parameters::lightcurves_radii.size(); ++i) {
			fprintf(fd, "%.20e\t", luminosity_values[i]*units::power.get_cgs_factor());
		}
		fprintf(fd,"\n");

		// close file
		fclose(fd);

		// write dissipation
		static bool fd_created_dissipation = false;

		if (asprintf(&fd_filename, "%s%s", OUTPUTDIR, "dissipation.dat") == -1) {
			logging::print_master(LOG_ERROR "Not enough memory for string buffer.\n");
			PersonalExit(1);
		}

		// check if file exists and we restarted
		if ((options::restart) && !(fd_created_dissipation)) {
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
			logging::print_master(LOG_ERROR "Can't write 'dissipation.dat' file. Aborting.\n");
			PersonalExit(1);
		}

		free(fd_filename);

		if (!fd_created_dissipation) {
			// print header
			fprintf(fd,"# PhysicalTime\tdissipation\n");
			fd_created_dissipation = true;
		}

		fprintf(fd,"%.20e\t",PhysicalTime);

		for (unsigned int i = 1; i < parameters::lightcurves_radii.size(); ++i) {
			fprintf(fd, "%.20e\t", dissipation_values[i]*units::power.get_cgs_factor());
		}
		fprintf(fd,"\n");

		// close file
		fclose(fd);
	}

    delete[] luminosity_values;
    delete[] dissipation_values;
}

/**
Write for each coarse output step the corresponding fine grained output number and the simulation time in cgs units.
*/
void write_coarse_time(unsigned int coarseOutputNumber, unsigned int fineOutputNumber)
{
	FILE* fd = 0;
	char* fd_filename;
	static bool fd_created = false;

	if (CPU_Master) {

		if (asprintf(&fd_filename, "%s%s", OUTPUTDIR, "timeCoarse.dat") == -1) {
			logging::print_master(LOG_ERROR "Not enough memory for string buffer.\n");
			PersonalExit(1);
		}
		// check if file exists and we restarted
		if ((options::restart) && !(fd_created)) {
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
			logging::print_master(LOG_ERROR "Can't write 'timeCoarse.dat' file. Aborting.\n");
			PersonalExit(1);
		}

		free(fd_filename);

		if (!fd_created) {
			// print header
			fprintf(fd,"# Time log for course output.\n# Syntax: coarse output step <tab> fine output step <tab> physical time (cgs)\n");
			fd_created = true;
		}
	}

	if (CPU_Master) {
		fprintf(fd, "%u\t%u\t%#.16e\n"
				, coarseOutputNumber
				, fineOutputNumber
				, PhysicalTime*units::time);
		fclose(fd);
	}
}

} // close namespace
