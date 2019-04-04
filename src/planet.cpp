#include "planet.h"
#include "util.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "LowTasks.h"
#include "output.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <map>

// define the variables in the planet data file
const std::map<const std::string, const int> planet_file_column_v1 = {
	{ "time step", 0 }
	,{ "x", 1 }
	,{ "y", 2 }
	,{ "vx", 3 }
	,{ "vy", 4 }
	,{ "mass", 5 }
	,{ "lost mass", 6 }
	,{ "physical time", 7 }
	,{ "omega frame", 8 }
	,{ "mdcp", 9 }
	,{ "exces mdcp", 10 }
	,{ "eccentricity calculated", 11 }
	,{ "angular momentum", 12 }
	,{ "semi-major axis", 13 }
	,{ "omega", 14 } };

// file version 2
const std::map<const std::string, const int> planet_file_column_v2 = {
	{ "time step", 0 }
	,{ "x", 1 }
	,{ "y", 2 }
	,{ "vx", 3 }
	,{ "vy", 4 }
	,{ "mass", 5 }
	,{ "physical time", 6 }
	,{ "omega frame", 7 }
	,{ "mdcp", 8 }
	,{ "exces mdcp", 9 }
	,{ "eccentricity calculated", 10 }
	,{ "angular momentum", 11 }
	,{ "semi-major axis", 12 }
	,{ "omega", 13 } };

auto planet_files_column = planet_file_column_v2;

const std::map<const std::string, const std::string> variable_units = {
	{ "time step", "1" }
	,{ "x", "length" }
	,{ "y", "length" }
	,{ "vx", "velocity" }
	,{ "vy", "velocity" }
	,{ "mass", "mass" }
	,{ "lost mass", "mass" }
	,{ "physical time", "time" }
	,{ "omega frame", "frequency" }
	,{ "mdcp", "mass" }
	,{ "exces mdcp", "mass" }
	,{ "eccentricity calculated", "1" }
	,{ "angular momentum", "angular_momentum" }
	,{ "semi-major axis", "length" }
	,{ "omega", "frequency" } };

/**
	set name of planet
*/
void t_planet::set_name(const char* name)
{
	// delete old name
	delete [] m_name;

	// aquire space for new name
	m_name = new char[strlen(name)+1];

	strcpy(m_name, name);
}

/**
	get planet distance to host star
*/
double t_planet::get_distance()
{
	return sqrt(pow2(get_x())+pow2(get_y()));
}

/**
	get planet semi major axis
*/
double t_planet::get_semi_major_axis()
{
	return pow2(get_angular_momentum()/get_mass()) / (constants::G * (M+get_mass())) / (1.0 - pow2(get_eccentricity()));
}

/**
	get planet period T
*/
double t_planet::get_period()
{
	return 2.0*PI*sqrt(pow3(get_semi_major_axis())/((M+get_mass())*constants::G));
}

/**
	get omega_kepler at current planet location
*/
double t_planet::get_omega()
{
	return sqrt(((M+get_mass())*constants::G)/pow3(get_distance()));
}

/**
	get angular momentum of planet
*/
double t_planet::get_angular_momentum()
{
	// j = r x p = r x mv
	return get_mass()* get_x() * get_vy() - get_mass() * get_y() * get_vx();
}

/**
	get planets eccentricity
*/
double t_planet::get_eccentricity()
{
	// distance
	double d = sqrt(pow2(get_x())+pow2(get_y()));
	// Runge-Lenz vector A = (p x L) - m * G * m * M * r/|r|;
	double A_x =  get_angular_momentum()/get_mass() * (M+get_mass()) * get_vy() - constants::G * M * pow2(M+get_mass()) * get_x()/d;
	double A_y = -get_angular_momentum()/get_mass() * (M+get_mass()) * get_vx() - constants::G * M * pow2(M+get_mass()) * get_y()/d;
	// eccentricity
	return sqrt(pow2(A_x) + pow2(A_y))/(constants::G*M*pow2(M+get_mass()));
}

void t_planet::create_planet_file()
{
	if (!CPU_Master)
		return;

	FILE *fd;
	char *filename = 0;

	std::string header_variable_description = output::text_file_variable_description(planet_files_column, variable_units);

	// create normal file
	if (asprintf(&filename, "%splanet%u.dat", OUTPUTDIR, get_planet_number()) == -1) {
		logging::print(LOG_ERROR "Not enough memory!\n");
		PersonalExit(1);
	}

	fd = fopen(filename, "w");
	free(filename);

	if (fd == NULL) {
		logging::print(LOG_ERROR "Can't write %s file. Aborting.\n", filename);
		PersonalExit(1);
	}
	fprintf(fd,"#FargoCPT planet file\n");
	fprintf(fd,"#version: 2\n");
	fprintf(fd, "%s", header_variable_description.c_str());
	fclose(fd);

	// create big file
	if (asprintf(&filename, "%sbigplanet%u.dat", OUTPUTDIR, get_planet_number()) == -1) {
		logging::print(LOG_ERROR "Not enough memory!\n");
		PersonalExit(1);
	}

	fd = fopen(filename, "w");
	free(filename);

	if (fd == NULL) {
		logging::print(LOG_ERROR "Can't write %s file. Aborting.\n", filename);
		PersonalExit(1);
	}

	fprintf(fd,"#FargoCPT planet file\n");
	fprintf(fd,"#version: 2\n");
	fprintf(fd, "%s", header_variable_description.c_str());

	fclose(fd);
}

void t_planet::write(unsigned int timestep, bool big_file)
{
	if (!CPU_Master)
		return;

	FILE *fd;
	char *filename = 0;

	// create filename
	if (asprintf(&filename, big_file ? "%sbigplanet%u.dat" : "%splanet%u.dat", OUTPUTDIR, get_planet_number()) == -1) {
		logging::print(LOG_ERROR "Not enough memory!\n");
		PersonalExit(1);
	}

	// open file
	fd = fopen(filename, "a");
	if (fd == NULL) {
		logging::print(LOG_ERROR "Can't write %s file. Aborting.\n", filename);
		PersonalExit(1);
	}
	free(filename);

	fprintf(fd, "%d\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n",
		timestep,
		get_x(),
		get_y(),
		get_vx(),
		get_vy(),
		get_mass(),
		PhysicalTime,
		OmegaFrame,
		mdcp,
		exces_mdcp,
		get_eccentricity(),
		get_angular_momentum(),
		get_semi_major_axis(),
		get_omega());

	// close file
	fclose(fd);
}

void t_planet::restart(unsigned int timestep)
{
	m_x = get_value_from_file(timestep, "x");
	m_y = get_value_from_file(timestep, "y");
	m_vx = get_value_from_file(timestep, "vx");
	m_vy = get_value_from_file(timestep, "vy");
	m_mass = get_value_from_file(timestep, "mass");
}

double t_planet::get_value_from_file(unsigned int timestep, std::string variable_name)
{
	double value;
	int column = -1;

	std::string filename = std::string(OUTPUTDIR) + "planet"
	    + std::to_string(get_planet_number()) + ".dat";

	std::string version = output::get_version(filename);
	auto variable_columns = planet_file_column_v2;

	if (version == "1") {
		variable_columns = planet_file_column_v1;
    } else if (version == "2") {
		variable_columns = planet_file_column_v2;
	} else {
		std::cerr << "Unknown version '" << version << "'for planet.dat file!" << std::endl;
	    PersonalExit(1);
	}

	// check whether the column map contains the variable
	auto iter = variable_columns.find(variable_name);
	if (iter != variable_columns.end() ) {
		column = iter->second;
	} else {
		std::cerr << "Unknown variable '" << variable_name << "' for planet.dat file version '" << version << "'\n" << std::endl;
	}

	if (column == -1) {
		std::cerr << "Something went wring in obtaining the column for " << variable_name << " for the planet data file." << std::endl;
		PersonalExit(1);
	}
	value = output::get_from_ascii_file(filename, timestep, column);
	return value;
}
