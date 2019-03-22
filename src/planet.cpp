#include "planet.h"
#include "util.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "LowTasks.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

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
	fprintf(fd,"#content: TimeStep\tXplanet\tYplanet\tVXplanet\tVYplanet\tMplanetVirtual\tPhysicalTime\tOmegaFrame\tmdcp\texces_mdcp\teccentricity_calculated\tj\tsemi_major_axis\tomega\n");

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
	m_x = get_value_from_file(timestep, 2);
	m_y = get_value_from_file(timestep, 3);
	m_vx = get_value_from_file(timestep, 4);
	m_vy = get_value_from_file(timestep, 5);
	m_mass = get_value_from_file(timestep, 6);
}

double t_planet::get_value_from_file(unsigned int timestep, unsigned int column)
{
	FILE *fd;
	char *filename = 0;
	char buffer[256];
	char *ptr;
	unsigned int line_timestep;
	double value;

	// create filename
	if (asprintf(&filename, "%splanet%u.dat", OUTPUTDIR, get_planet_number()) == -1) {
		logging::print(LOG_ERROR "Not enough memory!\n");
		PersonalExit(1);
	}

	// open file
	fd = fopen(filename, "r");
	if (fd == NULL) {
		logging::print(LOG_ERROR "Can't read %s file. Aborting.\n", filename);
		PersonalExit(1);
	}
	free(filename);

	// read file until line with correct timestep
	do {
		ptr = fgets(buffer, sizeof(buffer), fd);
		sscanf(buffer, "%u", &line_timestep);
	} while ((line_timestep != timestep) && (ptr != NULL));

	if (ptr == NULL) {
		logging::print_master(LOG_ERROR "Can't read entry %u in 'planet%u.dat'!\n", timestep, get_planet_number());
		PersonalExit(1);
	}

	// close file
	fclose(fd);

	// move ptr until correct column
	while (column > 1) {
		ptr += strspn(ptr, "eE0123456789-.");
		ptr += strspn(ptr, "\t :=>_");
		column--;
	}

	sscanf(ptr, "%lf", &value);

	return value;
}
