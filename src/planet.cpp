#include "planet.h"
#include "LowTasks.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "output.h"
#include "util.h"
#include <cstdio>
#include <iostream>
#include <map>
#include <math.h>
#include <stdio.h>
#include <string.h>

// define the variables in the planet data file
const std::map<const std::string, const int> planet_file_column_v1 = {
    {"time step", 0},
    {"x", 1},
    {"y", 2},
    {"vx", 3},
    {"vy", 4},
    {"mass", 5},
    {"lost mass", 6},
    {"physical time", 7},
    {"omega frame", 8},
    {"mdcp", 9},
    {"exces mdcp", 10},
    {"eccentricity", 11},
    {"angular momentum", 12},
    {"semi-major axis", 13},
    {"omega", 14}};

// file version 2
const std::map<const std::string, const int> planet_file_column_v2 = {
    {"time step", 0},
    {"x", 1},
    {"y", 2},
    {"vx", 3},
    {"vy", 4},
    {"mass", 5},
    {"physical time", 6},
    {"omega frame", 7},
    {"mdcp", 8},
    {"exces mdcp", 9},
    {"eccentricity", 10},
    {"angular momentum", 11},
    {"semi-major axis", 12},
    {"omega", 13}};

// file version 2.1
const std::map<const std::string, const int> planet_file_column_v2_1 = {
    {"time step", 0},
    {"x", 1},
    {"y", 2},
    {"vx", 3},
    {"vy", 4},
    {"mass", 5},
    {"physical time", 6},
    {"omega frame", 7},
    {"mdcp", 8},
    {"exces mdcp", 9},
    {"eccentricity", 10},
    {"angular momentum", 11},
    {"semi-major axis", 12},
    {"omega kepler", 13},
    {"mean anomaly", 14},
    {"eccentric anomaly", 15},
    {"true anomaly", 16},
    {"pericenter angle", 17}};

auto planet_files_column = planet_file_column_v2_1;

const std::map<const std::string, const std::string> variable_units = {
    {"time step", "1"},
    {"x", "length"},
    {"y", "length"},
    {"vx", "velocity"},
    {"vy", "velocity"},
    {"mass", "mass"},
    {"lost mass", "mass"},
    {"physical time", "time"},
    {"omega frame", "frequency"},
    {"mdcp", "mass"},
    {"exces mdcp", "mass"},
    {"eccentricity", "1"},
    {"angular momentum", "angular_momentum"},
    {"semi-major axis", "length"},
    {"mean anomaly", "1"},
    {"eccentric anomaly", "1"},
    {"true anomaly", "1"},
    {"pericenter angle", "1"},
    {"omega", "frequency"},
    {"omega kepler", "frequency"}};



t_planet::~t_planet()
{
	delete[] m_name;
}

/**
	set name of planet
*/
void t_planet::set_name(const char *name)
{
    // delete old name
    delete[] m_name;

    // aquire space for new name
    m_name = new char[strlen(name) + 1];

    strcpy(m_name, name);
}

/**
 * @brief t_planet::get_angle get phi coordinate
 * @return
 */
double t_planet::get_phi() const { return atan2(m_y, m_x); }

/**
	get planet distance to coordinate center
*/
double t_planet::get_r() const { return sqrt(pow2(m_x) + pow2(m_y)); }

/**
	get ramp up mass of the planet
*/
double t_planet::get_rampup_mass()
{
    double ramping = 1.0;
    if (get_rampuptime() > 0) {
	if (PhysicalTime < get_rampuptime() * DT) {
	    ramping =
		1.0 -
		pow2(cos(PhysicalTime * PI / 2.0 / (get_rampuptime() * DT)));
	}
    }
    return get_mass() * ramping;
}

/**
	get planet period T
*/
double t_planet::get_period()
{
    return 2.0 * PI *
	   sqrt(pow3(get_semi_major_axis()) /
		((hydro_center_mass + get_mass()) * constants::G));
}

/**
	get omega_kepler at current planet location
*/
double t_planet::get_omega()
{
    double distance = get_r();
    if (!is_distance_zero(distance)) {
	return sqrt(((hydro_center_mass + get_mass()) * constants::G) / pow3(distance));
    } else {
	return 0.0;
    }
}

/**
	get hill radius at current planet location
*/
double t_planet::get_rhill()
{
    const double r = get_r();
	const double Mp = get_mass();
	const double Mstar = hydro_center_mass;
	const double rhill = pow(Mp/(3*Mstar), 1.0/3.0)*r;\
	return rhill;
}


/**
	get angular momentum of planet
*/
double t_planet::get_angular_momentum()
{
    // j = r x p = r x mv
    return get_mass() * get_x() * get_vy() - get_mass() * get_y() * get_vx();
}

void t_planet::create_planet_file()
{
    if (!CPU_Master)
	return;

    FILE *fd;
    char *filename = 0;

    std::string header_variable_description =
	output::text_file_variable_description(planet_files_column,
					       variable_units);

    // create normal file
    if (asprintf(&filename, "%splanet%u.dat", OUTPUTDIR, get_planet_number()) ==
	-1) {
	logging::print(LOG_ERROR "Not enough memory!\n");
	PersonalExit(1);
    }

    fd = fopen(filename, "w");
    free(filename);

    if (fd == NULL) {
	logging::print(LOG_ERROR "Can't write %s file. Aborting.\n", filename);
	PersonalExit(1);
    }
    fprintf(fd, "#FargoCPT planet file\n");
    fprintf(fd, "#version: 2\n");
    fprintf(fd, "%s", header_variable_description.c_str());
    fclose(fd);

    // create big file
    if (asprintf(&filename, "%sbigplanet%u.dat", OUTPUTDIR,
		 get_planet_number()) == -1) {
	logging::print(LOG_ERROR "Not enough memory!\n");
	PersonalExit(1);
    }

    fd = fopen(filename, "w");
    free(filename);

    if (fd == NULL) {
	logging::print(LOG_ERROR "Can't write %s file. Aborting.\n", filename);
	PersonalExit(1);
    }

    fprintf(fd, "#FargoCPT planet file\n");
    fprintf(fd, "#version: 2\n");
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
    if (asprintf(&filename, big_file ? "%sbigplanet%u.dat" : "%splanet%u.dat",
		 OUTPUTDIR, get_planet_number()) == -1) {
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

    fprintf(
	fd,
	"%d\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n",
	timestep, get_x(), get_y(), get_vx(), get_vy(), get_mass(),
	PhysicalTime, OmegaFrame, mdcp, exces_mdcp, get_eccentricity(),
	get_angular_momentum(), get_semi_major_axis(), get_omega(),
	get_mean_anomaly(), get_eccentric_anomaly(), get_true_anomaly(),
	get_pericenter_angle());

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

double t_planet::get_value_from_file(unsigned int timestep,
				     std::string variable_name)
{
    double value;
    int column = -1;

    std::string filename = std::string(OUTPUTDIR) + "planet" +
			   std::to_string(get_planet_number()) + ".dat";

    std::string version = output::get_version(filename);
    auto variable_columns = planet_file_column_v2;

    if (version == "1") {
	variable_columns = planet_file_column_v1;
    } else if (version == "2") {
	variable_columns = planet_file_column_v2;
    } else {
	std::cerr << "Unknown version '" << version << "'for planet.dat file!"
		  << std::endl;
	PersonalExit(1);
    }

    // check whether the column map contains the variable
    auto iter = variable_columns.find(variable_name);
    if (iter != variable_columns.end()) {
	column = iter->second;
    } else {
	std::cerr << "Unknown variable '" << variable_name
		  << "' for planet.dat file version '" << version << "'\n"
		  << std::endl;
    }

    if (column == -1) {
	std::cerr << "Something went wring in obtaining the column for "
		  << variable_name << " for the planet data file." << std::endl;
	PersonalExit(1);
    }
    value = output::get_from_ascii_file(filename, timestep, column);
    return value;
}

void t_planet::set_orbital_elements_zero()
{
    m_semi_major_axis = 0.0;
    m_eccentricity = 0.0;
    m_mean_anomaly = 0.0;
    m_true_anomaly = 0.0;
    m_eccentric_anomaly = 0.0;
    m_pericenter_angle = 0.0;
}

void t_planet::calculate_orbital_elements(double x, double y, double vx,
					  double vy, double com_mass)
{
    // mass of reference (primary for default star and sum of inner planet mass
    // otherwise)
    double Ax, Ay, e, d, h, a, E, M, V;
    double PerihelionPA;
    double temp;
    double m = com_mass + get_mass();

    h = x * vy - y * vx;
    d = sqrt(x * x + y * y);
    if (is_distance_zero(d) || h == 0.0) {
	set_orbital_elements_zero();
	return;
    }
    Ax = x * vy * vy - y * vx * vy - constants::G * m * x / d;
    Ay = y * vx * vx - x * vx * vy - constants::G * m * y / d;
    e = sqrt(Ax * Ax + Ay * Ay) / constants::G / m;
    a = h * h / constants::G / m / (1.0 - e * e);

    if (e != 0.0) {
	temp = (1.0 - d / a) / e;
	if (temp > 1.0) {
	    // E = acos(1)
	    E = 0.0;
	} else if (temp < -1.0) {
	    // E = acos(-1)
	    E = PI;
	} else {
	    E = acos(temp);
	}
    } else {
	E = 0.0;
    }

    if ((x * y * (vy * vy - vx * vx) + vx * vy * (x * x - y * y)) < 0) {
	E = -E;
    }

    M = E - e * sin(E);

    if (e != 0.0) {
	temp = (a * (1.0 - e * e) / d - 1.0) / e;
	if (temp > 1.0) {
	    // V = acos(1)
	    V = 0.0;
	} else if (temp < -1.0) {
	    // V = acos(-1)
	    V = PI;
	} else {
	    V = acos(temp);
	}
    } else {
	V = 0.0;
    }

    if (E < 0.0) {
	V = -V;
    }

    if (e != 0.0) {
	PerihelionPA = atan2(Ay, Ax);
    } else {
	PerihelionPA = atan2(y, x);
    }

    m_semi_major_axis = a;
    m_eccentricity = e;
    m_mean_anomaly = M;
    m_true_anomaly = V;
    m_eccentric_anomaly = V;
    m_pericenter_angle = PerihelionPA;
}
