#include "planet.h"
#include "LowTasks.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "output.h"
#include "util.h"
#include <cstdio>
#include <iostream>
#include <fstream>
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

// file version 2.2
const std::map<const std::string, const int> planet_file_column_v2_2 = {
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
    {"pericenter angle", 17},
    {"torque", 18}};

// file version 2.2
const std::map<const std::string, const int> planet_file_column_v2_3 = {
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
	{"pericenter angle", 17},
	{"torque", 18},
	{"accreted mass", 19}};

auto planet_files_column = planet_file_column_v2_3;

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
    {"omega kepler", "frequency"},
	{"torque", "torque"},
	{"accreted mass", "mass"}};

t_planet::~t_planet() { }

t_planet::t_planet() {
	 m_mass = 0.0;
	 m_x = 0.0;
	 m_y = 0.0;
	 m_vx = 0.0;
	 m_vy = 0.0;

	 m_acc = 0.0;
	 m_accreted_mass = 0.0;
	 m_name = "";

	 m_planet_number = 0;
	 m_temperature = 0.0;
	 m_radius = 0.0;
	 m_irradiate = false;
	 m_rampuptime = 0.0;
	 m_disk_on_planet_acceleration = {0.0, 0.0};
	 m_nbody_on_planet_acceleration = {0.0, 0.0};
	 m_semi_major_axis = 0.0;
	 m_eccentricity = 0.0;
	 m_mean_anomaly = 0.0;
	 m_true_anomaly = 0.0;
	 m_eccentric_anomaly = 0.0;
	 m_pericenter_angle = 0.0;
	 m_torque = 0.0;
}

/**
	set name of planet
*/
void t_planet::set_name(std::string name)
{
	m_name = name;
}

/**
 * @brief t_planet::get_angle get phi coordinate
 * @return
 */
double t_planet::get_phi() const { return std::atan2(m_y, m_x); }

/**
	get planet distance to coordinate center
*/
double t_planet::get_r() const
{
    return std::sqrt(std::pow(m_x, 2) + std::pow(m_y, 2));
}

/**
	get ramp up mass of the planet
*/
double t_planet::get_rampup_mass()
{
    double ramping = 1.0;
    if (get_rampuptime() > 0) {
	if (PhysicalTime < get_rampuptime() * get_period()) {
		ramping = 1.0 - std::pow(std::cos(PhysicalTime * M_PI_2 /
						  (get_rampuptime() * get_period())),
				     2);
	}
    }
    return get_mass() * ramping;
}

/**
	get planet period T
*/
double t_planet::get_period() const
{
    return 2.0 * M_PI *
	   std::sqrt(std::pow(get_semi_major_axis(), 3) /
		     ((hydro_center_mass + get_mass()) * constants::G));
}

/**
	get omega_kepler at current planet location
*/
double t_planet::get_omega() const
{
    double distance = get_r();
    if (!is_distance_zero(distance)) {
	return std::sqrt(((hydro_center_mass + get_mass()) * constants::G) /
			 std::pow(distance, 3));
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
    const double rhill = std::pow(Mp / (3 * Mstar), 1.0 / 3.0) * r;
    return rhill;
}

/**
	get angular momentum of planet
*/
double t_planet::get_angular_momentum() const
{
    // j = r x p = r x mv
    return get_mass() * get_x() * get_vy() - get_mass() * get_y() * get_vx();
}

void t_planet::copy(const planet_member_variables &other){
	 m_mass = other.m_mass;
	 m_x = other.m_x;
	 m_y = other.m_y;
	 m_vx = other.m_vx;
	 m_vy = other.m_vy;

	 m_acc = other.m_acc;
	 m_accreted_mass = other.m_accreted_mass;

	m_planet_number = other.m_planet_number;
	m_temperature = other.m_temperature;
	m_radius = other.m_radius;
	m_irradiate = other.m_irradiate;
	m_rampuptime = other.m_rampuptime;
	m_disk_on_planet_acceleration = other.m_disk_on_planet_acceleration;
	m_nbody_on_planet_acceleration = other.m_nbody_on_planet_acceleration;

	/// orbital elements
	 m_semi_major_axis = other.m_semi_major_axis;
	 m_eccentricity = other.m_eccentricity;
	 m_mean_anomaly = other.m_mean_anomaly;
	 m_true_anomaly = other.m_true_anomaly;
	 m_eccentric_anomaly = other.m_eccentric_anomaly;
	 m_pericenter_angle = other.m_pericenter_angle;

	 m_torque = other.m_torque;

}

void t_planet::create_planet_file(bool debug_output)
{
    if (!CPU_Master)
	return;

    FILE *fd;
    char *filename = 0;

    std::string header_variable_description =
	output::text_file_variable_description(planet_files_column,
					       variable_units);

	if(debug_output){
		// create normal file
		if (asprintf(&filename, "%sdebugplanet%u.bin", OUTPUTDIR, get_planet_number()) ==
		-1) {
		logging::print(LOG_ERROR "Not enough memory!\n");
		PersonalExit(1);
		}
	} else {
		// create normal file
		if (asprintf(&filename, "%splanet%u.bin", OUTPUTDIR, get_planet_number()) ==
		-1) {
		logging::print(LOG_ERROR "Not enough memory!\n");
		PersonalExit(1);
		}
	}

    // create big file
    if (asprintf(&filename, "%sbigplanet%u.dat", OUTPUTDIR,
		 get_planet_number()) == -1) {
	logging::print(LOG_ERROR "Not enough memory!\n");
	PersonalExit(1);
    }

    fd = fopen(filename, "w");

    if (fd == NULL) {
	logging::print(LOG_ERROR "Can't write %s file. Aborting.\n", filename);
	PersonalExit(1);
    }

    fprintf(fd, "#FargoCPT planet file\n");
    fprintf(fd, "#version: 2\n");
    fprintf(fd, "%s", header_variable_description.c_str());

	free(filename);
    fclose(fd);
}

void t_planet::write(const unsigned int timestep, const unsigned int file_type)
{
    if (!CPU_Master)
	return;

    char *filename = 0;

    // create filename
	switch (file_type){
		case 0:
			if(asprintf(&filename, "%splanet%u.bin",
					 OUTPUTDIR, get_planet_number()) == -1) {
				logging::print(LOG_ERROR "Not enough memory!\n");
				PersonalExit(1);
			}
			write_binary(filename, timestep);
			break;
		case 1:
			if(asprintf(&filename, "%sbigplanet%u.dat",
					 OUTPUTDIR, get_planet_number()) == -1) {
				logging::print(LOG_ERROR "Not enough memory!\n");
				PersonalExit(1);
			}
			write_ascii(filename, timestep);
			// reset accreted mass, total accreted mass should be computed in postprocessing
			// don't forget to only use unique values to not overestimate mass accreted.
			// otherwise restarting simulations causes confusion as accreted mass would start from 0 again
			m_accreted_mass = 0.0;
			break;
		case 2:
			if(asprintf(&filename, "%sdebugplanet%u.bin",
					 OUTPUTDIR, get_planet_number()) == -1) {
				logging::print(LOG_ERROR "Not enough memory!\n");
				PersonalExit(1);
			}
			write_binary(filename, timestep);
			break;
		default:
			die("Bad file_type value for writing planet files!\n");
	}

	free(filename);
}

void t_planet::write_ascii(const char *filename, const unsigned int timestep) const
{
	// open file
	FILE *fd = fopen(filename, "a");
	if (fd == NULL) {
	logging::print(LOG_ERROR "Can't write %s file. Aborting.\n", filename);
	PersonalExit(1);
	}

	fprintf(
	fd,
	"%d\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\t%#.18g\n",
	timestep, get_x(), get_y(), get_vx(), get_vy(), get_mass(),
	PhysicalTime, OmegaFrame, mdcp, exces_mdcp, get_eccentricity(),
	get_angular_momentum(), get_semi_major_axis(), get_omega(),
	get_mean_anomaly(), get_eccentric_anomaly(), get_true_anomaly(),
	get_pericenter_angle(), get_torque(), get_accreted_mass());

	// close file
	fclose(fd);
}

void t_planet::write_binary(const char *filename, const unsigned int timestep) const
{

	std::ofstream wf(filename, std::ios::out | std::ios::binary | std::ios::app);
	if(!wf) {
	   logging::print(LOG_ERROR "Can't write %s file. Aborting.\n", filename);
	   die("End\n");
	}

	planet_member_variables pl;
	pl.timestep = timestep;
	pl.m_mass = m_mass;
	pl.m_x = m_x;
	pl.m_y = m_y;
	pl.m_vx = m_vx;
	pl.m_vy = m_vy;

	pl.m_acc = m_acc;
	pl.m_accreted_mass = m_accreted_mass;
	pl.m_planet_number = m_planet_number;
	pl.m_temperature = m_temperature;
	pl.m_radius = m_radius;
	pl.m_irradiate = m_irradiate;
	pl.m_rampuptime = m_rampuptime;
	pl.m_disk_on_planet_acceleration = m_disk_on_planet_acceleration;
	pl.m_nbody_on_planet_acceleration = m_nbody_on_planet_acceleration;

	/// orbital elements
	pl.m_semi_major_axis = m_semi_major_axis;
	pl.m_eccentricity = m_eccentricity;
	pl.m_mean_anomaly = m_mean_anomaly;
	pl.m_true_anomaly = m_true_anomaly;
	pl.m_eccentric_anomaly = m_eccentric_anomaly;
	pl.m_pericenter_angle = m_pericenter_angle;

	pl.m_torque = m_torque;

	wf.write((char*) (&pl), sizeof(planet_member_variables));
	wf.close();
}

void t_planet::restart(unsigned int timestep, bool debug)
{

	std::string filename;
	if(debug){
		filename = std::string(OUTPUTDIR) + "debugplanet" +
			   std::to_string(get_planet_number()) + ".bin";
	} else {
		filename = std::string(OUTPUTDIR) + "planet" +
			   std::to_string(get_planet_number()) + ".bin";
	}

	std::ifstream rf(filename, std::ofstream::binary | std::ios::in);

	if (!rf.is_open()) {
		logging::print_master(LOG_ERROR
				  "Can't read '%s' file. Aborting.\n", filename.c_str());
		PersonalExit(1);
	}

	//rf.ignore(timestep * sizeof(planet_member_variables));

	planet_member_variables pl;
	rf.read((char *) &pl, sizeof(planet_member_variables));

	while(pl.timestep != timestep)
	{
		printf("pl.timestep = %d	%d\n", pl.timestep, timestep);
		if(rf.eof())
		{
			logging::print_master(LOG_ERROR "Could not read timestep %d in %s\n", timestep, filename.c_str());
			die("End\n");
		}
		rf.read((char *) &pl, sizeof(planet_member_variables));
	}

	copy(pl);
	rf.close();
}

double t_planet::get_value_from_file(unsigned int timestep,
					 std::string variable_name, bool debug)
{
    double value;
    int column = -1;

	std::string filename;
	if(debug){
		filename = std::string(OUTPUTDIR) + "debugplanet" +
			   std::to_string(get_planet_number()) + ".dat";
	} else {
		filename = std::string(OUTPUTDIR) + "planet" +
			   std::to_string(get_planet_number()) + ".dat";
	}

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
	value = output::get_from_ascii_file(filename, timestep, column, debug);
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

	m_torque = 0.0;
	m_accreted_mass = 0.0;
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
    d = std::sqrt(x * x + y * y);
    if (is_distance_zero(d) || h == 0.0) {
	set_orbital_elements_zero();
	return;
    }
    Ax = x * vy * vy - y * vx * vy - constants::G * m * x / d;
    Ay = y * vx * vx - x * vx * vy - constants::G * m * y / d;
    e = std::sqrt(Ax * Ax + Ay * Ay) / constants::G / m;
    a = h * h / constants::G / m / (1.0 - e * e);

    if (e != 0.0) {
	temp = (1.0 - d / a) / e;
	if (temp > 1.0) {
	    // E = acos(1)
	    E = 0.0;
	} else if (temp < -1.0) {
	    // E = acos(-1)
	    E = M_PI;
	} else {
	    E = std::acos(temp);
	}
    } else {
	E = 0.0;
    }

    if ((x * y * (vy * vy - vx * vx) + vx * vy * (x * x - y * y)) < 0) {
	E = -E;
    }

    M = E - e * std::sin(E);

    if (e != 0.0) {
	temp = (a * (1.0 - e * e) / d - 1.0) / e;
	if (temp > 1.0) {
	    // V = acos(1)
	    V = 0.0;
	} else if (temp < -1.0) {
	    // V = acos(-1)
	    V = M_PI;
	} else {
	    V = std::acos(temp);
	}
    } else {
	V = 0.0;
    }

    if (E < 0.0) {
	V = -V;
    }

    if (e != 0.0) {
	PerihelionPA = std::atan2(Ay, Ax);
    } else {
	PerihelionPA = std::atan2(y, x);
    }

    m_semi_major_axis = a;
    m_eccentricity = e;
    m_mean_anomaly = M;
    m_true_anomaly = V;
    m_eccentric_anomaly = V;
    m_pericenter_angle = PerihelionPA;
}
