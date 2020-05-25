#include "planetary_system.h"
#include "LowTasks.h"
#include "constants.h"
#include "fpe.h"
#include "global.h"
#include "logging.h"
#include "parameters.h"
#include "types.h"
#include <cstring>
#include <ctype.h>
#include <math.h>
#include <stdio.h>

#include <experimental/filesystem> /// To check wether we restart from a legacy simulation


extern boolean CICPlanet;
extern int Corotating;

t_planetary_system::t_planetary_system() {
	planet_restart_legacy = false;
}

t_planetary_system::~t_planetary_system()
{
    for (unsigned int i = 0; i < m_planets.size(); ++i) {
	delete m_planets.at(i);
    }

    m_planets.clear();
}

void t_planetary_system::init_rebound()
{
    m_rebound = reb_create_simulation();
    m_rebound->G = constants::G;
    m_rebound->dt = 1e-6;
    m_rebound->softening = 0.0; // 5e-4; // Jupiter radius in au
    m_rebound->integrator = reb_simulation::REB_INTEGRATOR_IAS15;
    // m_rebound->integrator = reb_simulation::REB_INTEGRATOR_MERCURIUS;
    // m_rebound->integrator = reb_simulation::REB_INTEGRATOR_WHFAST; // crashes
    for (unsigned int i = 0; i < get_number_of_planets(); ++i) {
	auto &planet = get_planet(i);
	struct reb_particle p;
	p.x = planet.get_x();
	p.y = planet.get_y();
	p.z = 0;
	p.vx = planet.get_vx();
	p.vy = planet.get_vy();
	p.vz = 0;
	p.ax = 0;
	p.ay = 0;
	p.az = 0;
	p.m = planet.get_mass();
	p.r = 0;
	p.lastcollision = 0;
	p.c = nullptr;
	p.hash = 0;
	p.ap = nullptr;
	p.sim = nullptr;
	reb_add(m_rebound, p);
    }
}

void t_planetary_system::initialize_default_star()
{
    t_planet *planet = new t_planet();
    initialize_planet_legacy(planet, M, 0.0, 0.0, 0.0);

    planet->set_name("Default Star");
    planet->set_acc(0.0);

    planet->set_radius(parameters::star_radius);
    planet->set_temperature(parameters::star_temperature);
    planet->set_irradiate(false);
    planet->set_rampuptime(0.0);

    planet->set_disk_on_planet_acceleration(Pair()); // initialize to zero
    planet->set_nbody_on_planet_acceleration(Pair());
    add_planet(planet);
}

void t_planetary_system::read_from_file(char *filename)
{
    FILE *fd;


	/// Note: techically we only need to do this if we are restarting the simulation. But this works anyway.
	/// Check if last planet file is missing
	/// If this is the case, we are restarting from a legacy simulation in which the central object did not have a planet file
   std::string last_filename = std::string(OUTPUTDIR) + "planet" +
			  std::to_string(get_number_of_planets()) + ".dat";
   planet_restart_legacy = !std::experimental::filesystem::exists(last_filename);


    if (parameters::default_star) {
	initialize_default_star();
    }

    // check if a filename was specified
    if (filename == NULL) {
	logging::print_master(LOG_INFO "No planetfile specified.\n");
    } else {
	// Only read from file if filename is given.

	// open fill
	fd = fopen(filename, "r");

	// check if file was readable
	if (fd == NULL) {
	    logging::print_master(LOG_ERROR "Error : can't find '%s'.\n",
				  filename);
	    PersonalExit(1);
	    return;
	}

	char buffer[512];

	// read line by line
	while (fgets(buffer, sizeof(buffer), fd) != NULL) {
		char name[256], feeldisk[8], feelother[8], irradiate[8];
	    double semi_major_axis, mass, acc, eccentricity = 0.0, temperature,
					       radius, phi, rampuptime;
	    int num_args;

	    // check if this line is a comment
	    if ((strlen(buffer) > 0) && (buffer[0] == '#'))
		continue;

	    // try to cut line into pieces
	    num_args = sscanf(
		buffer, "%255s %lf %lf %lf %7s %7s %lf %lf %lf %7s %lf %lf",
		name, &semi_major_axis, &mass, &acc, feeldisk, feelother,
		&eccentricity, &radius, &temperature, irradiate, &phi,
		&rampuptime);
	    if (num_args < 6)
		continue;

	    if (num_args < 7) {
		eccentricity = 0.0;
	    }

	    if (num_args < 8) {
		radius = 0.009304813;
	    }

	    if (num_args < 9) {
		temperature = 5778.0;
	    }

	    if (num_args < 10) {
		irradiate[0] = 'n';
	    }

	    if (num_args < 11) {
		phi = 0.0;
	    }

	    if (num_args < 12) {
		rampuptime = 0;
	    }

	    if (CICPlanet) {
		// initialization puts centered-in-cell planets (with
		// excentricity = 0 only)
		unsigned int j = 0;
		while (GlobalRmed[j] < semi_major_axis)
		    j++;
		semi_major_axis = Radii[j + 1];
	    }

	    t_planet *planet = new t_planet();

	    if (parameters::default_star) {
		initialize_planet_legacy(planet, mass, semi_major_axis,
					 eccentricity, phi);
	    } else {
		// planets starts at Periastron
		double nu = 0.0;

		double pericenter_angle = phi;
		if (get_number_of_planets() < 2) {
		    initialize_planet_jacobi_adjust_first_two(
			planet, mass, semi_major_axis, eccentricity,
			pericenter_angle, nu);
		} else {
		    initialize_planet_jacobi(planet, mass, semi_major_axis,
					     eccentricity, pericenter_angle,
					     nu);
		}
	    }
	    planet->set_name(name);
	    planet->set_acc(acc);

	    logging::print_master(
		LOG_WARNING,
		"Warning: feeldisk flag is deprecated. Interaction is now set globally by the DiskFeedback flag. Value is ignored!\n");
	    logging::print_master(
		LOG_WARNING,
		"Warning: feelother flag is deprecated. Interaction is now set globally by the DiskFeedback flag. Value is ignored!\n");

	    planet->set_radius(radius);
	    planet->set_temperature(temperature / units::temperature);
	    planet->set_irradiate(tolower(irradiate[0]) == 'y');
	    planet->set_rampuptime(rampuptime);

	    planet->set_disk_on_planet_acceleration(
		Pair()); // initialize to zero
	    planet->set_nbody_on_planet_acceleration(Pair());

	    add_planet(planet);
	}
	// close file
	fclose(fd);

	logging::print_master(LOG_INFO "%d planet(s) found.\n",
			      get_number_of_planets());
    }

    if (get_number_of_planets() > 0) {
	HillRadius =
	    get_planet(0).get_x() * pow(get_planet(0).get_mass() / 3., 1. / 3.);
    } else {
	HillRadius = 0;
    }

    // set up hydro frame center
    if (get_number_of_planets() == 0) {
	die("No stars or planets!");
    }
    if (parameters::n_bodies_for_hydroframe_center == 0) {
	// use all bodies to calculate hydro frame center
	parameters::n_bodies_for_hydroframe_center = get_number_of_planets();
    }
    if (parameters::n_bodies_for_hydroframe_center > get_number_of_planets()) {
	// use as many bodies to calculate hydro frame center as possible
	parameters::n_bodies_for_hydroframe_center = get_number_of_planets();
    }
    logging::print_master(
	LOG_INFO
	"The first %d planets are used to calculate the hydro frame center.\n",
	parameters::n_bodies_for_hydroframe_center);

    move_to_hydro_frame_center();

    if (Corotating == YES &&
	parameters::corotation_reference_body > get_number_of_planets() - 1) {
	die("Id of reference planet for corotation is not valid. Is '%d' but must be <= '%d'.",
	    parameters::corotation_reference_body, get_number_of_planets() - 1);
    }

    update_global_hydro_frame_center_mass();
    logging::print_master(
	LOG_INFO "The mass of the planets used as hydro frame center is %e.\n",
	hydro_center_mass);

	init_rebound();
}

void t_planetary_system::list_planets()
{
    if (!CPU_Master)
	return;

    calculate_orbital_elements();

    if (get_number_of_planets() == 0) {
	// logging::print(LOG_INFO "Planet overview: No planets specified.\n");
	return;
    }

    logging::print(LOG_INFO "Planet overview:\n");
    logging::print(LOG_INFO "\n");
    logging::print(
	LOG_INFO
	" #   | name                    | mass [m0]  | x [l0]     | y [l0]     | vx         | vy         |\n");
    logging::print(
	LOG_INFO
	"-----+-------------------------+------------+------------+------------+------------+------------+\n");

    for (unsigned int i = 0; i < get_number_of_planets(); ++i) {
	logging::print(
	    LOG_INFO
	    " %3i | %-23s | % 10.7g | % 10.7g | % 10.7g | % 10.7g | % 10.7g |\n",
	    i, get_planet(i).get_name(), get_planet(i).get_mass(),
	    get_planet(i).get_x(), get_planet(i).get_y(),
	    get_planet(i).get_vx(), get_planet(i).get_vy());
    }

    logging::print(LOG_INFO "\n");
    logging::print(
	LOG_INFO
	" #   | e          | a          | T [t0]     | T [a]      | accreting  | feels disk | feels plan.|\n");
    logging::print(
	LOG_INFO
	"-----+------------+------------+------------+------------+------------+------------+------------+\n");

    for (unsigned int i = 0; i < get_number_of_planets(); ++i) {
	logging::print(
	    LOG_INFO
	    " %3i | % 10.7g | % 10.7g | % 10.7g | % 10.6g | % 10.7g |          %c |          %c |\n",
	    i, get_planet(i).get_eccentricity(),
	    get_planet(i).get_semi_major_axis(), get_planet(i).get_period(),
	    get_planet(i).get_period() * units::time.get_cgs_factor() /
		(24 * 60 * 60 * 365.2425),
		get_planet(i).get_acc(), '-', '-');
    }

    logging::print(LOG_INFO "\n");
    logging::print(
	LOG_INFO
	" #   | Temp [K]   | R [l0]     | irradiates | rampuptime |\n");
    logging::print(
	LOG_INFO
	"-----+------------+------------+------------+------------+\n");

    for (unsigned int i = 0; i < get_number_of_planets(); ++i) {
	logging::print(LOG_INFO
		       " %3i | % 10.7g | % 10.7g |          %c | % 10.7g |\n",
		       i, get_planet(i).get_temperature() * units::temperature,
		       get_planet(i).get_radius(),
		       (get_planet(i).get_irradiate()) ? 'X' : '-',
		       get_planet(i).get_rampuptime());
    }

    logging::print(LOG_INFO "\n");
}

void t_planetary_system::rotate(double angle)
{
    for (unsigned int i = 0; i < get_number_of_planets(); ++i) {
	// rotate positions
	double old_x = get_planet(i).get_x();
	double old_y = get_planet(i).get_y();
	get_planet(i).set_x(old_x * cos(angle) + old_y * sin(angle));
	get_planet(i).set_y(-old_x * sin(angle) + old_y * cos(angle));

	// rotate velocities
	double old_vx = get_planet(i).get_vx();
	double old_vy = get_planet(i).get_vy();
	get_planet(i).set_vx(old_vx * cos(angle) + old_vy * sin(angle));
	get_planet(i).set_vy(-old_vx * sin(angle) + old_vy * cos(angle));
    }
}

void t_planetary_system::restart(unsigned int timestep)
{
    for (unsigned int i = 0; i < get_number_of_planets(); ++i) {
		if(planet_restart_legacy && i == 0) { /// Legacy restart means that the central star has no file, thus we skip restarting it
			continue;
		}
	get_planet(i).restart(timestep);
    }
	m_rebound->t = PhysicalTime;
}

void t_planetary_system::create_planet_files()
{
    for (unsigned int i = 0; i < get_number_of_planets(); ++i) {
	get_planet(i).create_planet_file();
    }
}

void t_planetary_system::write_planets(unsigned int timestep, bool big_file)
{
    for (unsigned int i = 0; i < get_number_of_planets(); ++i) {
		if(planet_restart_legacy && i == 0) { /// Legacy restart means that the central star has no file, thus we skip writing it
			continue;
		}
	get_planet(i).write(timestep, big_file);
    }
}

/**
   Initialize the planets position and velocity in the legacy way
*/
void t_planetary_system::initialize_planet_legacy(t_planet *planet, double mass,
						  double semi_major_axis,
						  double eccentricity,
						  double phi)
{
    planet->set_mass(mass);
    // planets starts at Apastron
    double r = semi_major_axis * (1.0 + eccentricity);
    planet->set_x(r * cos(phi));
    planet->set_y(r * sin(phi));
    double v = 0.0;
    if (semi_major_axis != 0.0) {
	v = sqrt(constants::G * (1.0 + mass) / semi_major_axis) *
	    sqrt((1.0 - eccentricity) / (1.0 + eccentricity));
    }
    planet->set_vx(-v * sin(phi));
    planet->set_vy(v * cos(phi));
}

/**
   Set the position and velocity of the first planet according
   to the second planets position and velocity.

   This is done because the orbital elements of a single object are meaningless
   and the first two bodies have to be initialized together.
*/
void t_planetary_system::initialize_planet_jacobi_adjust_first_two(
    t_planet *planet, double mass, double semi_major_axis, double eccentricity,
    double omega, double true_anomaly)
{
    if (get_number_of_planets() == 0) {
	// first planet always goes to origin
	planet->set_mass(mass);
	planet->set_x(0.0);
	planet->set_y(0.0);
	planet->set_vx(0.0);
	planet->set_vy(0.0);
    } else {
	// initialize the second planet around the origin, such that the two
	// bodies have the correct separation
	initialize_planet_jacobi(planet, mass, semi_major_axis, eccentricity,
				 omega, true_anomaly);
	t_planet *planet1 = m_planets[0];
	t_planet *planet2 = planet;

	double m1 = planet1->get_mass();
	double m2 = planet2->get_mass();

	// values are now such that the distance between the first two bodies
	// match

	double x = planet2->get_x();
	double y = planet2->get_y();
	double vx = planet2->get_vx();
	double vy = planet2->get_vy();

	// move both bodies into barycenter

	double k1 = m2 / (m1 + m2);
	planet1->set_x(-k1 * x);
	planet1->set_y(-k1 * y);
	planet1->set_vx(-k1 * vx);
	planet1->set_vy(-k1 * vy);

	double k2 = m1 / (m1 + m2);
	planet2->set_x(k2 * x);
	planet2->set_y(k2 * y);
	planet2->set_vx(k2 * vx);
	planet2->set_vy(k2 * vy);
    }
}

/**
   Initialize the planets position and velocity using jacobian coordinates
*/
void t_planetary_system::initialize_planet_jacobi(t_planet *planet, double mass,
						  double semi_major_axis,
						  double eccentricity,
						  double omega,
						  double true_anomaly)
{
    planet->set_mass(mass);
    Pair com = get_center_of_mass(); // of all previously added planets
    double com_mass = get_mass();    // of all previously added planets

    // some temporary variables for optimization and legibility
    double cos_ota = cos(omega + true_anomaly);
    double sin_ota = sin(omega + true_anomaly);
    double cos_o = cos(omega);
    double sin_o = sin(omega);
    double cos_ta = cos(true_anomaly);
    double sin_ta = sin(true_anomaly);

    double r = semi_major_axis * (1 - eccentricity * eccentricity) /
	       (1 + eccentricity * cos_ta);
    double x = com.x + r * cos_ota;
    double y = com.y + r * sin_ota;

    double v = 0.0;
    if (semi_major_axis > 0.0) {
	v = sqrt(constants::G * (com_mass + mass) /
		 (semi_major_axis * (1 - eccentricity * eccentricity)));
    }

    double vx = v * (-cos_o * sin_ta - sin_o * (eccentricity + cos_ta));
    double vy = v * (-sin_o * sin_ta + cos_o * (eccentricity + cos_ta));

    planet->set_x(x);
    planet->set_y(y);
    planet->set_vx(vx);
    planet->set_vy(vy);
}

/**
   Get the sum of masses of the first n particles
*/
double t_planetary_system::get_mass(unsigned int n)
{
    double mass = 0.0;
    for (unsigned int i = 0; i < n; i++) {
	mass += get_planet(i).get_mass();
    }
    return mass;
}

/**
   Get the sum of masses of all particles
*/
double t_planetary_system::get_mass()
{
    return get_mass(get_number_of_planets());
}

/**
   Get the center of mass of the first n particles
*/
Pair t_planetary_system::get_center_of_mass(unsigned int n)
{
    double x = 0.0;
    double y = 0.0;
    double mass = 0.0;
    for (unsigned int i = 0; i < n; i++) {
	t_planet &planet = get_planet(i);
	mass += planet.get_mass();
	x += planet.get_x() * planet.get_mass();
	y += planet.get_y() * planet.get_mass();
    }
    Pair com;
    if (mass > 0) {
	com.x = x / mass;
	com.y = y / mass;
    } else {
	com.x = 0.0;
	com.y = 0.0;
    }
    return com;
}

/**
   Get the velocity of the center of mass of the first n particles
*/
Pair t_planetary_system::get_center_of_mass_velocity(unsigned int n)
{
    double vx = 0.0;
    double vy = 0.0;
    double mass = 0.0;
    for (unsigned int i = 0; i < n; i++) {
	t_planet &planet = get_planet(i);
	mass += planet.get_mass();
	vx += planet.get_vx() * planet.get_mass();
	vy += planet.get_vy() * planet.get_mass();
    }
    Pair vcom;
    if (mass > 0) {
	vcom.x = vx / mass;
	vcom.y = vy / mass;
    } else {
	vcom.x = 0.0;
	vcom.y = 0.0;
    }
    return vcom;
}

/**
   Get the center of mass of all particles
*/
Pair t_planetary_system::get_center_of_mass()
{
    return get_center_of_mass(get_number_of_planets());
}

/**
   Get the center of coordinate system as chosen by
   parameters::n_bodies_for_hydroframe_center.
   This is the function to be called later in the code
   whenever the distance to the center is used.
*/
Pair t_planetary_system::get_hydro_frame_center_position()
{
    return get_center_of_mass(parameters::n_bodies_for_hydroframe_center);
}

/**
   Analogous to get_hydro_frame_center but returns its velocity.
*/

Pair t_planetary_system::get_hydro_frame_center_velocity()
{
    return get_center_of_mass_velocity(
	parameters::n_bodies_for_hydroframe_center);
}

/**
   Analogous to get_hydro_frame_center but returns its mass.
*/

double t_planetary_system::get_hydro_frame_center_mass()
{
    return get_mass(parameters::n_bodies_for_hydroframe_center);
}

/**
   Update the global variable hydro_center_mass.
*/

void t_planetary_system::update_global_hydro_frame_center_mass()
{
    hydro_center_mass = get_hydro_frame_center_mass();
}

/**
   Move the planetary system to the chosen frame center.
 */
void t_planetary_system::move_to_hydro_frame_center()
{
    Pair center = get_hydro_frame_center_position();
    Pair vcenter = get_hydro_frame_center_velocity();
    for (unsigned int i = 0; i < get_number_of_planets(); i++) {
	t_planet &planet = get_planet(i);
	double x = planet.get_x();
	double y = planet.get_y();
	double vx = planet.get_vx();
	double vy = planet.get_vy();

	planet.set_x(x - center.x);
	planet.set_y(y - center.y);
	planet.set_vx(vx - vcenter.x);
	planet.set_vy(vy - vcenter.y);
    }
}

/**
   Calculate orbital elements of all planets.
 */
void t_planetary_system::calculate_orbital_elements()
{
    double x, y, vx, vy, M;
    for (unsigned int i = 0; i < get_number_of_planets(); i++) {
	auto &planet = get_planet(i);
	if (i == 0) {
		get_planet(0).set_orbital_elements_zero();
		continue;
	}

	Pair com_pos = get_center_of_mass(i);
	Pair com_vel = get_center_of_mass_velocity(i);

	M = get_mass(i);
	x = planet.get_x() - com_pos.x;
	y = planet.get_y() - com_pos.y;
	vx = planet.get_vx() - com_vel.x;
	vy = planet.get_vy() - com_vel.y;

	planet.calculate_orbital_elements(x, y, vx, vy, M);
    }
}


/**
   Copy positions, velocities and masses
   from planetary system to rebound.
*/
void t_planetary_system::copy_data_to_rebound() {
	for (unsigned int i = 0; i < get_number_of_planets(); i++) {
	    auto &planet = get_planet(i);
	    m_rebound->particles[i].x = planet.get_x();
	    m_rebound->particles[i].y = planet.get_y();
	    m_rebound->particles[i].vx = planet.get_vx();
	    m_rebound->particles[i].vy = planet.get_vy();
	    m_rebound->particles[i].m = planet.get_mass();
	}
}

/**
   Copy positions, velocities and masses back
   from rebound to planetary system.
*/
void t_planetary_system::copy_data_from_rebound() {
		for (unsigned int i = 0; i < get_number_of_planets(); i++) {
	    auto &planet = get_planet(i);
	    planet.set_x(m_rebound->particles[i].x);
	    planet.set_y(m_rebound->particles[i].y);
	    planet.set_vx(m_rebound->particles[i].vx);
	    planet.set_vy(m_rebound->particles[i].vy);
	}
}

/**
   Integrate the nbody system forward in time using rebound.
*/
void t_planetary_system::integrate(double time, double dt) {
    if (get_number_of_planets() < 2) {
	// don't integrate a single particle that doesn't move
		return;
    }

	copy_data_to_rebound();
	m_rebound->t = time;

	disable_trap_fpe_gnu();
	reb_integrate(m_rebound, time + dt);
	enable_trap_fpe_gnu();

	copy_data_from_rebound();

	move_to_hydro_frame_center();
}


/**
	Updates planets velocities due to disk influence if "DiskFeedback" is
   set.
*/
void t_planetary_system::correct_velocity_for_disk_accel()
{

	if (!parameters::disk_feedback){
		return;
	}

	for (unsigned int k = 0; k < get_number_of_planets(); k++) {

		/*
		 * from centrifugal balance folows
		 * v_new**2 / r = a_disk + v_old**2 / r
		 * v_new = sqrt(v_old**2 - r * a_disk)
		 */
		t_planet &planet = get_planet(k);

		const Pair gas_accel = planet.get_disk_on_planet_acceleration();
		const double vx_old = planet.get_vx();
		const double vy_old = planet.get_vy();
		const double v_old = std::sqrt(std::pow(vx_old, 2.0) + std::pow(vy_old, 2.0));

		if(v_old == 0.0)
		{
			continue;
		}

		const double x = planet.get_x();
		const double y = planet.get_y();
		const double specific_torque_gas = gas_accel.x * x + gas_accel.y * y; // = a_disk * r

		if(specific_torque_gas > std::pow(v_old, 2.0))
		{
			continue;
		}

		const double v_new = std::sqrt(std::pow(v_old, 2.0) - specific_torque_gas);

		const double new_vx = v_new / v_old * vx_old;
		const double new_vy = v_new / v_old * vy_old;

		planet.set_vx(new_vx);
		planet.set_vy(new_vy);
	}
}
