/**
	\file particles.cpp
	\author Tobias Mueller <Tobias_Mueller@twam.info>
*/

#include "particles.h"
#include "Force.h"
#include "LowTasks.h"
#include "Pframeforce.h"
#include "SourceEuler.h"
#include "commbound.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "parameters.h"
#include "selfgravity.h"
#include "util.h"
#include <math.h>
#include <mpi.h>
#include <random>
#include <stdlib.h>
#include <vector>
#include "find_cell_id.h"

extern Pair IndirectTerm;

namespace particles
{

/// local particle storage
t_particle *particles;

/// current size of local particle storage
unsigned int particles_size;

/// current number of particles on this node
unsigned int local_number_of_particles;

/// global number of particles;
unsigned int global_number_of_particles;

double local_r_min;
double local_r_max;

static MPI_Datatype mpi_particle;

// inverse of the Cumulative distribution function for a slope proportional to
// r^n
static double f(double x, double n)
{
    double rmin = parameters::particle_minimum_radius;
    double rmax = parameters::particle_maximum_radius;

    // make sure particles are not initialized on an orbit that moves out of
    // particle-bounds
    if (RMAX < parameters::particle_maximum_radius *
		   (1.0 + parameters::particle_eccentricity))
	rmax = parameters::particle_maximum_radius /
	       (1.0 + parameters::particle_eccentricity);

    if (RMIN > parameters::particle_minimum_radius *
		   (1.0 - parameters::particle_eccentricity))
	rmin = parameters::particle_minimum_radius /
	       (1.0 - parameters::particle_eccentricity);

    if (n == -1) {
	return exp(log(rmax / rmin) * x) * rmin;
    } else {
	return exp(
	    log(x * (-pow(rmin, n + 1) + pow(rmax, n + 1)) + pow(rmin, n + 1)) /
	    (n + 1));
    }
}

static double check_angle(double &phi)
{
    if (phi >= 2.0 * PI) {
	return phi -= 2.0 * PI;
    } else {
	if (phi < 0.0) {
	    return phi += 2.0 * PI;
	} else {
	    return phi;
	}
    }
}

static void transformCartCyl(double *cart, double *cyl, double r, double phi)
{
    (void)r;
    cyl[0] = cart[0] * cos(phi);
    cyl[0] += cart[1] * sin(phi);

    cyl[1] = -cart[0] * sin(phi);
    cyl[1] += cart[1] * cos(phi);
}

void find_nearest(unsigned int &n_radial_a_minus,
		  unsigned int &n_radial_a_plus,
		  unsigned int &n_azimuthal_a_minus,
		  unsigned int &n_azimuthal_a_plus,
		  unsigned int &n_radial_b_minus, unsigned int &n_radial_b_plus,
		  unsigned int &n_azimuthal_b_minus,
		  unsigned int &n_azimuthal_b_plus, const double r,
		  const double phi)
{
	n_radial_a_minus = get_rinf_id(parameters::radial_grid_type, r);
	n_radial_a_plus = n_radial_a_minus + 1;

	n_radial_b_minus = get_rmed_id(parameters::radial_grid_type, r);
	n_radial_b_plus = n_radial_b_minus + 1;

	n_azimuthal_a_minus = clamp_phi_id_to_grid(get_inf_azimuthal_id(phi));
	n_azimuthal_a_plus = get_next_azimuthal_id(n_azimuthal_a_minus);

	n_azimuthal_b_minus = clamp_phi_id_to_grid(get_med_azimuthal_id(phi));
	n_azimuthal_b_plus = get_next_azimuthal_id(n_azimuthal_b_minus);


}

static double interpolate_bilinear_sg(const double *array1D,
				      const unsigned int n_radial_minus,
				      const unsigned int n_radial_plus,
				      const unsigned int n_azimuthal_minus,
				      const unsigned int n_azimuthal_plus,
				      const double r, double phi)
{

    const unsigned int last_index = NAzimuthal - 1;
    double dphi = 2.0 * PI / (double)NAzimuthal;

    // values at corners
    double Qmm = array1D[n_radial_minus * NAzimuthal + n_azimuthal_minus];
    double Qpm = array1D[n_radial_plus * NAzimuthal + n_azimuthal_minus];
    double Qmp = array1D[n_radial_minus * NAzimuthal + n_azimuthal_plus];
    double Qpp = array1D[n_radial_plus * NAzimuthal + n_azimuthal_plus];

    double rm = Rb[n_radial_minus];
    double rp = Rb[n_radial_plus];

    double phim;
    double phip;

    if (n_azimuthal_minus == last_index) {
	if (phi < PI) {
	    // particle is inside cell with n_azimuthal = 0
	    // previous cell is at N_azimuthal_max - 1, but measure distance
	    // from -0.5*dphi
	    phim = -0.5 * dphi;
	    phip = 0.5 * dphi;
	} else // phi > PI
	{
	    // particle is inside cell with n_azimuthal = N_azimuthal_max - 1
	    // next cell is at N_azimuthal = 0, but measure distance from
	    // 2*PI+0.5*dphi
	    phim = (n_azimuthal_minus + 0.5) * dphi;
	    phip = (n_azimuthal_minus + 1.5) * dphi;
	}
    } else {
	phim = (n_azimuthal_minus + 0.5) * dphi;
	phip = (n_azimuthal_plus + 0.5) * dphi;
    }

    if ((phim > 2.0 * PI) || (phip > 2.0 * PI)) {
	phim -= PI;
	phip -= PI;
	phi -= PI;
    }

    double Qm = ((rm * phip - rm * phi) * Qmm + (rm * phi - rm * phim) * Qmp) /
		(rm * phip - rm * phim);
    double Qp = ((rp * phip - rp * phi) * Qpm + (rp * phi - rp * phim) * Qpp) /
		(rp * phip - rp * phim);

    return ((rp - r) * Qm + (r - rm) * Qp) / (rp - rm);
}

/*
static double get_particle_eccentricity(int particle_index)
{
	// used for debugging the integrators
	const int i = particle_index;

	// Runge-Lenz vector A = (p x L) - m * G * m * hydro_center_mass * r/|r|
	const double kappa = constants::G * hydro_center_mass *
particles[i].mass; const double reduced_mass =
hydro_center_mass*particles[i].mass / (hydro_center_mass + particles[i].mass);
	const double A_r = reduced_mass * pow3(particles[i].r) *
pow2(particles[i].phi_dot) - kappa; const double A_phi = reduced_mass *
pow2(particles[i].r) * particles[i].r_dot * particles[i].phi_dot; const double A
= sqrt(A_r*A_r + A_phi*A_phi); const double eccentricity = A / kappa;

	return eccentricity;
}
*/

/**
 * @brief init_particle_timestep init particle timestep for adaptive explicit
 * integrator
 * @param data
 */
static void init_particle_timestep(t_data &data)
{

    for (unsigned int i = 0; i < local_number_of_particles; ++i) {
	double sk, h1, der2, der12, sqr;

	int iord = 5;
	double dnf = 0.0;
	double dny = 0.0;
	double atoli = 1e-12;
	double rtoli = 0.0;
	particles[i].facold = 1.0e-4;

	double temp_ar, temp_aphi, temp_r, temp_phi, temp_r_dot, temp_phi_dot;

	temp_r = particles[i].r;
	temp_phi = particles[i].phi;
	temp_r_dot = particles[i].r_dot;
	temp_phi_dot = particles[i].phi_dot;

	// Cash–Karp method
	// (http://en.wikipedia.org/wiki/Cash%E2%80%93Karp_method) cartesian
	// coordinates are written inside the polar coordinates
	if (CartesianParticles) {
	    calculate_accelerations_from_star_and_planets_cart(
		temp_ar, temp_aphi, temp_r, temp_phi, data);
	} else {
	    calculate_accelerations_from_star_and_planets(
		temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi, temp_phi_dot,
		data);
	}

	// calculate k1
	double k1_r, k1_phi, k1_r_dot, k1_phi_dot;
	double k2_r, k2_phi, k2_r_dot, k2_phi_dot;

	k1_r_dot = temp_ar;
	k1_phi_dot = temp_aphi;
	k1_r = temp_r_dot;
	k1_phi = temp_phi_dot;

	sk = atoli + rtoli * fabs(temp_r);
	sqr = k1_r / sk;
	dnf += sqr * sqr;
	sqr = temp_r / sk;
	dny += sqr * sqr;

	sk = atoli + rtoli * fabs(temp_phi);
	sqr = k1_phi / sk;
	dnf += sqr * sqr;
	sqr = temp_phi / sk;
	dny += sqr * sqr;

	sk = atoli + rtoli * fabs(temp_r_dot);
	sqr = k1_r_dot / sk;
	dnf += sqr * sqr;
	sqr = temp_r_dot / sk;
	dny += sqr * sqr;

	sk = atoli + rtoli * fabs(temp_phi_dot);
	sqr = k1_phi_dot / sk;
	dnf += sqr * sqr;
	sqr = temp_phi_dot / sk;
	dny += sqr * sqr;

	if ((dnf <= 1.0E-10) || (dny <= 1.0E-10))
	    particles[i].timestep = 1.0E-6;
	else
	    particles[i].timestep = sqrt(dny / dnf) * 0.01;

	// perform an explicit Euler step
	temp_r = particles[i].r + particles[i].timestep * k1_r;
	temp_phi = particles[i].phi + particles[i].timestep * k1_phi;

	temp_r_dot = particles[i].r_dot + particles[i].timestep * k1_r_dot;
	temp_phi_dot =
	    particles[i].phi_dot + particles[i].timestep * k1_phi_dot;

	if (CartesianParticles) {
	    calculate_accelerations_from_star_and_planets_cart(
		temp_ar, temp_aphi, temp_r, temp_phi, data);
	} else {
	    calculate_accelerations_from_star_and_planets(
		temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi, temp_phi_dot,
		data);
	}

	// calculate k2
	k2_r_dot = temp_ar;
	k2_phi_dot = temp_aphi;
	k2_r = temp_r_dot;
	k2_phi = temp_phi_dot;

	// estimate the second derivative of the solution
	der2 = 0.0;

	sk = atoli + rtoli * fabs(temp_r);
	sqr = (k2_r - k1_r) / sk;
	der2 += sqr * sqr;

	sk = atoli + rtoli * fabs(temp_phi);
	sqr = (k2_phi - k1_phi) / sk;
	der2 += sqr * sqr;

	sk = atoli + rtoli * fabs(temp_r_dot);
	sqr = (k2_r_dot - k1_r_dot) / sk;
	der2 += sqr * sqr;

	sk = atoli + rtoli * fabs(temp_phi_dot);
	sqr = (k2_phi_dot - k1_phi_dot) / sk;
	der2 += sqr * sqr;

	der2 = sqrt(der2) / particles[i].timestep;

	// step size is computed such that h**iord * max_d(norm(f0),norm(der2))
	// = 0.01
	der12 = fmax(fabs(der2), sqrt(dnf));
	if (der12 <= 1.0E-15)
	    h1 = fmax(1.0E-6, fabs(particles[i].timestep) * 1.0E-3);
	else
	    h1 = pow(0.01 / der12, 1.0 / (double)iord);
	particles[i].timestep = fmin(100.0 * particles[i].timestep, h1);
    }
}

static void correct_for_self_gravity(t_data &data, const unsigned int i)
{

    const double r = particles[i].get_distance_to_star();
    const double phi = particles[i].get_angle();

    const double eccentricity =
	particles[i].phi_ddot; // reused variable for storage, will be
			       // initialized correctly in this function to 0.0
    const double semi_major_axis = r / (1.0 + eccentricity);

    double sg_radial = 0.0;
    if (parameters::particle_disk_gravity_enabled) {

	// find nearest cells
	unsigned int n_radial_a_minus = 0, n_radial_a_plus = 1;
	unsigned int n_radial_b_minus = 0, n_radial_b_plus = 1;
	unsigned int n_azimuthal_a_minus = 0, n_azimuthal_a_plus = 0;
	unsigned int n_azimuthal_b_minus = 0, n_azimuthal_b_plus = 0;
	find_nearest(n_radial_a_minus, n_radial_a_plus,
		     n_azimuthal_a_minus, n_azimuthal_a_plus, n_radial_b_minus,
		     n_radial_b_plus, n_azimuthal_b_minus, n_azimuthal_b_plus,
		     r, phi);

	sg_radial = interpolate_bilinear_sg(
	    selfgravity::g_radial, n_radial_a_minus, n_radial_a_plus,
	    n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
    }

    double v = sqrt(constants::G * (hydro_center_mass + particles[i].mass) /
			semi_major_axis -
		    sg_radial * semi_major_axis) *
	       sqrt((1.0 - eccentricity) / (1.0 + eccentricity));

    // only need to update velocities, set accelerations to 0 anyway
    if (CartesianParticles) {
	// Beware: cartesian particles still use the names of polar coordinates
	// vx = r_dot
	// vy = phi_dot

	particles[i].r_dot = -v * sin(phi);
	particles[i].phi_dot = v * cos(phi);
    } else {
	particles[i].r_dot = 0.0;
	particles[i].phi_dot = v / r;
    }

    particles[i].r_ddot = 0.0;
    particles[i].phi_ddot = 0.0;
}

static void
init_particle(const unsigned int &i, const unsigned int &id_offset,
	      std::mt19937 &generator,
	      std::uniform_real_distribution<double> &dis_one,
	      std::uniform_real_distribution<double> &dis_twoPi,
	      std::uniform_real_distribution<double> &dis_eccentricity)
{
    double semi_major_axis = f(dis_one(generator), parameters::particle_slope);
    double phi = dis_twoPi(generator);
    double eccentricity = dis_eccentricity(generator);

    /*
    // debug setup
    const int num_particles_per_ring = 10;
    int global_id = i + id_offset;

    double rmin = parameters::particle_minimum_radius;
    double rmax = parameters::particle_maximum_radius;

    // make sure particles are not initialized on an orbit that moves out of
    particle-bounds if(RMAX <
    parameters::particle_maximum_radius*(1.0+parameters::particle_eccentricity))
	    rmax =
    parameters::particle_maximum_radius/(1.0+parameters::particle_eccentricity);

    if(RMIN >
    parameters::particle_minimum_radius*(1.0-parameters::particle_eccentricity))
	    rmin =
    parameters::particle_minimum_radius/(1.0-parameters::particle_eccentricity);


    phi = 2.0*PI / num_particles_per_ring * (global_id/num_particles_per_ring);
    semi_major_axis = (rmax - rmin) *
    double(global_id%num_particles_per_ring)/double(num_particles_per_ring) +
    rmin; eccentricity = 0.0;
    */

    particles[i].radius = parameters::particle_radius;

    if (true) {
	const int particle_type = i % 4; // 4 different particle sizes
	switch (particle_type) {
	case 0: // very small particles
	    particles[i].radius *= 1.0;
	    break;
	case 1: // small particles
	    particles[i].radius *= 5.0;
	    break;
	case 2: // medium particles
	    particles[i].radius *= 10.0;
	    break;
	case 3: // very large particles
	    particles[i].radius *= 50.0;
	    break;
	}
    }

    double volume = 4.0 / 3.0 * PI * pow3(particles[i].radius);

    particles[i].mass = volume * parameters::particle_density;

    double r = semi_major_axis * (1.0 + eccentricity);
    double v = sqrt(constants::G * (hydro_center_mass + particles[i].mass) /
		    semi_major_axis) *
	       sqrt((1.0 - eccentricity) / (1.0 + eccentricity));

    if (CartesianParticles) {
	// Beware: cartesian particles still use the names of polar coordinates
	// x = r
	// y = phi
	// vx = r_dot
	// vy = phi_dot
	particles[i].r = r * cos(phi);
	particles[i].phi = r * sin(phi);

	particles[i].r_dot = -v * sin(phi);
	particles[i].phi_dot = v * cos(phi);
    } else {
	particles[i].r = r;
	particles[i].phi = phi;

	particles[i].r_dot = 0.0;
	particles[i].phi_dot = v / r;
    }

    particles[i].r_ddot = v;
    particles[i].phi_ddot = eccentricity;

    // only needed for adaptive integrator, initialized later
    particles[i].timestep = 0.0;
    particles[i].facold = 0.0;

    particles[i].id = id_offset + i;
}

void init(t_data &data)
{
    // calculate number of initial local particles
    global_number_of_particles = parameters::number_of_particles;
    local_number_of_particles = global_number_of_particles / CPU_Number;

    if ((unsigned int)CPU_Rank <
	global_number_of_particles - CPU_Number * local_number_of_particles) {
	local_number_of_particles++;
    }

    local_r_min = Ra[CPU_Rank == 0 ? GHOSTCELLS_A - 1 : CPUOVERLAP];
    local_r_max = Ra[CPU_Rank == CPU_Highest ? NRadial - GHOSTCELLS_A + 1
					     : NRadial - CPUOVERLAP];

    // create storage
    particles_size = local_number_of_particles;
    particles = (t_particle *)malloc(sizeof(t_particle) * particles_size);

    const unsigned int seed = parameters::random_seed;
    logging::print(LOG_DEBUG "random generator seed: %u\n", seed);

    // random generator and distributions
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> dis_one(0.0, 1.0);
    std::uniform_real_distribution<double> dis_twoPi(0.0, 2.0 * PI);

    std::uniform_real_distribution<double> dis_eccentricity(
	0.0,
	parameters::particle_eccentricity); // for generating eccentricities
					    // same as Marzari & Scholl 2000

    // get number of local particles from all nodes to compute correct offsets
    std::vector<unsigned int> nodes_number_of_particles(CPU_Number);
    MPI_Allgather(&local_number_of_particles, 1, MPI_UNSIGNED,
		  &nodes_number_of_particles[0], 1, MPI_UNSIGNED,
		  MPI_COMM_WORLD);

    // compute local offset
    unsigned int local_offset = 0;
    for (int cpu = 0; cpu < CPU_Rank; ++cpu) {
	local_offset += nodes_number_of_particles[cpu];
    }

    for (unsigned int i = 0; i < local_offset; ++i) {
	init_particle(
	    0, local_offset, generator, dis_one, dis_twoPi,
	    dis_eccentricity); // bring random generator to correct position
    }

    for (unsigned int i = 0; i < local_number_of_particles; ++i) {
	init_particle(i, local_offset, generator, dis_one, dis_twoPi,
		      dis_eccentricity);
    }

    // create MPI datatype
    int mpi_particle_count = 11;
    int mpi_particle_lengths[11] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Aint mpi_particle_offsets[11];
    MPI_Datatype mpi_particle_types[11] = {MPI_UNSIGNED, MPI_DOUBLE, MPI_DOUBLE,
					   MPI_DOUBLE,   MPI_DOUBLE, MPI_DOUBLE,
					   MPI_DOUBLE,   MPI_DOUBLE, MPI_DOUBLE,
					   MPI_DOUBLE,   MPI_DOUBLE};

    // calculate offsets
    MPI_Aint base;
    MPI_Get_address(particles, &base);
    MPI_Get_address(&particles[0].id, mpi_particle_offsets);
    MPI_Get_address(&particles[0].r, mpi_particle_offsets + 1);
    MPI_Get_address(&particles[0].phi, mpi_particle_offsets + 2);
    MPI_Get_address(&particles[0].r_dot, mpi_particle_offsets + 3);
    MPI_Get_address(&particles[0].phi_dot, mpi_particle_offsets + 4);
    MPI_Get_address(&particles[0].r_ddot, mpi_particle_offsets + 5);
    MPI_Get_address(&particles[0].phi_ddot, mpi_particle_offsets + 6);
    MPI_Get_address(&particles[0].mass, mpi_particle_offsets + 7);
    MPI_Get_address(&particles[0].radius, mpi_particle_offsets + 8);
    MPI_Get_address(&particles[0].timestep, mpi_particle_offsets + 9);
    MPI_Get_address(&particles[0].facold, mpi_particle_offsets + 10);

    for (int i = 0; i < 11; ++i) {
	mpi_particle_offsets[i] -= base;
    }

    MPI_Type_create_struct(mpi_particle_count, mpi_particle_lengths,
			   mpi_particle_offsets, mpi_particle_types,
			   &mpi_particle);
    MPI_Type_commit(&mpi_particle);

    // particles might be created on the wrong node, so move them to the correct
    // one :)
    for (int i = 0; i < CPU_Number; ++i) {
	move();
    }

    if (parameters::particle_disk_gravity_enabled) {
	for (unsigned int i = 0; i < local_number_of_particles; ++i) {
	    correct_for_self_gravity(data, i);
	}
    }
    check_tstop(data);

    if (parameters::integrator == parameters::integrator_adaptive)
	init_particle_timestep(data);
}

void restart(unsigned int timestep)
{

    logging::print_master(
	LOG_WARNING
	"Beware: when restarting particles, the user is responsible that the loaded particle file is written with the same coordinate system as the simulation is running on!\n\n");

    FILE *fd;
    char *filename = 0;
    // create filename
    if (asprintf(&filename, "%sparticles%u.dat", OUTPUTDIR, timestep) < 0) {
	logging::print(LOG_ERROR "Not enough memory!\n");
	PersonalExit(1);
    }

    fd = fopen(filename, "r");
    if (fd == nullptr) {
	logging::print_master(
	    LOG_INFO
	    "Can't find file particles%d.dat. Using generated particles.\n",
	    timestep);
	return;
    }
    fseek(fd, 0L, SEEK_END);
    long int size = ftell(fd);
    fclose(fd); // go back to beginning of file

    unsigned int num_particles_in_file = size / sizeof(t_particle);

    if (num_particles_in_file < parameters::number_of_particles) {
	logging::print_master(
	    LOG_ERROR
	    "Warning: Simulation runs with %d particles but only %d particles can be loaded from file!\n\n",
	    parameters::number_of_particles, num_particles_in_file);
	global_number_of_particles = num_particles_in_file;
	parameters::number_of_particles = num_particles_in_file;
	local_number_of_particles = global_number_of_particles / CPU_Number;
	free(particles);

	if ((unsigned int)CPU_Rank <
	    global_number_of_particles -
		CPU_Number * local_number_of_particles) {
	    local_number_of_particles++;
	}

	local_r_min = Ra[CPU_Rank == 0 ? GHOSTCELLS_A - 1 : CPUOVERLAP];
	local_r_max = Ra[CPU_Rank == CPU_Highest ? NRadial - GHOSTCELLS_A + 1
						 : NRadial - CPUOVERLAP];

	// create storage
	particles_size = local_number_of_particles;
	particles = (t_particle *)malloc(sizeof(t_particle) * particles_size);
    }
    if (num_particles_in_file > parameters::number_of_particles) {
	logging::print_master(
	    LOG_WARNING
	    "Warning: particle File contains %d particles but only %d particles can be loaded into Simulation!\n\n",
	    num_particles_in_file, parameters::number_of_particles);
    }

    MPI_File fh;
    MPI_Status status;
    int error, error_class, error_length;
    char error_string[MPI_MAX_ERROR_STRING + 1];

    // try to open file
    error = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,
			  MPI_INFO_NULL, &fh);
    if (error != MPI_SUCCESS) {
	logging::print_master(
	    LOG_ERROR
	    "Error while reading to file '%s'. Check file permissions and IO support of MPI library\n",
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

    logging::print_master(LOG_INFO "Reading file '%s' with %u bytes.\n",
			  filename, size);
    free(filename);

    // get number of local particles from all nodes to compute correct offsets
    std::vector<unsigned int> nodes_number_of_particles(CPU_Number);
    MPI_Allgather(&local_number_of_particles, 1, MPI_UNSIGNED,
		  &nodes_number_of_particles[0], 1, MPI_UNSIGNED,
		  MPI_COMM_WORLD);

    // compute local offset
    unsigned int local_offset = 0;
    for (int cpu = 0; cpu < CPU_Rank; ++cpu) {
	local_offset += nodes_number_of_particles[cpu];
    }

    MPI_File_set_view(fh, 0, mpi_particle, mpi_particle,
		      const_cast<char *>("native"), MPI_INFO_NULL);
    MPI_File_seek(fh, local_offset, MPI_SEEK_SET);
    MPI_File_read_all(fh, particles, local_number_of_particles, mpi_particle,
		      &status);

    // close file
    MPI_File_close(&fh);

    // particles might be loaded to the wrong node, so move them to the correct
    // one
    for (int i = 0; i < CPU_Number; ++i)
	move();
}

void calculate_accelerations_from_star_and_planets(
    double &ar, double &aphi, const double r, const double r_dot,
    const double phi, const double phi_dot, t_data &data)
{
    (void)r_dot;
    (void)phi_dot;

    constexpr double epsilon = 0.005;
    constexpr double epsilon_sq = epsilon * epsilon;

    ar = 0.0;
    aphi = 0.0;

    // planets
    for (unsigned int k = 0;
	 k < data.get_planetary_system().get_number_of_planets(); ++k) {
	const t_planet &planet = data.get_planetary_system().get_planet(k);

	const double r_planet = planet.get_r();
	const double phi_planet = planet.get_phi();

	const double delta_phi = phi - phi_planet;
	double sin_delta_phi;
	double cos_delta_phi;
	sincos(delta_phi, &sin_delta_phi, &cos_delta_phi);

	const double distance_to_planet = sqrt(
	    r * r + r_planet * r_planet - 2 * r * r_planet * cos_delta_phi);
	const double distance_to_planet_pow2_smoothed =
	    distance_to_planet * distance_to_planet + epsilon_sq;

	const double factor =
	    constants::G * planet.get_mass() /
	    (distance_to_planet * distance_to_planet_pow2_smoothed);

	ar -= factor * (r - r_planet * cos_delta_phi);
	aphi -= factor * r_planet * sin_delta_phi / r;
    }
}

void calculate_accelerations_from_star_and_planets_cart(double &ax, double &ay,
							const double x,
							const double y,
							t_data &data)
{
    ax = 0;
    ay = 0;

    constexpr double epsilon = 0.005;
    constexpr double epsilon_sq = epsilon * epsilon;

    // planets
    for (unsigned int k = 0;
	 k < data.get_planetary_system().get_number_of_planets(); ++k) {
	t_planet &planet = data.get_planetary_system().get_planet(k);
	double r2 = pow2(planet.get_x() - x) + pow2(planet.get_y() - y);
	const double r = sqrt(r2);
	r2 += epsilon_sq;
	const double factor = constants::G * planet.get_mass() / (r * r2);

	ax += factor * (planet.get_x() - x);
	ay += factor * (planet.get_y() - y);
    }
}

void calculate_derivitives_from_star_and_planets(double &grav_r_ddot,
						 double &minus_grav_l_dot,
						 const double r,
						 const double phi, t_data &data)
{
    constexpr double epsilon = 0.005;
    constexpr double epsilon_sq = epsilon * epsilon;

    grav_r_ddot = 0.0;
    minus_grav_l_dot = 0.0;

    // planets
    for (unsigned int k = 0;
	 k < data.get_planetary_system().get_number_of_planets(); ++k) {
	const t_planet &planet = data.get_planetary_system().get_planet(k);

	const double r_planet = planet.get_r();
	const double phi_planet = planet.get_phi();
	const double planet_mass = planet.get_mass();

	const double delta_phi = phi - phi_planet;
	double sin_delta_phi;
	double cos_delta_phi;
	sincos(delta_phi, &sin_delta_phi, &cos_delta_phi);

	const double distance_to_planet = sqrt(
	    r * r + r_planet * r_planet - 2.0 * r * r_planet * cos_delta_phi);
	const double distance_to_planet_smoothed_pow2 =
	    distance_to_planet * distance_to_planet + epsilon_sq;

	// direct term
	grav_r_ddot -= constants::G * planet_mass *
		       (r - r_planet * cos_delta_phi) /
		       (distance_to_planet_smoothed_pow2 * distance_to_planet);
	minus_grav_l_dot -=
	    constants::G * planet_mass * r * r_planet * sin_delta_phi /
	    (distance_to_planet_smoothed_pow2 * distance_to_planet);
    }
}

void calculate_derivitives_from_star_and_planets_in_cart(
    double &grav_r_ddot, double &minus_grav_l_dot, const double r,
    const double phi, t_data &data)
{
    constexpr double epsilon = 0.005;
    constexpr double epsilon_sq = epsilon * epsilon;

    double acart[2];
    double acyl[2];

    acart[0] = 0.0;
    acart[1] = 0.0;

    const double x = r * cos(phi);
    const double y = r * sin(phi);

    // planets
    for (unsigned int k = 0;
	 k < data.get_planetary_system().get_number_of_planets(); ++k) {
	const t_planet &planet = data.get_planetary_system().get_planet(k);

	const double x_planet = planet.get_x();
	const double y_planet = planet.get_y();
	const double planet_mass = planet.get_mass();

	const double x_dist = x - x_planet;
	const double y_dist = y - y_planet;

	double dist2 = x_dist * x_dist + y_dist * y_dist;
	const double dist = sqrt(dist2);
	dist2 += epsilon_sq;

	// direct term
	acart[0] += -constants::G * planet_mass * x_dist / (dist * dist2);
	acart[1] += -constants::G * planet_mass * y_dist / (dist * dist2);
    }

    transformCartCyl(acart, acyl, r, phi);
    grav_r_ddot = acyl[0];
    minus_grav_l_dot = acyl[1] * r;
}

static double
interpolate_bilinear(t_polargrid &quantity, bool radial_a_grid,
		     bool azimuthal_a_grid, unsigned int n_radial_minus,
		     unsigned int n_radial_plus, unsigned int n_azimuthal_minus,
		     unsigned int n_azimuthal_plus, double r, double phi)
{
    const double dphi = 2.0 * PI / (double)quantity.get_size_azimuthal();
    const unsigned int last_index = quantity.get_max_azimuthal();

    // values at corners
    double Qmm = quantity(n_radial_minus, n_azimuthal_minus);
    double Qpm = quantity(n_radial_plus, n_azimuthal_minus);
    double Qmp = quantity(n_radial_minus, n_azimuthal_plus);
    double Qpp = quantity(n_radial_plus, n_azimuthal_plus);

    double rm = radial_a_grid ? Ra[n_radial_minus] : Rb[n_radial_minus];
    double rp = radial_a_grid ? Ra[n_radial_plus] : Rb[n_radial_plus];

    double phim;
    double phip;

    if (azimuthal_a_grid) {

	phim = n_azimuthal_minus * dphi;

	// particle is in cell with n_azimuthal = N_azimuthal_max - 1
	// next cells is at n_azimuthal = 0, but measure distance from 2*PI
	// instead
	if (n_azimuthal_plus == 0) {
	    phip = double(n_azimuthal_minus + 1) * dphi;
	} else {
	    phip = n_azimuthal_plus * dphi;
	}
    } else {
	if (n_azimuthal_minus == last_index) {
	    if (phi < PI) {
		// particle is inside cell with n_azimuthal = 0
		// previous cell is at N_azimuthal_max - 1, but measure distance
		// from -0.5*dphi
		phim = -0.5 * dphi;
		phip = 0.5 * dphi;
	    } else {
		// particle is inside cell with n_azimuthal = N_azimuthal_max -
		// 1 next cell is at N_azimuthal = 0, but measure distance from
		// 2*PI+0.5*dphi
		phim = (n_azimuthal_minus + 0.5) * dphi;
		phip = (n_azimuthal_minus + 1.5) * dphi;
	    }
	} else {
	    phim = (n_azimuthal_minus + 0.5) * dphi;
	    phip = (n_azimuthal_plus + 0.5) * dphi;
	}
    }

    // double Qm =
    // ((rm*phip-rm*phi)*Qmm+(rm*phi-rm*phim)*Qmp)/(rm*phip-rm*phim); double Qp
    // =
    // ((rp*phip-rp*phi)*Qpm+(rp*phi-rp*phim)*Qpp)/(rp*phip-rp*phim);
    // mathematically the same, but less computations
    const double Qm = ((phip - phi) * Qmm + (phi - phim) * Qmp) / dphi;
    const double Qp = ((phip - phi) * Qpm + (phi - phim) * Qpp) / dphi;

    const double Q = ((rp - r) * Qm + (r - rm) * Qp) / (rp - rm);

    return Q;
}

static void calculate_tstop(const double r, const double phi,
			    const double r_dot, const double phi_dot,
			    t_data &data, const double radius,
			    double &minus_r_dotel_r, double &minus_l_rel,
			    double &tstop)
{

    unsigned int n_radial_a_minus = 0, n_radial_a_plus = 1,
		 n_radial_b_minus = 0, n_radial_b_plus = 1;
    unsigned int n_azimuthal_a_minus = 0, n_azimuthal_a_plus = 0,
		 n_azimuthal_b_minus = 0, n_azimuthal_b_plus = 0;
	find_nearest(n_radial_a_minus, n_radial_a_plus,
		 n_azimuthal_a_minus, n_azimuthal_a_plus, n_radial_b_minus,
		 n_radial_b_plus, n_azimuthal_b_minus, n_azimuthal_b_plus, r,
		 phi);

    // calculate gas quantities at the particle location
    const double rho = interpolate_bilinear(
	data[t_data::RHO], false, false, n_radial_b_minus, n_radial_b_plus,
	n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
    const double temperature = interpolate_bilinear(
	data[t_data::TEMPERATURE], false, false, n_radial_b_minus,
	n_radial_b_plus, n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
    const double vg_radial = interpolate_bilinear(
	data[t_data::V_RADIAL], true, false, n_radial_a_minus, n_radial_a_plus,
	n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
    const double vg_azimuthal = interpolate_bilinear(
	data[t_data::V_AZIMUTHAL], false, true, n_radial_b_minus,
	n_radial_b_plus, n_azimuthal_a_minus, n_azimuthal_a_plus, r, phi);
    const double m0 = parameters::MU * constants::m_u.get_code_value();
    const double vthermal =
	sqrt(8.0 * constants::k_B.get_code_value() * temperature / (PI * m0));

    // calculate relative velocities
    minus_r_dotel_r = vg_radial - r_dot;
    const double vrel_phi = vg_azimuthal - phi_dot * r;
    minus_l_rel = r * vrel_phi;

    const double vrel = sqrt(pow2(minus_r_dotel_r) + pow2(vrel_phi));

    // a0 = 1.5e-8 cm for molecular hydrogen
    double sigma = (PI * pow2(1.5e-8 / units::length.get_cgs_factor()));
    double nu = 1.0 / 3.0 * m0 * vthermal / sigma;

    // calculate Reynolds number
    double reynolds = 2.0 * rho * radius * vrel / nu;

    // calculate coefficient
    double Cd;
    if (reynolds < 1.0) {
	Cd = 24.0 / reynolds;
    } else if (reynolds < 800.0) {
	Cd = 24.0 * pow(reynolds, -0.6);
    } else {
	Cd = 0.44;
    }

    // mean free path for molecular hydrogen (see Haghighipour & Boss, 2003 eq.
    // 20)
    double l = 4.72e-9 / rho;
    double f = radius / (radius + l);

    // ***************************************************************************
    // Combined Epstein + Stokes drag regimes (see Haghighipour & Boss, 2003 eq.
    // 8)
    double term = rho * ((1.0 - f) * vthermal + 3.0 / 8.0 * f * Cd * vrel);
    // ***************************************************************************
    // Epstein only drag regime
    // double term = rho*vthermal;
    // ***************************************************************************
    // Stokes only drag regime
    // double term = rho*vthermal*9.0/4.0*l/particles[i].radius;
    // ***************************************************************************

    // Stopping time
    tstop = radius * parameters::particle_density / term;
}

static void calculate_tstop2(const double r, const double phi,
			     const double r_dot, const double phi_dot,
			     t_data &data, const double radius,
			     double &minus_r_dotel_r, double &minus_l_rel,
			     double &tstop, const double r0, const double l0)
{

    // r has been updated since last  move(), so we need to confirm that the
    // particle is still inside the domain the gas values are determined at the
    // edge of the domain if particle is outside the domain
    double r_tmp = fmax(r, local_r_min);
    r_tmp = fmin(r, local_r_max);

    unsigned int n_radial_a_minus = 0, n_radial_a_plus = 1,
		 n_radial_b_minus = 0, n_radial_b_plus = 1;
    unsigned int n_azimuthal_a_minus = 0, n_azimuthal_a_plus = 0,
		 n_azimuthal_b_minus = 0, n_azimuthal_b_plus = 0;
	find_nearest(n_radial_a_minus, n_radial_a_plus,
		 n_azimuthal_a_minus, n_azimuthal_a_plus, n_radial_b_minus,
		 n_radial_b_plus, n_azimuthal_b_minus, n_azimuthal_b_plus,
		 r_tmp, phi);

    // calculate gas quantities at the particle location
    const double rho = interpolate_bilinear(
	data[t_data::RHO], false, false, n_radial_b_minus, n_radial_b_plus,
	n_azimuthal_b_minus, n_azimuthal_b_plus, r_tmp, phi);
    const double temperature = interpolate_bilinear(
	data[t_data::TEMPERATURE], false, false, n_radial_b_minus,
	n_radial_b_plus, n_azimuthal_b_minus, n_azimuthal_b_plus, r_tmp, phi);
    const double vg_radial = interpolate_bilinear(
	data[t_data::V_RADIAL], true, false, n_radial_a_minus, n_radial_a_plus,
	n_azimuthal_b_minus, n_azimuthal_b_plus, r_tmp, phi);
    const double vg_azimuthal = interpolate_bilinear(
	data[t_data::V_AZIMUTHAL], false, true, n_radial_b_minus,
	n_radial_b_plus, n_azimuthal_a_minus, n_azimuthal_a_plus, r_tmp, phi);
    const double m0 = parameters::MU * constants::m_u.get_code_value();
    const double vthermal =
	sqrt(8.0 * constants::k_B.get_code_value() * temperature / (PI * m0));

    // calculate relative velocities
    minus_r_dotel_r = vg_radial - r_dot;
    const double vrel_phi = vg_azimuthal - phi_dot * r0;

    const double vrel = sqrt(pow2(minus_r_dotel_r) + pow2(vrel_phi));

    minus_l_rel = r * vg_azimuthal - l0;

    // a0 = 1.5e-8 cm for molecular hydrogen
    double sigma = (PI * pow2(1.5e-8 / units::length.get_cgs_factor()));
    double nu = 1.0 / 3.0 * m0 * vthermal / sigma;

    // calculate Reynolds number
    double reynolds = 2.0 * rho * radius * vrel / nu;

    // calculate coefficient
    double Cd;
    if (reynolds < 1.0) {
	Cd = 24.0 / reynolds;
    } else if (reynolds < 800.0) {
	Cd = 24.0 * pow(reynolds, -0.6);
    } else {
	Cd = 0.44;
    }

    // mean free path for molecular hydrogen (see Haghighipour & Boss, 2003 eq.
    // 20)
    double l = 4.72e-9 / rho;
    double f = radius / (radius + l);

    // ***************************************************************************
    // Combined Epstein + Stokes drag regimes (see Haghighipour & Boss, 2003 eq.
    // 8)
    double term = rho * ((1.0 - f) * vthermal + 3.0 / 8.0 * f * Cd * vrel);
    // ***************************************************************************
    // Epstein only drag regime
    // double term = rho*vthermal;
    // ***************************************************************************
    // Stokes only drag regime
    // double term = rho*vthermal*9.0/4.0*l/particles[i].radius;
    // ***************************************************************************

    // Stopping time
    tstop = radius * parameters::particle_density / term;
}

// confirm that tstop < dt/10, else make user change integrator
// tstop > dt/10 can cause numerical instabilities for the explicit integrator
void check_tstop(t_data &data)
{

    double dt = DT;

    double local_gas_time_step_cfl = 1.0;
    double global_gas_time_step_cfl = 1.0;
    CommunicateBoundaries(&data[t_data::DENSITY], &data[t_data::V_RADIAL],
			  &data[t_data::V_AZIMUTHAL], &data[t_data::ENERGY]);
    local_gas_time_step_cfl =
	condition_cfl(data, data[t_data::V_RADIAL], data[t_data::V_AZIMUTHAL],
		      data[t_data::SOUNDSPEED], DT);

    MPI_Allreduce(&local_gas_time_step_cfl, &global_gas_time_step_cfl, 1,
		  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    dt = DT / global_gas_time_step_cfl;

    compute_rho(data, true);
    compute_temperature(data, true);
    for (unsigned int i = 0; i < local_number_of_particles; ++i) {
	const double radius = particles[i].radius;
	const double r = particles[i].get_distance_to_star();
	const double phi = particles[i].get_angle();
	const double r_dot = particles[i].get_r_dot();
	const double phi_dot = particles[i].get_phi_dot();

	unsigned int n_radial_a_minus = 0, n_radial_a_plus = 1,
		     n_radial_b_minus = 0, n_radial_b_plus = 1;
	unsigned int n_azimuthal_a_minus = 0, n_azimuthal_a_plus = 0,
		     n_azimuthal_b_minus = 0, n_azimuthal_b_plus = 0;

	find_nearest(n_radial_a_minus, n_radial_a_plus,
		     n_azimuthal_a_minus, n_azimuthal_a_plus, n_radial_b_minus,
		     n_radial_b_plus, n_azimuthal_b_minus, n_azimuthal_b_plus,
		     r, phi);

	// calculate gas quantities at the particle location
	const double rho = interpolate_bilinear(
	    data[t_data::RHO], false, false, n_radial_b_minus, n_radial_b_plus,
	    n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
	const double temperature = interpolate_bilinear(
	    data[t_data::TEMPERATURE], false, false, n_radial_b_minus,
	    n_radial_b_plus, n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
	const double vg_radial = interpolate_bilinear(
	    data[t_data::V_RADIAL], true, false, n_radial_a_minus,
	    n_radial_a_plus, n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
	const double vg_azimuthal = interpolate_bilinear(
	    data[t_data::V_AZIMUTHAL], false, true, n_radial_b_minus,
	    n_radial_b_plus, n_azimuthal_a_minus, n_azimuthal_a_plus, r, phi);
	const double m0 = parameters::MU * constants::m_u.get_code_value();
	const double vthermal = sqrt(8.0 * constants::k_B.get_code_value() *
				     temperature / (PI * m0));

	// calculate relative velocities
	double minus_r_dotel_r = vg_radial - r_dot;
	const double vrel_phi = vg_azimuthal - phi_dot * r;

	const double vrel = sqrt(pow2(minus_r_dotel_r) + pow2(vrel_phi));

	// a0 = 1.5e-8 cm for molecular hydrogen
	double sigma = (PI * pow2(1.5e-8 / units::length.get_cgs_factor()));
	double nu = 1.0 / 3.0 * m0 * vthermal / sigma;

	// calculate Reynolds number
	double reynolds = 2.0 * rho * radius * vrel / nu;

	// calculate coefficient
	double Cd;
	if (reynolds < 1.0) {
	    Cd = 24.0 / reynolds;
	} else if (reynolds < 800.0) {
	    Cd = 24.0 * pow(reynolds, -0.6);
	} else {
	    Cd = 0.44;
	}

	// mean free path for molecular hydrogen (see Haghighipour & Boss, 2003
	// eq. 20)
	double l = 4.72e-9 / rho;
	double f = radius / (radius + l);

	// ***************************************************************************
	// Combined Epstein + Stokes drag regimes (see Haghighipour & Boss, 2003
	// eq. 8)
	double term = rho * ((1.0 - f) * vthermal + 3.0 / 8.0 * f * Cd * vrel);
	// ***************************************************************************
	// Epstein only drag regime
	// double term = rho*vthermal;
	// ***************************************************************************
	// Stokes only drag regime
	// double term = rho*vthermal*9.0/4.0*l/particles[i].radius;
	// ***************************************************************************

	// Stopping time
	double tstop = radius * parameters::particle_density / term;

	if (tstop < 10.0 * dt && parameters::particle_gas_drag_enabled &&
	    (parameters::integrator == parameters::integrator_explicit ||
	     parameters::integrator == parameters::integrator_adaptive)) {
	    logging::print_master(
		LOG_ERROR
		"Particle stopping time too small for explicit integrator! Use the semiimplicit or implicit integrator instead!");
	    PersonalExit(1);
	}

	if (tstop < 0.005 * dt && parameters::particle_gas_drag_enabled &&
	    parameters::integrator == parameters::integrator_semiimplicit) {
	    logging::print_master(
		LOG_ERROR
		"Particle stopping time too small for semiimplicit integrator! Use the implicit integrator instead!");
	    PersonalExit(1);
	}
    }
}

// apply disk feedback on primary onto the particles
void update_velocities_from_indirect_term(const double dt)
{
    // Naming of r and phi weird!!!!!!!
    for (unsigned int i = 0; i < local_number_of_particles; ++i) {

	double indirect_q1_dot;
	double indirect_q2_dot;
	if (CartesianParticles) {
	    indirect_q1_dot = IndirectTerm.x * dt + particles[i].r_dot;
	    indirect_q2_dot = IndirectTerm.y * dt + particles[i].phi_dot;
	} else {
	    double r = particles[i].r;
	    double phi = particles[i].phi;

	    const double r_dot = particles[i].r_dot;
	    const double phi_dot = particles[i].phi_dot;

	    indirect_q1_dot = r_dot + dt * (IndirectTerm.x * cos(phi) +
					    IndirectTerm.y * sin(phi));
	    indirect_q2_dot =
		phi_dot +
		dt * (-IndirectTerm.x * sin(phi) + IndirectTerm.y * cos(phi)) /
		    r;
	}

	particles[i].r_dot = indirect_q1_dot;
	particles[i].phi_dot = indirect_q2_dot;
    }
}

void update_velocities_from_gas_drag_cart(t_data &data, double dt)
{

    compute_rho(data, true);

    for (unsigned int i = 0; i < local_number_of_particles; ++i) {
	double r = particles[i].get_distance_to_star();
	double phi = particles[i].get_angle();

	// check if particle has left disc
	if ((r < RMIN) || (r > RMAX)) {
	    continue;
	}

	// find nearest cells
	unsigned int n_radial_a_minus = 0, n_radial_a_plus = 1;
	unsigned int n_radial_b_minus = 0, n_radial_b_plus = 1;
	unsigned int n_azimuthal_a_minus = 0, n_azimuthal_a_plus = 0;
	unsigned int n_azimuthal_b_minus = 0, n_azimuthal_b_plus = 0;
	find_nearest(n_radial_a_minus, n_radial_a_plus,
		     n_azimuthal_a_minus, n_azimuthal_a_plus, n_radial_b_minus,
		     n_radial_b_plus, n_azimuthal_b_minus, n_azimuthal_b_plus,
		     r, phi);

	// calculate quantities needed
	double rho = interpolate_bilinear(
	    data[t_data::RHO], false, false, n_radial_b_minus, n_radial_b_plus,
	    n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
	double vg_radial = interpolate_bilinear(
	    data[t_data::V_RADIAL], true, false, n_radial_a_minus,
	    n_radial_a_plus, n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
	double vg_azimuthal = interpolate_bilinear(
	    data[t_data::V_AZIMUTHAL], false, true, n_radial_b_minus,
	    n_radial_b_plus, n_azimuthal_a_minus, n_azimuthal_a_plus, r, phi);
	double temperature = interpolate_bilinear(
	    data[t_data::TEMPERATURE], false, false, n_radial_b_minus,
	    n_radial_b_plus, n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);

	// calculate gas velocities in cartesian coordinates
	double vg_x =
	    cos(phi) * vg_radial - sin(phi) * (vg_azimuthal + OmegaFrame * r);
	double vg_y =
	    sin(phi) * vg_radial + cos(phi) * (vg_azimuthal + OmegaFrame * r);

	// particles store cartesian data
	double &vx = particles[i].r_dot;
	double &vy = particles[i].phi_dot;

	// calculate relative velocities
	double vrel_x = vx - vg_x;
	double vrel_y = vy - vg_y;
	double vrel = sqrt(pow2(vrel_x) + pow2(vrel_y));

	double m0 = parameters::MU * constants::m_u.get_code_value();
	double vthermal = sqrt(8.0 * constants::k_B.get_code_value() *
			       temperature / (PI * m0));
	// a0 = 1.5e-8 cm for molecular hydrogen
	double sigma = (PI * pow2(1.5e-8 / units::length.get_cgs_factor()));
	double nu = 1.0 / 3.0 * m0 * vthermal / sigma;

	// calculate Reynolds number
	double reynolds = 2.0 * rho * particles[i].radius * vrel / nu;

	// calculate coefficient
	double Cd;
	if (reynolds < 1.0) {
	    Cd = 24.0 / reynolds;
	} else if (reynolds < 800.0) {
	    Cd = 24.0 * pow(reynolds, -0.6);
	} else {
	    Cd = 0.44;
	}

	double fdrag_temp =
	    -0.5 * Cd * PI * pow2(particles[i].radius) * rho * vrel;
	double fdrag_x = fdrag_temp * vrel_x;
	double fdrag_y = fdrag_temp * vrel_y;

	// update velocities
	vx += dt * fdrag_x / particles[i].mass;
	vy += dt * fdrag_y / particles[i].mass;

	if (parameters::particle_disk_gravity_enabled) {
	    update_velocity_from_disk_gravity_cart(
		n_radial_a_minus, n_radial_a_plus, n_azimuthal_b_minus,
		n_azimuthal_b_plus, r, phi, i, dt);
	}
    }

    /*
    // For testing purpose only, very slow
    if (parameters::particle_disk_gravity_enabled) {
	    update_velocity_from_disk_gravity_cart_old(data, dt);
    }
    */
}

void update_velocities_from_gas_drag(t_data &data, double dt)
{

    compute_rho(data, true);

    for (unsigned int i = 0; i < local_number_of_particles; ++i) {
	double r = particles[i].get_distance_to_star();
	double phi = particles[i].get_angle();

	// check if particle has left disc
	if ((r < RMIN) || (r > RMAX)) {
	    continue;
	}

	// find nearest cells
	unsigned int n_radial_a_minus = 0, n_radial_a_plus = 1;
	unsigned int n_radial_b_minus = 0, n_radial_b_plus = 1;
	unsigned int n_azimuthal_a_minus = 0, n_azimuthal_a_plus = 0;
	unsigned int n_azimuthal_b_minus = 0, n_azimuthal_b_plus = 0;
	find_nearest(n_radial_a_minus, n_radial_a_plus,
		     n_azimuthal_a_minus, n_azimuthal_a_plus, n_radial_b_minus,
		     n_radial_b_plus, n_azimuthal_b_minus, n_azimuthal_b_plus,
		     r, phi);

	// calculate quantities needed
	const double rho = interpolate_bilinear(
	    data[t_data::RHO], false, false, n_radial_b_minus, n_radial_b_plus,
	    n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
	const double vg_radial = interpolate_bilinear(
	    data[t_data::V_RADIAL], true, false, n_radial_a_minus,
	    n_radial_a_plus, n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
	const double vg_azimuthal = interpolate_bilinear(
	    data[t_data::V_AZIMUTHAL], false, true, n_radial_b_minus,
	    n_radial_b_plus, n_azimuthal_a_minus, n_azimuthal_a_plus, r, phi);
	const double temperature = interpolate_bilinear(
	    data[t_data::TEMPERATURE], false, false, n_radial_b_minus,
	    n_radial_b_plus, n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);

	// calculate relative velocities
	const double vrel_r = particles[i].r_dot - vg_radial;
	const double vrel_phi =
	    r * (particles[i].phi_dot - OmegaFrame) - vg_azimuthal;
	const double vrel = sqrt(pow2(vrel_r) + pow2(vrel_phi));

	const double m0 = parameters::MU * constants::m_u.get_code_value();
	const double vthermal = sqrt(8.0 * constants::k_B.get_code_value() *
				     temperature / (PI * m0));
	// a0 = 1.5e-8 cm for molecular hydrogen
	const double sigma =
	    (PI * pow2(1.5e-8 / units::length.get_cgs_factor()));
	const double nu = 1.0 / 3.0 * m0 * vthermal / sigma;

	// calculate Reynolds number
	const double reynolds = 2.0 * rho * particles[i].radius * vrel / nu;

	// calculate coefficient
	double Cd;
	if (reynolds < 1.0) {
	    Cd = 24.0 / reynolds;
	} else if (reynolds < 800.0) {
	    Cd = 24.0 * pow(reynolds, -0.6);
	} else {
	    Cd = 0.44;
	}

	const double fdrag_temp = -0.5 * Cd * PI * pow2(particles[i].radius) *
				  rho * vrel / particles[i].mass;

	const double fdrag_r = dt * fdrag_temp * vrel_r;
	const double fdrag_phi = dt * fdrag_temp * vrel_phi / r;

	// update velocities
	particles[i].r_dot += fdrag_r;
	particles[i].phi_dot += fdrag_phi;

	if (parameters::particle_disk_gravity_enabled) {
	    update_velocity_from_disk_gravity(
		n_radial_a_minus, n_radial_a_plus, n_azimuthal_b_minus,
		n_azimuthal_b_plus, r, phi, i, dt);
	}
    }
}

void update_velocity_from_disk_gravity(const int n_radial_a_minus,
				       const int n_radial_a_plus,
				       const int n_azimuthal_b_minus,
				       const int n_azimuthal_b_plus,
				       const double r, const double phi,
				       const int particle_id, const double dt)
{

    const double sg_radial = interpolate_bilinear_sg(
	selfgravity::g_radial, n_radial_a_minus, n_radial_a_plus,
	n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
    const double sg_azimuthal = interpolate_bilinear_sg(
	selfgravity::g_azimuthal, n_radial_a_minus, n_radial_a_plus,
	n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);

    particles[particle_id].r_dot += dt * sg_radial;
    particles[particle_id].phi_dot += dt * sg_azimuthal / r;
}

void integrate(t_data &data, const double dt)
{
    switch (parameters::integrator) {
    case parameters::integrator_explicit: {
	integrate_explicit(data, dt);
	break;
    }
    case parameters::integrator_adaptive: {
	integrate_explicit_adaptive(data, dt);
	break;
    }
    case parameters::integrator_semiimplicit: {
	integrate_semiimplicit(data, dt);
	break;
    }
    case parameters::integrator_implicit: {
	integrate_implicit(data, dt);
	break;
    }
    default: {
	die("No particle integrator");
    }
    }
}

void update_velocity_from_disk_gravity_cart(
    const int n_radial_a_minus, const int n_radial_a_plus,
    const int n_azimuthal_b_minus, const int n_azimuthal_b_plus, const double r,
    const double phi, const int particle_id, const double dt)
{

    const double sg_radial = interpolate_bilinear_sg(
	selfgravity::g_radial, n_radial_a_minus, n_radial_a_plus,
	n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
    const double sg_azimuthal = interpolate_bilinear_sg(
	selfgravity::g_azimuthal, n_radial_a_minus, n_radial_a_plus,
	n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);

    double &vx = particles[particle_id].r_dot;
    double &vy = particles[particle_id].phi_dot;

    vx += dt * (cos(phi) * sg_radial - sin(phi) * sg_azimuthal);
    vy += dt * (sin(phi) * sg_radial + cos(phi) * sg_azimuthal);
}

void update_velocity_from_disk_gravity_cart_old(t_data &data, double dt)
{
    int *number_of_particles = (int *)malloc(sizeof(int) * CPU_Number);
    int *particle_offsets = (int *)malloc(sizeof(int) * CPU_Number);
    t_particle *all_particles =
	(t_particle *)malloc(sizeof(t_particle) * global_number_of_particles);
    double *force_x =
	(double *)malloc(sizeof(double) * global_number_of_particles);
    double *force_y =
	(double *)malloc(sizeof(double) * global_number_of_particles);

    // get particle amounts from all nodes
    MPI_Allgather(&local_number_of_particles, 1, MPI_INT, number_of_particles,
		  1, MPI_INT, MPI_COMM_WORLD);

    // calculate offsets
    for (int i = 0; i < CPU_Number; ++i) {
	particle_offsets[i] = 0;
	for (int j = 0; j < i; ++j) {
	    particle_offsets[i] += number_of_particles[j];
	}
    }

    // get all particles
    MPI_Allgatherv(particles, local_number_of_particles, mpi_particle,
		   all_particles, number_of_particles, particle_offsets,
		   mpi_particle, MPI_COMM_WORLD);

    for (unsigned int i = 0; i < global_number_of_particles; ++i) {
	force_x[i] = 0;
	force_y[i] = 0;
    }

    double dphi = 2.0 * PI / (double)data[t_data::DENSITY].get_size_azimuthal();
    for (unsigned int n_radial = Zero_or_active; n_radial < Max_or_active;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
	     ++n_azimuthal) {
	    double cell_angle = n_azimuthal * dphi;
	    double cell_x = Rmed[n_radial] * cos(cell_angle);
	    double cell_y = Rmed[n_radial] * sin(cell_angle);
	    double cell_mass =
		Surf[n_radial] * data[t_data::DENSITY](n_radial, n_azimuthal);
	    for (unsigned int i = 0; i < global_number_of_particles; ++i) {
		const double &x = all_particles[i].r;
		const double &y = all_particles[i].phi;
		double smoothing = compute_smoothing_isothermal(
		    all_particles[i].get_distance_to_star());
		double d_x = cell_x - x;
		double d_y = cell_y - y;
		double invdist3 =
		    pow(pow2(d_x) + pow2(d_y) + pow2(smoothing), -3.0 / 2.0);

		force_x[i] += constants::G * cell_mass * d_x * invdist3;
		force_y[i] += constants::G * cell_mass * d_y * invdist3;
	    }
	}
    }

    // send force back
    MPI_Allreduce(MPI_IN_PLACE, force_x, global_number_of_particles, MPI_DOUBLE,
		  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, force_y, global_number_of_particles, MPI_DOUBLE,
		  MPI_SUM, MPI_COMM_WORLD);

    // update particles
    unsigned int offset = particle_offsets[CPU_Rank];
    for (unsigned int i = 0; i < local_number_of_particles; ++i) {
	double &vx = particles[i].r_dot;
	double &vy = particles[i].phi_dot;
	vx += dt * force_x[offset + i];
	vy += dt * force_y[offset + i];
    }

    free(number_of_particles);
    free(particle_offsets);
    free(all_particles);
    free(force_x);
    free(force_y);
}

void integrate_semiimplicit(t_data &data, const double dt)
{

    if (parameters::disk_feedback) {
	update_velocities_from_indirect_term(dt);
    }

    // Semi implicit integrator in cylindrical coordinates (see Zhu et al. 2014,
    // eqs. A4-A12)
    if (parameters::particle_gas_drag_enabled) {
	compute_rho(data, true);
    }

    for (unsigned int i = 0; i < local_number_of_particles; ++i) {

	// initialize
	const double r0 = particles[i].r;
	const double phi0 = particles[i].phi;

	const double r_dot0 = particles[i].r_dot;
	const double phi_dot0 = particles[i].phi_dot;

	const double l0 = r0 * r0 * phi_dot0;

	const double hfdt = 0.5 * dt;
	double tstop = dt;
	double dt1 = dt;
	double vrel_r;

	// Half-drift
	const double r_dot1 = r_dot0;
	const double l1 = l0;
	const double phi_dot1 = phi_dot0; // follows from  l1 = l0

	const double r1 = r0 + r_dot0 * hfdt;
	double phi1 = phi0 + 0.5 * (l0 / pow2(r0) + l0 / pow2(r1)) * hfdt;
	check_angle(phi1);

	// Kick
	const double r2 = r1;
	const double phi2 = phi1;
	double minus_l_rel;

	if (parameters::particle_gas_drag_enabled) {
	    calculate_tstop2(r1, phi1, r_dot1, phi_dot1, data,
			     particles[i].radius, vrel_r, minus_l_rel, tstop,
			     r0, l0);
	    dt1 = dt / (1.0 + hfdt / tstop);
	}
	double grav_r_ddot;
	double minus_grav_l_dot;
	if (ParticlesInCartesian) {
	    calculate_derivitives_from_star_and_planets_in_cart(
		grav_r_ddot, minus_grav_l_dot, r1, phi1, data);
	} else {

	    calculate_derivitives_from_star_and_planets(
		grav_r_ddot, minus_grav_l_dot, r1, phi1, data);
	}

	double l2 = l1 + minus_grav_l_dot * dt1;
	if (parameters::particle_gas_drag_enabled) {
	    l2 += minus_l_rel / tstop * dt1;
	}

	// Updating radial velocity
	double r_dot2 = r_dot1;
	r_dot2 += grav_r_ddot * dt1;
	r_dot2 += 0.5 * (l1 * l1 + l2 * l2) / (r1 * r1 * r1) * dt1;

	if (parameters::particle_gas_drag_enabled) {
	    r_dot2 += vrel_r / tstop * dt1;
	}

	// Half-drift
	const double r3 = r1 + r_dot2 * hfdt;
	const double phi3 = phi2 + 0.5 * (l2 / pow2(r2) + l2 / pow2(r3)) * hfdt;

	particles[i].r_dot = r_dot2;
	particles[i].r = r3;
	// particles[i].phi  = phi1 + 0.5*(phi_dot0 + l2/pow2(particles[i].r)) *
	// hfdt; particles[i].phi  = phi1 +
	// 0.5*(l1/pow2(r1)+l2/pow2(particles[i].r)) * hfdt;
	particles[i].phi = phi3;
	check_angle(particles[i].phi);
	particles[i].phi_dot = l2 / pow2(particles[i].r);
    }
    move();
}

void integrate_implicit(t_data &data, const double dt)
{
    if (parameters::disk_feedback) {
	update_velocities_from_indirect_term(dt);
    }

    // Fully implicit integrator in cylindrical coordinates (see Zhu et al.
    // 2014, eqs. A15-A18)
    double r1, phi1, l1, r_dot1, tstop1, tstop0, dt0, dt1;
    double r_ddot0, minus_l_dot0, r_ddot1, minus_l_dot1, hfdt, minus_r_dot_rel0,
	minus_l_rel0, minus_r_dot_rel1, minus_l_rel1;

    // initialize with failsave values to suppress compiler warning
    // "-Wmaybe-uninitialized"
    minus_l_rel0 = 0.0;
    minus_r_dot_rel0 = 0.0;
    tstop0 = 1e+300;
    if (parameters::particle_gas_drag_enabled) {
	compute_rho(data, true);
    }

    for (unsigned int i = 0; i < local_number_of_particles; ++i) {
	const double r0 = particles[i].r;
	const double phi0 = particles[i].phi;
	const double r_dot0 = particles[i].r_dot;
	const double phi_dot0 = particles[i].phi_dot;
	const double l0 = r0 * r0 * phi_dot0;
	dt1 = dt;
	dt0 = dt1;
	hfdt = 0.5 * dt;

	// Half-drift
	r1 = r0 + r_dot0 * dt;
	phi1 = phi0 + 0.5 * (l0 / pow2(r0) + l0 / pow2(r1)) * dt;
	check_angle(phi1);

	// Kick
	// Current position
	if (parameters::particle_gas_drag_enabled) {
	    calculate_tstop(r0, phi0, r_dot0, phi_dot0, data,
			    particles[i].radius, minus_r_dot_rel0, minus_l_rel0,
			    tstop0);
	    dt0 = 1.0 + dt / tstop0;
	}

	if (ParticlesInCartesian) {
	    calculate_derivitives_from_star_and_planets_in_cart(
		r_ddot0, minus_l_dot0, r0, phi0, data);
	} else {
	    calculate_derivitives_from_star_and_planets(r_ddot0, minus_l_dot0,
							r0, phi0, data);
	}
	// Predicted position
	if (parameters::particle_gas_drag_enabled) {
	    calculate_tstop2(r1, phi1, r_dot0, phi_dot0, data,
			     particles[i].radius, minus_r_dot_rel1,
			     minus_l_rel1, tstop1, r0, l0);
	    dt1 = hfdt / (1.0 + hfdt / tstop0 + hfdt / tstop1 +
			  hfdt * dt / (tstop0 * tstop1));
	}

	if (ParticlesInCartesian) {
	    calculate_derivitives_from_star_and_planets_in_cart(
		r_ddot1, minus_l_dot1, r1, phi1, data);
	} else {
	    calculate_derivitives_from_star_and_planets(r_ddot1, minus_l_dot1,
							r1, phi1, data);
	}

	l1 = l0;
	l1 += (minus_l_dot0 + minus_l_dot1 * dt0) * dt1;
	if (parameters::particle_gas_drag_enabled) {
	    l1 += (minus_l_rel0 / tstop0 + minus_l_rel1 / tstop1 * dt0) * dt1;
	}

	r_dot1 = r_dot0;
	r_dot1 += (r_ddot0 + r_ddot1 * dt0) * dt1;
	r_dot1 += (pow2(l0) / pow3(r0) + pow2(l1) / pow3(r1) * dt0) * dt1;
	if (parameters::particle_gas_drag_enabled) {
	    r_dot1 +=
		(minus_r_dot_rel0 / tstop0 + minus_r_dot_rel1 / tstop1 * dt0) *
		dt1;
	}

	// Half-drift
	const double r2 = r0 + (r_dot0 + r_dot1) * hfdt;
	const double phi2 = phi0 + (l0 / pow2(r0) + l1 / pow2(r2)) * hfdt;
	particles[i].r_dot = r_dot1;
	particles[i].r = r2;
	particles[i].phi = phi2;
	check_angle(particles[i].phi);
	particles[i].phi_dot = l1 / pow2(particles[i].r);
    }

    move();
}

void integrate_explicit_adaptive(t_data &data, const double dt)
{

    if (CartesianParticles) {
	// disk gravity on particles is inside gas_drag function
	if (parameters::particle_gas_drag_enabled)
	    update_velocities_from_gas_drag_cart(data, dt);
    } else {
	// disk gravity on particles is inside gas_drag function
	if (parameters::particle_gas_drag_enabled)
	    update_velocities_from_gas_drag(data, dt);
    }

    if (parameters::disk_feedback) {
	update_velocities_from_indirect_term(dt);
    }

    // as particles move independent of each other, we can integrate one after
    // one
    for (unsigned int i = 0; i < local_number_of_particles; ++i) {

	double fac, fac11;
	double atoli, rtoli, err, sk;
	double sqr;

	constexpr double beta = 0.04;
	constexpr double fac1 = 0.2;
	constexpr double fac2 = 10.0;
	constexpr double safe = 0.9;

	constexpr double expo1 = 0.2 - beta * 0.75;
	constexpr double facc1 = 1.0 / fac1;
	constexpr double facc2 = 1.0 / fac2;

	constexpr double e1 = 37.0 / 378.0 - 2825.0 / 27648.0;
	constexpr double e3 = 250.0 / 621.0 - 18575.0 / 48384.0;
	constexpr double e4 = 125.0 / 594.0 - 13525 / 55296.0;
	constexpr double e5 = -277.0 / 14336.0;
	constexpr double e6 = 512.0 / 1771.0 - 1.0 / 4.0;

	// initial preparations
	atoli = 1e-14;
	rtoli = 1e-12;
	bool last = false;

	bool reject = false;

	// basic integration step
	double dt_particle_temp = 0.0;

	while (!last) {
	    double timestep = particles[i].timestep;
	    const double &old_timestep = particles[i].timestep;
	    bool update_timestep = true;

	    if (dt_particle_temp + timestep * 1.01 > dt) {
		timestep = dt - dt_particle_temp;
		last = true;
		update_timestep = false;
	    }

	    double k1_r, k1_phi, k1_r_dot, k1_phi_dot;
	    double k2_r, k2_phi, k2_r_dot, k2_phi_dot;
	    double k3_r, k3_phi, k3_r_dot, k3_phi_dot;
	    double k4_r, k4_phi, k4_r_dot, k4_phi_dot;
	    double k5_r, k5_phi, k5_r_dot, k5_phi_dot;
	    double k6_r, k6_phi, k6_r_dot, k6_phi_dot;

	    double temp_ar, temp_aphi, temp_r, temp_phi, temp_r_dot,
		temp_phi_dot;

	    temp_r = particles[i].r;
	    temp_phi = particles[i].phi;
	    temp_r_dot = particles[i].r_dot;
	    temp_phi_dot = particles[i].phi_dot;

	    // Cash–Karp method
	    // (http://en.wikipedia.org/wiki/Cash%E2%80%93Karp_method) cartesian
	    // coordinates are written inside the polar coordinates
	    if (CartesianParticles) {
		calculate_accelerations_from_star_and_planets_cart(
		    temp_ar, temp_aphi, temp_r, temp_phi, data);
	    } else {
		calculate_accelerations_from_star_and_planets(
		    temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi,
		    temp_phi_dot, data);
	    }

	    // calculate k1
	    k1_r_dot = temp_ar;
	    k1_phi_dot = temp_aphi;
	    k1_r = temp_r_dot;
	    k1_phi = temp_phi_dot;

	    temp_r = particles[i].r + timestep * 0.2 * k1_r;
	    temp_phi = particles[i].phi + timestep * 0.2 * k1_phi;

	    temp_r_dot = particles[i].r_dot + timestep * 0.2 * k1_r_dot;
	    temp_phi_dot = particles[i].phi_dot + timestep * 0.2 * k1_phi_dot;

	    if (CartesianParticles) {
		calculate_accelerations_from_star_and_planets_cart(
		    temp_ar, temp_aphi, temp_r, temp_phi, data);
	    } else {
		calculate_accelerations_from_star_and_planets(
		    temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi,
		    temp_phi_dot, data);
	    }
	    // calculate k2
	    k2_r_dot = temp_ar;
	    k2_phi_dot = temp_aphi;
	    k2_r = temp_r_dot;
	    k2_phi = temp_phi_dot;

	    temp_r = particles[i].r + timestep * (0.075 * k1_r + 0.225 * k2_r);
	    temp_phi =
		particles[i].phi + timestep * (0.075 * k1_phi + 0.225 * k2_phi);

	    temp_r_dot = particles[i].r_dot +
			 timestep * (0.075 * k1_r_dot + 0.225 * k2_r_dot);
	    temp_phi_dot = particles[i].phi_dot +
			   timestep * (0.075 * k1_phi_dot + 0.225 * k2_phi_dot);

	    if (CartesianParticles) {
		calculate_accelerations_from_star_and_planets_cart(
		    temp_ar, temp_aphi, temp_r, temp_phi, data);
	    } else {
		calculate_accelerations_from_star_and_planets(
		    temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi,
		    temp_phi_dot, data);
	    }
	    // calculate k3
	    k3_r_dot = temp_ar;
	    k3_phi_dot = temp_aphi;
	    k3_r = temp_r_dot;
	    k3_phi = temp_phi_dot;

	    temp_r = particles[i].r +
		     timestep * (0.3 * k1_r - 0.9 * k2_r + 1.2 * k3_r);
	    temp_phi = particles[i].phi +
		       timestep * (0.3 * k1_phi - 0.9 * k2_phi + 1.2 * k3_phi);

	    temp_r_dot =
		particles[i].r_dot +
		timestep * (0.3 * k1_r_dot - 0.9 * k2_r_dot + 1.2 * k3_r_dot);
	    temp_phi_dot = particles[i].phi_dot +
			   timestep * (0.3 * k1_phi_dot - 0.9 * k2_phi_dot +
				       1.2 * k3_phi_dot);

	    if (CartesianParticles) {
		calculate_accelerations_from_star_and_planets_cart(
		    temp_ar, temp_aphi, temp_r, temp_phi, data);
	    } else {
		calculate_accelerations_from_star_and_planets(
		    temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi,
		    temp_phi_dot, data);
	    }
	    // calculate k4
	    k4_r_dot = temp_ar;
	    k4_phi_dot = temp_aphi;
	    k4_r = temp_r_dot;
	    k4_phi = temp_phi_dot;

	    temp_r = particles[i].r +
		     timestep * (-11.0 / 54.0 * k1_r + 2.5 * k2_r -
				 70.0 / 27.0 * k3_r + 35.0 / 27.0 * k4_r);
	    temp_phi = particles[i].phi +
		       timestep * (-11.0 / 54.0 * k1_phi + 2.5 * k2_phi -
				   70.0 / 27.0 * k3_phi + 35.0 / 27.0 * k4_phi);

	    temp_r_dot =
		particles[i].r_dot +
		timestep * (-11.0 / 54.0 * k1_r_dot + 2.5 * k2_r_dot -
			    70.0 / 27.0 * k3_r_dot + 35.0 / 27.0 * k4_r_dot);
	    temp_phi_dot =
		particles[i].phi_dot +
		timestep *
		    (-11.0 / 54.0 * k1_phi_dot + 2.5 * k2_phi_dot -
		     70.0 / 27.0 * k3_phi_dot + 35.0 / 27.0 * k4_phi_dot);

	    if (CartesianParticles) {
		calculate_accelerations_from_star_and_planets_cart(
		    temp_ar, temp_aphi, temp_r, temp_phi, data);
	    } else {
		calculate_accelerations_from_star_and_planets(
		    temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi,
		    temp_phi_dot, data);
	    }
	    // calculate k5
	    k5_r_dot = temp_ar;
	    k5_phi_dot = temp_aphi;
	    k5_r = temp_r_dot;
	    k5_phi = temp_phi_dot;

	    temp_r =
		particles[i].r +
		timestep * (1631.0 / 55296.0 * k1_r + 175.0 / 512.0 * k2_r +
			    575.0 / 13824.0 * k3_r + 44275.0 / 110592.0 * k4_r +
			    253.0 / 4096.0 * k5_r);
	    temp_phi =
		particles[i].phi +
		timestep *
		    (1631.0 / 55296.0 * k1_phi + 175.0 / 512.0 * k2_phi +
		     575.0 / 13824.0 * k3_phi + 44275.0 / 110592.0 * k4_phi +
		     253.0 / 4096.0 * k5_phi);

	    temp_r_dot =
		particles[i].r_dot +
		timestep *
		    (1631.0 / 55296.0 * k1_r_dot + 175.0 / 512.0 * k2_r_dot +
		     575.0 / 13824.0 * k3_r_dot +
		     44275.0 / 110592.0 * k4_r_dot + 253.0 / 4096.0 * k5_r_dot);
	    temp_phi_dot = particles[i].phi_dot +
			   timestep * (1631.0 / 55296.0 * k1_phi_dot +
				       175.0 / 512.0 * k2_phi_dot +
				       575.0 / 13824.0 * k3_phi_dot +
				       44275.0 / 110592.0 * k4_phi_dot +
				       253.0 / 4096.0 * k5_phi_dot);

	    if (CartesianParticles) {
		calculate_accelerations_from_star_and_planets_cart(
		    temp_ar, temp_aphi, temp_r, temp_phi, data);
	    } else {
		calculate_accelerations_from_star_and_planets(
		    temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi,
		    temp_phi_dot, data);
	    }
	    // calculate k6
	    k6_r_dot = temp_ar;
	    k6_phi_dot = temp_aphi;
	    k6_r = temp_r_dot;
	    k6_phi = temp_phi_dot;

	    // estimate new positions and velocities
	    const double r_new =
		particles[i].r +
		timestep * (37.0 / 378.0 * k1_r + 250.0 / 621.0 * k3_r +
			    125.0 / 594.0 * k4_r + 512.0 / 1771.0 * k6_r);
	    double phi_new =
		particles[i].phi +
		timestep * (37.0 / 378.0 * k1_phi + 250.0 / 621.0 * k3_phi +
			    125.0 / 594.0 * k4_phi + 512.0 / 1771.0 * k6_phi);

	    if (!CartesianParticles) {
		check_angle(phi_new);
	    }

	    const double r_ddot_new =
		(37.0 / 378.0 * k1_r_dot + 250.0 / 621.0 * k3_r_dot +
		 125.0 / 594.0 * k4_r_dot + 512.0 / 1771.0 * k6_r_dot);
	    const double phi_ddot_new =
		(37.0 / 378.0 * k1_phi_dot + 250.0 / 621.0 * k3_phi_dot +
		 125.0 / 594.0 * k4_phi_dot + 512.0 / 1771.0 * k6_phi_dot);

	    const double r_dot_new = particles[i].r_dot + timestep * r_ddot_new;
	    const double phi_dot_new =
		particles[i].phi_dot + timestep * phi_ddot_new;

	    // error estimates, reuse k2 variables
	    k2_r = timestep *
		   (e1 * k1_r + e3 * k3_r + e4 * k4_r + e5 * k5_r + e6 * k6_r);
	    k2_phi = timestep * (e1 * k1_phi + e3 * k3_phi + e4 * k4_phi +
				 e5 * k5_phi + e6 * k6_phi);
	    k2_r_dot =
		timestep * (e1 * k1_r_dot + e3 * k3_r_dot + e4 * k4_r_dot +
			    e5 * k5_r_dot + e6 * k6_r_dot);
	    k2_phi_dot = timestep *
			 (e1 * k1_phi_dot + e3 * k3_phi_dot + e4 * k4_phi_dot +
			  e5 * k5_phi_dot + e6 * k6_phi_dot);

	    // error estimation
	    err = 0.0;

	    sk = atoli + rtoli * fmax(fabs(particles[i].r), fabs(r_new));
	    sqr = k2_r / sk;
	    err += sqr * sqr;

	    sk = atoli + rtoli * fmax(fabs(particles[i].phi), fabs(phi_new));
	    sqr = k2_phi / sk;
	    err += sqr * sqr;

	    sk =
		atoli + rtoli * fmax(fabs(particles[i].r_dot), fabs(r_dot_new));
	    sqr = k2_r_dot / sk;
	    err += sqr * sqr;

	    sk = atoli +
		 rtoli * fmax(fabs(particles[i].phi_dot), fabs(phi_dot_new));
	    sqr = k2_phi_dot / sk;
	    err += sqr * sqr;

	    err = sqrt(err / (double)(4));
	    // computation of hnew
	    fac11 = pow(err, expo1);
	    // Lund-stabilization
	    fac = fac11 / pow(particles[i].facold, beta);
	    // we require fac1 <= hnew/h <= fac2
	    fac = fmax(facc2, fmin(facc1, fac / safe));

	    // if timestep is reduced to reach full dt, don't make it larger
	    if (!update_timestep)
		fac = fmax(fac, 1.0);

	    particles[i].timestep = old_timestep / fac;

	    if (err <= 1.0) {
		// step accepted
		particles[i].facold = fmax(err, 1.0e-4);
		dt_particle_temp += timestep;

		// update position & velocity
		particles[i].r = r_new;
		particles[i].phi = phi_new;

		particles[i].r_dot = r_dot_new;
		particles[i].phi_dot = phi_dot_new;

		particles[i].r_ddot = r_ddot_new;
		particles[i].phi_ddot = phi_ddot_new;

		if (reject) // last try was rejected
		{
		    particles[i].timestep =
			fmin(fabs(particles[i].timestep), fabs(old_timestep));
		}
		reject = false;
	    } else // timestep rejected
	    {
		particles[i].timestep =
		    old_timestep / fmin(facc1, fac11 / safe);
		reject = true;
		last = false;
	    }
	}
    }
    move();
}

void integrate_explicit(t_data &data, const double dt)
{

    if (CartesianParticles) {
	// disk gravity on particles is inside gas_drag function
	if (parameters::particle_gas_drag_enabled)
	    update_velocities_from_gas_drag_cart(data, dt);
    } else {
	// disk gravity on particles is inside gas_drag function
	if (parameters::particle_gas_drag_enabled)
	    update_velocities_from_gas_drag(data, dt);
    } // else {
    // disk gravity on particles is inside gas_drag function
    if (parameters::disk_feedback) {
	update_velocities_from_indirect_term(dt);
    }

    // as particles move independent of each other, we can integrate one after
    // one
    for (unsigned int i = 0; i < local_number_of_particles; ++i) {

	double k1_r, k1_phi, k1_r_dot, k1_phi_dot;
	double k2_r, k2_phi, k2_r_dot, k2_phi_dot;
	double k3_r, k3_phi, k3_r_dot, k3_phi_dot;
	double k4_r, k4_phi, k4_r_dot, k4_phi_dot;
	double k5_r, k5_phi, k5_r_dot, k5_phi_dot;
	double k6_r, k6_phi, k6_r_dot, k6_phi_dot;

	double temp_ar, temp_aphi, temp_r, temp_phi, temp_r_dot, temp_phi_dot;

	temp_r = particles[i].r;
	temp_phi = particles[i].phi;
	temp_r_dot = particles[i].r_dot;
	temp_phi_dot = particles[i].phi_dot;

	// Cash–Karp method
	// (http://en.wikipedia.org/wiki/Cash%E2%80%93Karp_method) cartesian
	// coordinates are written inside the polar coordinates
	if (CartesianParticles) {
	    calculate_accelerations_from_star_and_planets_cart(
		temp_ar, temp_aphi, temp_r, temp_phi, data);
	} else {
	    calculate_accelerations_from_star_and_planets(
		temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi, temp_phi_dot,
		data);
	}

	// calculate k1
	k1_r_dot = temp_ar;
	k1_phi_dot = temp_aphi;
	k1_r = temp_r_dot;
	k1_phi = temp_phi_dot;

	temp_r = particles[i].r + dt * 0.2 * k1_r;
	temp_phi = particles[i].phi + dt * 0.2 * k1_phi;

	temp_r_dot = particles[i].r_dot + dt * 0.2 * k1_r_dot;
	temp_phi_dot = particles[i].phi_dot + dt * 0.2 * k1_phi_dot;

	if (CartesianParticles) {
	    calculate_accelerations_from_star_and_planets_cart(
		temp_ar, temp_aphi, temp_r, temp_phi, data);
	} else {
	    calculate_accelerations_from_star_and_planets(
		temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi, temp_phi_dot,
		data);
	}
	// calculate k2
	k2_r_dot = temp_ar;
	k2_phi_dot = temp_aphi;
	k2_r = temp_r_dot;
	k2_phi = temp_phi_dot;

	temp_r = particles[i].r + dt * (0.075 * k1_r + 0.225 * k2_r);
	temp_phi = particles[i].phi + dt * (0.075 * k1_phi + 0.225 * k2_phi);

	temp_r_dot =
	    particles[i].r_dot + dt * (0.075 * k1_r_dot + 0.225 * k2_r_dot);
	temp_phi_dot = particles[i].phi_dot +
		       dt * (0.075 * k1_phi_dot + 0.225 * k2_phi_dot);

	if (CartesianParticles) {
	    calculate_accelerations_from_star_and_planets_cart(
		temp_ar, temp_aphi, temp_r, temp_phi, data);
	} else {
	    calculate_accelerations_from_star_and_planets(
		temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi, temp_phi_dot,
		data);
	}
	// calculate k3
	k3_r_dot = temp_ar;
	k3_phi_dot = temp_aphi;
	k3_r = temp_r_dot;
	k3_phi = temp_phi_dot;

	temp_r = particles[i].r + dt * (0.3 * k1_r - 0.9 * k2_r + 1.2 * k3_r);
	temp_phi = particles[i].phi +
		   dt * (0.3 * k1_phi - 0.9 * k2_phi + 1.2 * k3_phi);

	temp_r_dot = particles[i].r_dot +
		     dt * (0.3 * k1_r_dot - 0.9 * k2_r_dot + 1.2 * k3_r_dot);
	temp_phi_dot =
	    particles[i].phi_dot +
	    dt * (0.3 * k1_phi_dot - 0.9 * k2_phi_dot + 1.2 * k3_phi_dot);

	if (CartesianParticles) {
	    calculate_accelerations_from_star_and_planets_cart(
		temp_ar, temp_aphi, temp_r, temp_phi, data);
	} else {
	    calculate_accelerations_from_star_and_planets(
		temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi, temp_phi_dot,
		data);
	}
	// calculate k4
	k4_r_dot = temp_ar;
	k4_phi_dot = temp_aphi;
	k4_r = temp_r_dot;
	k4_phi = temp_phi_dot;

	temp_r =
	    particles[i].r + dt * (-11.0 / 54.0 * k1_r + 2.5 * k2_r -
				   70.0 / 27.0 * k3_r + 35.0 / 27.0 * k4_r);
	temp_phi = particles[i].phi +
		   dt * (-11.0 / 54.0 * k1_phi + 2.5 * k2_phi -
			 70.0 / 27.0 * k3_phi + 35.0 / 27.0 * k4_phi);

	temp_r_dot = particles[i].r_dot +
		     dt * (-11.0 / 54.0 * k1_r_dot + 2.5 * k2_r_dot -
			   70.0 / 27.0 * k3_r_dot + 35.0 / 27.0 * k4_r_dot);
	temp_phi_dot =
	    particles[i].phi_dot +
	    dt * (-11.0 / 54.0 * k1_phi_dot + 2.5 * k2_phi_dot -
		  70.0 / 27.0 * k3_phi_dot + 35.0 / 27.0 * k4_phi_dot);

	if (CartesianParticles) {
	    calculate_accelerations_from_star_and_planets_cart(
		temp_ar, temp_aphi, temp_r, temp_phi, data);
	} else {
	    calculate_accelerations_from_star_and_planets(
		temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi, temp_phi_dot,
		data);
	}
	// calculate k5
	k5_r_dot = temp_ar;
	k5_phi_dot = temp_aphi;
	k5_r = temp_r_dot;
	k5_phi = temp_phi_dot;

	temp_r = particles[i].r +
		 dt * (1631.0 / 55296.0 * k1_r + 175.0 / 512.0 * k2_r +
		       575.0 / 13824.0 * k3_r + 44275.0 / 110592.0 * k4_r +
		       253.0 / 4096.0 * k5_r);
	temp_phi = particles[i].phi +
		   dt * (1631.0 / 55296.0 * k1_phi + 175.0 / 512.0 * k2_phi +
			 575.0 / 13824.0 * k3_phi +
			 44275.0 / 110592.0 * k4_phi + 253.0 / 4096.0 * k5_phi);

	temp_r_dot =
	    particles[i].r_dot +
	    dt * (1631.0 / 55296.0 * k1_r_dot + 175.0 / 512.0 * k2_r_dot +
		  575.0 / 13824.0 * k3_r_dot + 44275.0 / 110592.0 * k4_r_dot +
		  253.0 / 4096.0 * k5_r_dot);
	temp_phi_dot =
	    particles[i].phi_dot +
	    dt *
		(1631.0 / 55296.0 * k1_phi_dot + 175.0 / 512.0 * k2_phi_dot +
		 575.0 / 13824.0 * k3_phi_dot +
		 44275.0 / 110592.0 * k4_phi_dot + 253.0 / 4096.0 * k5_phi_dot);

	if (CartesianParticles) {
	    calculate_accelerations_from_star_and_planets_cart(
		temp_ar, temp_aphi, temp_r, temp_phi, data);
	} else {
	    calculate_accelerations_from_star_and_planets(
		temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi, temp_phi_dot,
		data);
	}
	// calculate k6
	k6_r_dot = temp_ar;
	k6_phi_dot = temp_aphi;
	k6_r = temp_r_dot;
	k6_phi = temp_phi_dot;

	// update position & velocity
	particles[i].r += dt * (37.0 / 378.0 * k1_r + 250.0 / 621.0 * k3_r +
				125.0 / 594.0 * k4_r + 512.0 / 1771.0 * k6_r);
	particles[i].phi +=
	    dt * (37.0 / 378.0 * k1_phi + 250.0 / 621.0 * k3_phi +
		  125.0 / 594.0 * k4_phi + 512.0 / 1771.0 * k6_phi);

	if (!CartesianParticles) {
	    check_angle(particles[i].phi);
	}

	particles[i].r_ddot =
	    (37.0 / 378.0 * k1_r_dot + 250.0 / 621.0 * k3_r_dot +
	     125.0 / 594.0 * k4_r_dot + 512.0 / 1771.0 * k6_r_dot);
	particles[i].phi_ddot =
	    (37.0 / 378.0 * k1_phi_dot + 250.0 / 621.0 * k3_phi_dot +
	     125.0 / 594.0 * k4_phi_dot + 512.0 / 1771.0 * k6_phi_dot);

	particles[i].r_dot += dt * particles[i].r_ddot;
	particles[i].phi_dot += dt * particles[i].phi_ddot;
    }

    move();
}

void move(void)
{
    // remove escaped particles
    for (unsigned int i = 0; i < local_number_of_particles; ++i) {
	if (particles[i].get_squared_distance_to_star() >
	    parameters::particle_maximum_escape_radius_sq) {
	    particles[i] = particles[local_number_of_particles - 1];
	    local_number_of_particles--;
	    i--; // new particle with id 'i' must be checked too
	} else if (particles[i].get_squared_distance_to_star() <
		   parameters::particle_minimum_escape_radius_sq) {
	    particles[i] = particles[local_number_of_particles - 1];
	    local_number_of_particles--;
	    i--; // new particle with id 'i' must be checked too
	}
    }

    // we don't need any move if we're running single core
    if (CPU_Number == 1)
	return;

    // check if particles are still on correct node
    double local_r_min_squared = pow2(local_r_min);
    double local_r_max_squared = pow2(local_r_max);

    // move particles inwards starting from the outer most node

    // receive particles from outer node first
    if (CPU_Rank < CPU_Highest) {
	MPI_Status status;
	MPI_Probe(CPU_Next, 0, MPI_COMM_WORLD, &status);
	int number;
	MPI_Get_count(&status, mpi_particle, &number);

	// check if array is large enough for new particles
	if (particles_size < local_number_of_particles + number) {
	    particles_size = local_number_of_particles + number;
	    particles = (t_particle *)realloc(particles, sizeof(t_particle) *
							     particles_size);
	}

	MPI_Recv(particles + local_number_of_particles, number, mpi_particle,
		 CPU_Next, 0, MPI_COMM_WORLD, &status);

	local_number_of_particles += number;
    }

    // check if we need to send particles to inner node
    std::vector<int> inward_offset(local_number_of_particles);
    std::vector<int> inward_size(local_number_of_particles);
    int inward_count = 0;
    for (unsigned int i = 0; i < local_number_of_particles; ++i) {
	if ((CPU_Rank > 0) && (particles[i].get_squared_distance_to_star() <
			       local_r_min_squared)) {
	    inward_offset[inward_count] = i;
	    inward_size[inward_count] = 1;
	    inward_count++;
	}
    }

    // send particles to inner node
    if (CPU_Rank > 0) {
	MPI_Datatype inward_type;
	MPI_Type_indexed(inward_count, &inward_size[0], &inward_offset[0],
			 mpi_particle, &inward_type);
	MPI_Type_commit(&inward_type);
	MPI_Send(particles, 1, inward_type, CPU_Prev, 0, MPI_COMM_WORLD);
	MPI_Type_free(&inward_type);
    }

    // delete particles sent to inner node
    for (int i = 0; i < inward_count; ++i) {
	// swap particle with last particle
	particles[inward_offset[i]] = particles[local_number_of_particles - 1];

	// check if the particle we switched is also in inward_offset
	for (int j = i + 1; j < inward_count; ++j) {
	    if (inward_offset[j] == (int)local_number_of_particles - 1)
		inward_offset[j] = inward_offset[i];
	}

	local_number_of_particles--;
    }

    // receive particles from inner node first
    if (CPU_Rank > 0) {
	MPI_Status status;
	MPI_Probe(CPU_Prev, 0, MPI_COMM_WORLD, &status);
	int number;
	MPI_Get_count(&status, mpi_particle, &number);

	// check if array is large enough for new particles
	if (particles_size < local_number_of_particles + number) {
	    particles_size = local_number_of_particles + number;
	    particles = (t_particle *)realloc(particles, sizeof(t_particle) *
							     particles_size);
	}

	MPI_Recv(particles + local_number_of_particles, number, mpi_particle,
		 CPU_Prev, 0, MPI_COMM_WORLD, &status);

	local_number_of_particles += number;
    }

    // move particles outwards starting from the inner most node

    // check if we need to send particles to outer node
    std::vector<int> outward_offset(local_number_of_particles);
    std::vector<int> outward_size(local_number_of_particles);
    int outward_count = 0;
    for (unsigned int i = 0; i < local_number_of_particles; ++i) {
	if ((CPU_Rank < CPU_Highest) &&
	    (particles[i].get_squared_distance_to_star() >
	     local_r_max_squared)) {
	    outward_offset[outward_count] = i;
	    outward_size[outward_count] = 1;
	    outward_count++;
	}
    }

    // send particles to outer node
    if (CPU_Rank < CPU_Highest) {
	MPI_Datatype outward_type;
	MPI_Type_indexed(outward_count, &outward_size[0], &outward_offset[0],
			 mpi_particle, &outward_type);
	MPI_Type_commit(&outward_type);
	MPI_Send(particles, 1, outward_type, CPU_Next, 0, MPI_COMM_WORLD);
	MPI_Type_free(&outward_type);
    }

    // delete particles sent to inner node
    for (int i = 0; i < outward_count; ++i) {
	// swap particle with last particle
	particles[outward_offset[i]] = particles[local_number_of_particles - 1];

	// check if the particle we switched is also in outward_offset
	for (int j = i + 1; j < outward_count; ++j) {
	    if (outward_offset[j] == (int)local_number_of_particles - 1)
		outward_offset[j] = outward_offset[i];
	}

	local_number_of_particles--;
    }

    // update global_number_of_particles
    MPI_Allreduce(&local_number_of_particles, &global_number_of_particles, 1,
		  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

void write(unsigned int timestep)
{
    MPI_File fh;
    MPI_Status status;
    int error, error_class, error_length;
    char *filename, error_string[MPI_MAX_ERROR_STRING + 1];

    if (asprintf(&filename, "%s/particles%i.dat", OUTPUTDIR, timestep) < 0) {
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

    // get number of local particles from all nodes to compute correct offsets
    std::vector<unsigned int> nodes_number_of_particles(CPU_Number);
    MPI_Allgather(&local_number_of_particles, 1, MPI_UNSIGNED,
		  &nodes_number_of_particles[0], 1, MPI_UNSIGNED,
		  MPI_COMM_WORLD);

    // compute local offset
    unsigned int local_offset = 0;
    for (int cpu = 0; cpu < CPU_Rank; ++cpu) {
	local_offset += nodes_number_of_particles[cpu];
    }

    MPI_File_set_view(fh, 0, mpi_particle, mpi_particle,
		      const_cast<char *>("native"), MPI_INFO_NULL);
    MPI_File_seek(fh, local_offset, MPI_SEEK_SET);
    MPI_File_write(fh, particles, local_number_of_particles, mpi_particle,
		   &status);

    // close file
    MPI_File_close(&fh);
}

void rotate(double Omega, double dt)
{
    for (unsigned int i = 0; i < local_number_of_particles; ++i) {
	// rotate positions
	if (CartesianParticles) {
	    const double x_old = particles[i].r;
	    const double y_old = particles[i].phi;

	    const double vx_old = particles[i].r_dot;
	    const double vy_old = particles[i].phi_dot;

	    double &x = particles[i].r;
	    double &y = particles[i].phi;

	    double &vx = particles[i].r_dot;
	    double &vy = particles[i].phi_dot;

	    // rotate positions
	    double angle = Omega * dt;
	    x = x_old * cos(angle) + y_old * sin(angle);
	    y = -x_old * sin(angle) + y_old * cos(angle);
	    // rotate velocities
	    vx = vx_old * cos(angle) + vy_old * sin(angle);
	    vy = -vx_old * sin(angle) + vy_old * cos(angle);
	} else {

	    particles[i].phi -= Omega * dt;
	    check_angle(particles[i].phi);
	}
    }
}

} // namespace particles
