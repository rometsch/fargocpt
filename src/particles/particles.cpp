/**
	\file particles.cpp
	\author Tobias Mueller <Tobias_Mueller@twam.info>
*/

#include "particles.h"
#include "dust_diffusion.h"
#include "../Force.h"
#include "../LowTasks.h"
#include "../Pframeforce.h"
#include "../SourceEuler.h"
#include "../commbound.h"
#include "../constants.h"
#include "../find_cell_id.h"
#include "../global.h"
#include "../logging.h"
#include "../mpi_utils.h"
#include "../parameters.h"
#include "../selfgravity.h"
#include "../util.h"
#include "../output.h"
#include "../Theo.h"
#include "../frame_of_reference.h"
#include "../cfl.h"
#include "../simulation.h"
#include <cstring>
#include <cmath>
#include <mpi.h>
#include <random>
#include <sstream>
#include <stdlib.h>
#include <vector>


namespace particles
{

/// local particle storage
std::vector<t_particle> particles;

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
static double power_law_distribution(double x, double n)
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
	return std::exp(std::log(rmax / rmin) * x) * rmin;
    } else {
	return std::exp(
	    std::log(x * (-std::pow(rmin, n + 1) + std::pow(rmax, n + 1)) +
		     std::pow(rmin, n + 1)) /
	    (n + 1));
    }
}

static double check_angle(double &phi)
{
    if (phi >= 2.0 * M_PI) {
	return phi -= 2.0 * M_PI;
    } else {
	if (phi < 0.0) {
	    return phi += 2.0 * M_PI;
	} else {
	    return phi;
	}
    }
}

static double corret_v_gas_azimuthal_omega_frame(const double v_az_gas,
						 const double r_particle)
{
    const double v_az_corrected = v_az_gas + r_particle * refframe::OmegaFrame;
    return v_az_corrected;
}

static void transform_cart_to_cyl(const double *const cart, double *const cyl,
				  const double r, const double phi)
{
    (void)r;
    cyl[0] = cart[0] * std::cos(phi);
    cyl[0] += cart[1] * std::sin(phi);

    cyl[1] = -cart[0] * std::sin(phi);
    cyl[1] += cart[1] * std::cos(phi);
}

static void
find_nearest(unsigned int &n_radial_a_minus, unsigned int &n_radial_a_plus,
	     unsigned int &n_azimuthal_a_minus,
	     unsigned int &n_azimuthal_a_plus, unsigned int &n_radial_b_minus,
	     unsigned int &n_radial_b_plus, unsigned int &n_azimuthal_b_minus,
	     unsigned int &n_azimuthal_b_plus, const double r, const double phi)
{
    n_radial_a_minus = get_rinf_id(r);
    n_radial_a_plus = n_radial_a_minus + 1;

    n_radial_b_minus = get_rmed_id(r);
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
    double dphi = 2.0 * M_PI / (double)NAzimuthal;

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
	if (phi < M_PI) {
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

    if ((phim > 2.0 * M_PI) || (phip > 2.0 * M_PI)) {
	phim -= M_PI;
	phip -= M_PI;
	phi -= M_PI;
    }

    double Qm = ((rm * phip - rm * phi) * Qmm + (rm * phi - rm * phim) * Qmp) /
		(rm * phip - rm * phim);
    double Qp = ((rp * phip - rp * phi) * Qpm + (rp * phi - rp * phim) * Qpp) /
		(rp * phip - rp * phim);

    return ((rp - r) * Qm + (r - rm) * Qp) / (rp - rm);
}

static double get_particle_eccentricity_polar(int particle_index)
{
    // used for debugging the integrators
    const int i = particle_index;

    // Runge-Lenz vector A = (p x L) - m * G * m * hydro_center_mass * r/|r|
    const double kappa = constants::G * hydro_center_mass * particles[i].mass;
    const double reduced_mass = hydro_center_mass * particles[i].mass /
				(hydro_center_mass + particles[i].mass);
    const double A_r = reduced_mass * std::pow(particles[i].r, 3) *
			   std::pow(particles[i].phi_dot, 2) -
		       kappa;
    const double A_phi = reduced_mass * std::pow(particles[i].r, 2) *
			 particles[i].r_dot * particles[i].phi_dot;
    const double A = std::sqrt(A_r * A_r + A_phi * A_phi);
    const double eccentricity = A / kappa;

    return eccentricity;
}

static double get_particle_eccentricity_cart(int particle_index)
{
    // used for debugging the integrators
    const int i = particle_index;

    // Runge-Lenz vector A = (p x L) - m * G * m * hydro_center_mass * r/|r|
    const double m = hydro_center_mass + particles[i].mass;

    const double x = particles[i].r;
    const double y = particles[i].phi;

    const double d = std::sqrt(x * x + y * y);

    const double vx = particles[i].r_dot;
    const double vy = particles[i].phi_dot;
    const double Ax = x * vy * vy - y * vx * vy - constants::G * m * x / d;
    const double Ay = y * vx * vx - x * vx * vy - constants::G * m * y / d;
    const double e = std::sqrt(Ax * Ax + Ay * Ay) / constants::G / m;

    return e;
}

[[maybe_unused]] static double get_particle_eccentricity(int particle_index)
{
    if (parameters::CartesianParticles) {
	return get_particle_eccentricity_cart(particle_index);
    } else {
	return get_particle_eccentricity_polar(particle_index);
    }
}

/**
 * @brief init_particle_timestep init particle timestep for adaptive explicit
 * integrator
 * @param data
 */
static void init_particle_timestep(t_data &data)
{

	#pragma omp parallel for
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

	const double rsmooth = 0.05;

	// Cash–Karp method
	// (http://en.wikipedia.org/wiki/Cash%E2%80%93Karp_method) cartesian
	// coordinates are written inside the polar coordinates
	if (parameters::CartesianParticles) {
	    calculate_accelerations_from_star_and_planets_cart(
		temp_ar, temp_aphi, temp_r, temp_phi, rsmooth, data);
	} else {
	    calculate_accelerations_from_star_and_planets(
		temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi, temp_phi_dot,
		rsmooth, data);
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
	    particles[i].timestep = std::sqrt(dny / dnf) * 0.01;

	// perform an explicit Euler step
	temp_r = particles[i].r + particles[i].timestep * k1_r;
	temp_phi = particles[i].phi + particles[i].timestep * k1_phi;

	temp_r_dot = particles[i].r_dot + particles[i].timestep * k1_r_dot;
	temp_phi_dot =
	    particles[i].phi_dot + particles[i].timestep * k1_phi_dot;

	if (parameters::CartesianParticles) {
	    calculate_accelerations_from_star_and_planets_cart(
		temp_ar, temp_aphi, temp_r, temp_phi, rsmooth, data);
	} else {
	    calculate_accelerations_from_star_and_planets(
		temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi, temp_phi_dot,
		rsmooth, data);
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

	der2 = std::sqrt(der2) / particles[i].timestep;

	// step size is computed such that h**iord * max_d(norm(f0),norm(der2))
	// = 0.01
	der12 = std::fmax(std::fabs(der2), std::sqrt(dnf));
	if (der12 <= 1.0E-15)
	    h1 = std::fmax(1.0E-6, std::fabs(particles[i].timestep) * 1.0E-3);
	else
	    h1 = std::pow(0.01 / der12, 1.0 / (double)iord);
	particles[i].timestep = fmin(100.0 * particles[i].timestep, h1);
    }
}

static void correct_for_self_gravity(const unsigned int i)
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
	find_nearest(n_radial_a_minus, n_radial_a_plus, n_azimuthal_a_minus,
		     n_azimuthal_a_plus, n_radial_b_minus, n_radial_b_plus,
		     n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);

	sg_radial = interpolate_bilinear_sg(
	    selfgravity::g_radial, n_radial_a_minus, n_radial_a_plus,
	    n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
    }

    double v = sqrt(constants::G * (hydro_center_mass + particles[i].mass) /
			semi_major_axis -
		    sg_radial * semi_major_axis) *
	       std::sqrt((1.0 - eccentricity) / (1.0 + eccentricity));

    // only need to update velocities, set accelerations to 0 anyway
    if (parameters::CartesianParticles) {
	// Beware: cartesian particles still use the names of polar coordinates
	// vx = r_dot
	// vy = phi_dot

	particles[i].r_dot = -v * std::sin(phi);
	particles[i].phi_dot = v * std::cos(phi);
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
    double semi_major_axis =
	power_law_distribution(dis_one(generator), parameters::particle_slope);
    double phi = dis_twoPi(generator);
    double eccentricity = dis_eccentricity(generator);

    /*
    // debug setup
////////////////////////////////////////////////////////////// const int
num_particles_per_ring = 10; int global_id = i + id_offset;

double rmin = parameters::particle_minimum_radius;
double rmax = parameters::particle_maximum_radius;

// make sure particles are not initialized on an orbit that moves out of
    //particle-bounds
    if(RMAX <
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
    rmin+0.3;
    eccentricity = 0.0;
    ///////////////////////////////////////////////////////////////////////////
    */

    particles[i].radius = parameters::particle_radius;

    const unsigned int particle_type = i % parameters::particle_species_number;
    particles[i].radius *=
	std::pow(parameters::particle_radius_increase_factor, particle_type);

    double volume = 4.0 / 3.0 * M_PI * std::pow(particles[i].radius, 3);

    particles[i].mass = volume * parameters::particle_density;

    double r = semi_major_axis * (1.0 + eccentricity);
	// TODO: adjust v to include smoothing!
	// Smoothing reduces the gravitational pull, so the velocity needs to be reduced.
	// Otherwise, the orbit of the particle gets eccentric.
    double v =
	std::sqrt(constants::G * (hydro_center_mass + particles[i].mass) /
		  semi_major_axis) *
	std::sqrt((1.0 - eccentricity) / (1.0 + eccentricity));

    if (parameters::CartesianParticles) {
	// Beware: cartesian particles still use the names of polar coordinates
	// x = r
	// y = phi
	// vx = r_dot
	// vy = phi_dot
	particles[i].r = r * std::cos(phi);
	particles[i].phi = r * std::sin(phi);

	particles[i].r_dot = -v * std::sin(phi);
	particles[i].phi_dot = v * std::cos(phi);
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

/**
	computes density rho
*/
void compute_rho(t_data &data, const double current_time)
{
	compute_scale_height(data, current_time);

	const unsigned int Nr = data[t_data::RHO].get_size_radial();
	const unsigned int Nphi = data[t_data::RHO].get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const double H = data[t_data::SCALE_HEIGHT](nr, naz);
		data[t_data::RHO](nr, naz) =
		data[t_data::SIGMA](nr, naz) /
		(parameters::density_factor * H);
	}
    }
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
    particles.resize(particles_size);

    const unsigned int seed = parameters::random_seed;
    logging::print(LOG_DEBUG "random generator seed: %u\n", seed);

    // random generator and distributions
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> dis_one(0.0, 1.0);
    std::uniform_real_distribution<double> dis_twoPi(0.0, 2.0 * M_PI);

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
	const int mpi_particle_count = 12;
    int mpi_particle_lengths[mpi_particle_count] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    MPI_Aint mpi_particle_offsets[mpi_particle_count];
    MPI_Datatype mpi_particle_types[mpi_particle_count] = {MPI_UNSIGNED, MPI_DOUBLE, MPI_DOUBLE,
					   MPI_DOUBLE,	 MPI_DOUBLE, MPI_DOUBLE,
					   MPI_DOUBLE,	 MPI_DOUBLE, MPI_DOUBLE,
					   MPI_DOUBLE,   MPI_DOUBLE, MPI_DOUBLE};

    // calculate offsets
    MPI_Aint base;
    MPI_Get_address(&particles[0], &base);
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
	MPI_Get_address(&particles[0].stokes, mpi_particle_offsets + 11);


    for (int i = 0; i < mpi_particle_count; ++i) {
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
	#pragma omp parallel for
	for (unsigned int i = 0; i < local_number_of_particles; ++i) {
	    correct_for_self_gravity(i);
	}
    }

	compute_rho(data, sim::PhysicalTime);
	compute_temperature(data);
    check_tstop(data);

	if (parameters::particle_integrator == parameters::integrator_adaptive) {
		init_particle_timestep(data);
	}

	if (parameters::particle_dust_diffusion) {
		dust_diffusion::init(data);
	}

}

void restart()
{

    logging::print_master(
	LOG_WARNING
	"Beware: when restarting particles, the user is responsible that the loaded particle file is written with the same coordinate system as the simulation is running on!\n\n");

    FILE *fd;
    const std::string filename = output::snapshot_dir + "/particles.dat";

    fd = fopen(filename.c_str(), "r");
    if (fd == nullptr) {
	logging::print_master(
	    LOG_INFO
	    "Can't find file particles.dat (%s). Using generated particles.\n",
	    filename.c_str());
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
	particles.resize(particles_size);
    }
    if (num_particles_in_file > parameters::number_of_particles) {
	logging::print_master(
	    LOG_WARNING
	    "Warning: particle File contains %d particles but only %d particles can be loaded into Simulation!\n\n",
	    num_particles_in_file, parameters::number_of_particles);
    }

    MPI_File fh;
    MPI_Status status;

    // try to open file

    mpi_error_check_file_read(
	MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY,
		      MPI_INFO_NULL, &fh),
	filename.c_str());

    logging::print_master(LOG_INFO "Reading file '%s' with %u bytes.\n",
			  filename.c_str(), size);

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
    MPI_File_read_all(fh, &particles[0], local_number_of_particles, mpi_particle,
		      &status);

    // close file
    MPI_File_close(&fh);

    // particles might be loaded to the wrong node, so move them to the correct
    // one
    for (int i = 0; i < CPU_Number; ++i)
	move();
}

/*
	Calculate the smoothing length from the dust scale height using
	Dubrulle et al. (1995) Eq. 39 H_d = H_g * sqrt(alpha/(alpha + St))
*/
double calculate_dust_smoothing(const double r, const double phi,
				const double stokes, t_data &data)
{

    // find indices of current cell
    const unsigned int n_rad = get_rmed_id(r);
    const unsigned int n_azi = clamp_phi_id_to_grid(get_inf_azimuthal_id(phi));
    // get gas scale height
    const double h_gas = data[t_data::ASPECTRATIO](n_rad, n_azi);
    const double alpha = parameters::ALPHAVISCOSITY;
    // compute dust scale height
    const double h_dust = h_gas * std::sqrt(alpha / (alpha + stokes));

    const double rsmooth = h_dust * r * parameters::thickness_smoothing;

    return rsmooth;
}

void calculate_accelerations_from_star_and_planets(
    double &ar, double &aphi, const double r, const double r_dot,
    const double phi, const double phi_dot, const double rsmooth, t_data &data)
{
    (void)r_dot;
    (void)phi_dot;

    const double epsilon = rsmooth;
    const double epsilon_sq = epsilon * epsilon;

    ar = r * phi_dot * phi_dot; // Centrifugal force
    aphi = -2.0 * r_dot / r * phi_dot;

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

	const double distance_to_planet =
	    std::sqrt(r * r + r_planet * r_planet -
		      2 * r * r_planet * cos_delta_phi + epsilon_sq);
	const double distance_to_planet_pow2_smoothed =
	    distance_to_planet * distance_to_planet;

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
							const double rsmooth,
							t_data &data)
{
    ax = 0;
    ay = 0;

    const double epsilon = rsmooth;
    const double epsilon_sq = epsilon * epsilon;

    // planets
    for (unsigned int k = 0;
	 k < data.get_planetary_system().get_number_of_planets(); ++k) {
	t_planet &planet = data.get_planetary_system().get_planet(k);
	double r2 =
	    std::pow(planet.get_x() - x, 2) + std::pow(planet.get_y() - y, 2);
	r2 += epsilon_sq;
	const double r = std::sqrt(r2);
	const double factor = constants::G * planet.get_mass() / (r * r2);

	ax += factor * (planet.get_x() - x);
	ay += factor * (planet.get_y() - y);
    }
}

void calculate_derivitives_from_star_and_planets(double &grav_r_ddot,
						 double &minus_grav_l_dot,
						 const double r,
						 const double phi, const double rsmooth, t_data &data)
{
    const double epsilon = rsmooth;
    const double epsilon_sq = epsilon * epsilon;

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

	const double distance_to_planet = std::sqrt(
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

static void calculate_derivitives_from_star_and_planets_in_cart(
    double &grav_r_ddot, double &minus_grav_l_dot, const double r,
    const double phi, const double rsmooth, t_data &data)
{
    const double epsilon = rsmooth;
    const double epsilon_sq = epsilon * epsilon;

    double acart[2];
    double acyl[2];

    acart[0] = 0.0;
    acart[1] = 0.0;

    const double x = r * std::cos(phi);
    const double y = r * std::sin(phi);

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
	dist2 += epsilon_sq;
	const double dist = std::sqrt(dist2);

	// direct term
	acart[0] += -constants::G * planet_mass * x_dist / (dist * dist2);
	acart[1] += -constants::G * planet_mass * y_dist / (dist * dist2);
    }

    transform_cart_to_cyl(acart, acyl, r, phi);
    grav_r_ddot = acyl[0];
    minus_grav_l_dot = acyl[1] * r;
}

static double
interpolate_bilinear(t_polargrid &quantity, bool radial_a_grid,
		     bool azimuthal_a_grid, unsigned int n_radial_minus,
		     unsigned int n_radial_plus, unsigned int n_azimuthal_minus,
		     unsigned int n_azimuthal_plus, double r, double phi)
{
    const double dphi = 2.0 * M_PI / (double)quantity.get_size_azimuthal();
    const unsigned int last_phi_index = quantity.get_max_azimuthal();

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
	if (n_azimuthal_minus == last_phi_index) {
	    if (phi < M_PI) {
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


static double calc_tstop_only(const double size, const double rho, const double vrel, const double temperature) {

    // From Giovanni Picogna, used in
    // https://www.aanda.org/articles/aa/pdf/2018/08/aa32523-17.pdf propably
    // adapted from https://www.aanda.org/articles/aa/pdf/2003/07/aah3912.pdf

    const double m0 = parameters::MU * constants::m_u.get_code_value();
    const double vthermal = std::sqrt(8.0 * constants::k_B.get_code_value() *
				      temperature / (M_PI * m0));

    if (vthermal < 1.e-20)
	die("Zero VT %e\n", vthermal);
    if (vthermal > 1.e20)
	die("Zero VT1 %e\n", vthermal);

    double sigma = M_PI * std::pow(1.5e-8 / units::length.get_cgs_factor(), 2);
    double nu = 1.0 / 3.0 * m0 * vthermal / sigma;
    
	if (nu < 1.e-20) {
		die("Zero nu %e\n", nu);
    } else if (nu > 1.e20) {
		die("Zero nu1 %e\n", nu);
    }

	double l = 4.72e-9 / rho;
    
	if (l < 1.e-20) {
		die("Zero l %e\n", l);
	} else if (l > 1.e20) {
		die("Zero l1 %e\n", l);
	}
    
	double c_s = vthermal * sqrt(M_PI / 8.0);
    
	if (c_s < 1.e-20) {
		die("Zero cs %e\n", c_s);
	} else if (c_s > 1.e20) {
		die("Zero cs1 %e\n", c_s);
	}
    
	double Kn = 0.5 * l / size;
    double Ma = vrel / c_s;

    if (Ma < 1.e-20) {
		die("Zero Ma %e\n", Ma);
	} else if (Ma > 1.e20) {
		die("Zero Ma1 %e\n", Ma);
	}
    
	double Re = 2.0 * size * rho * vrel / nu;
    
	double CdE = 2.0 * sqrt(Ma * Ma + 128.0 / 9.0 / M_PI);
	
	if (CdE < 1.e-20) {
		die("Zero CdE %e\n", CdE);
	} else if (CdE > 1.e20) {
		die("Zero CdE1 %e\n", CdE);
	}
    
	double CdS = 0.0;
    
	if (Re <= 1.e-3) {
		CdS = 24.0 * nu / (2.0 * size * rho * c_s) +
	      3.6 / c_s * pow(vrel, 0.687) * pow(2.0 * size * rho / nu, -0.313);
    } else if (Re <= 500.0) {
		CdS = 24.0 * Ma / Re + 3.6 * Ma * pow(Re, -0.313);
    } else if (Re <= 1500.0) {
		CdS = Ma * 9.5e-5 * pow(Re, 1.397);
    } else {
		CdS = Ma * 2.61;
    }
    
	if (CdS < 1.e-30) {
		die("Zero CdS %e\n", CdS);
	} else if (CdS > 1.e30) {
		die("Zero CdS1 %e\n", CdS);
	}

    double Cd = (9.0 * Kn * Kn * CdE + CdS) / (3.0 * Kn + 1.0) / (3.0 * Kn + 1.0);
    
	if (Cd < 1.e-20) {
		die("Zero Cd %e\n", Cd);
	} else if (Cd > 1.e20) {
		die("Zero Cd1 %e\n", Cd);
	}
    
	const double pdens = parameters::particle_density;
	const double tstop = 4.0 * l * pdens / (3.0 * rho * Cd * c_s * Kn);
    // printf("Cd = %.3e\n", Cd);
    // printf("Kn = %.3e\n", Kn);
    // printf("CdE = %.3e\n", CdE);
    // printf("CdS = %.3e\n", CdS);
    // printf("Ma = %.3e\n", Ma);
    // printf("particle density = %.3e\n", parameters::particle_density);

    // printf("tstop = %.3e, tstop_old = %.3e\n", tstop, tstop_old);
    return tstop;
}

static void interpolate_quantities(t_data& data, const double r, const double phi, double &rho, double &temperature, double &vg_radial, double &vg_azimuthal) {
	
	unsigned int n_radial_a_minus = 0, n_radial_a_plus = 1,
		 n_radial_b_minus = 0, n_radial_b_plus = 1;
    unsigned int n_azimuthal_a_minus = 0, n_azimuthal_a_plus = 0,
		 n_azimuthal_b_minus = 0, n_azimuthal_b_plus = 0;
    find_nearest(n_radial_a_minus, n_radial_a_plus, n_azimuthal_a_minus,
		 n_azimuthal_a_plus, n_radial_b_minus, n_radial_b_plus,
		 n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);

    // calculate gas quantities at the particle location
    rho = interpolate_bilinear(
	data[t_data::RHO], false, false, n_radial_b_minus, n_radial_b_plus,
	n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
    
	temperature = interpolate_bilinear(
	data[t_data::TEMPERATURE], false, false, n_radial_b_minus,
	n_radial_b_plus, n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
    
	vg_radial = interpolate_bilinear(
	data[t_data::V_RADIAL], true, false, n_radial_a_minus, n_radial_a_plus,
	n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
    
	const double vg_azimuthal_temp = interpolate_bilinear(
	data[t_data::V_AZIMUTHAL], false, true, n_radial_b_minus,
	n_radial_b_plus, n_azimuthal_a_minus, n_azimuthal_a_plus, r, phi);
    
	vg_azimuthal = corret_v_gas_azimuthal_omega_frame(vg_azimuthal_temp, r);
}

static void calculate_tstop(const double r, const double phi,
			    const double r_dot, const double phi_dot,
			    t_data &data, const double radius,
			    double &minus_r_dotel_r, double &minus_l_rel,
			    double &tstop)
{

    double rho, temperature, vg_radial, vg_azimuthal;
	interpolate_quantities(data, r, phi, rho, temperature, vg_radial, vg_azimuthal);

    // calculate relative velocities
    minus_r_dotel_r = vg_radial - r_dot;
    const double vrel_phi = vg_azimuthal - phi_dot * r;

    minus_l_rel = r * vrel_phi;

    const double vrel =
	std::sqrt(std::pow(minus_r_dotel_r, 2) + std::pow(vrel_phi, 2));

	tstop = calc_tstop_only(radius, rho, vrel, temperature);
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
    r_tmp = fmin(r_tmp, local_r_max);

	double rho, temperature, vg_radial, vg_azimuthal;
	interpolate_quantities(data, r_tmp, phi, rho, temperature, vg_radial, vg_azimuthal);
    
    // calculate relative velocities
    minus_r_dotel_r = vg_radial - r_dot;
    const double vrel_phi = vg_azimuthal - phi_dot * r0;

    const double vrel =
	std::sqrt(std::pow(minus_r_dotel_r, 2) + std::pow(vrel_phi, 2));

    minus_l_rel = r * vg_azimuthal - l0;

	tstop = calc_tstop_only(radius, rho, vrel, temperature);
}

// confirm that tstop < dt/10, else make user change integrator
// tstop > dt/10 can cause numerical instabilities for the explicit integrator
void check_tstop(t_data &data)
{
    double dt = sim::last_dt;

	#pragma omp parallel for
    for (unsigned int i = 0; i < local_number_of_particles; ++i) {
	const double radius = particles[i].radius;
	const double r = particles[i].get_distance_to_star();
	const double phi = particles[i].get_angle();
	const double r_dot = particles[i].get_r_dot();
	const double phi_dot = particles[i].get_phi_dot();

	double rho, temperature, vg_radial, vg_azimuthal;
	interpolate_quantities(data, r, phi, rho, temperature, vg_radial, vg_azimuthal);
	
	// calculate relative velocities
	double minus_r_dotel_r = vg_radial - r_dot;
	const double vrel_phi = vg_azimuthal - phi_dot * r;

	const double vrel = std::sqrt(std::pow(minus_r_dotel_r, 2) + std::pow(vrel_phi, 2));


	const double tstop = calc_tstop_only(radius, rho, vrel, temperature);
	const double stokes = tstop * calculate_omega_kepler(r);
	particles[i].stokes = stokes;

	if (tstop < 10.0 * dt && parameters::particle_gas_drag_enabled &&
		(parameters::particle_integrator == parameters::integrator_explicit ||
		 parameters::particle_integrator == parameters::integrator_adaptive)) {
	    logging::print_master(
		LOG_ERROR
		"Particle stopping time too small for explicit integrator! Use the semiimplicit or implicit integrator instead!\n");
	    PersonalExit(1);
	}

	if (tstop < 0.005 * dt && parameters::particle_gas_drag_enabled &&
		parameters::particle_integrator == parameters::integrator_semiimplicit) {
	    logging::print_master(
		LOG_ERROR
		"Particle stopping time too small for semiimplicit integrator! Use the implicit integrator instead!\n");
	    PersonalExit(1);
	}
    }
}

// apply disk feedback on primary onto the particles
void update_velocities_from_indirect_term(const double dt)
{
    // Naming of r and phi weird!!!!!!!
	#pragma omp parallel for
    for (unsigned int i = 0; i < local_number_of_particles; ++i) {

	double indirect_q1_dot;
	double indirect_q2_dot;
	if (parameters::CartesianParticles) {
	    indirect_q1_dot = refframe::IndirectTerm.x * dt + particles[i].r_dot;
	    indirect_q2_dot = refframe::IndirectTerm.y * dt + particles[i].phi_dot;
	} else {
	    double r = particles[i].r;
	    double phi = particles[i].phi;

	    const double r_dot = particles[i].r_dot;
	    const double phi_dot = particles[i].phi_dot;

	    indirect_q1_dot = r_dot + dt * (refframe::IndirectTerm.x * std::cos(phi) +
					    refframe::IndirectTerm.y * std::sin(phi));
	    indirect_q2_dot = phi_dot + dt *
					    (-refframe::IndirectTerm.x * std::sin(phi) +
					     refframe::IndirectTerm.y * std::cos(phi)) /
					    r;
	}

	particles[i].r_dot = indirect_q1_dot;
	particles[i].phi_dot = indirect_q2_dot;
    }
}

void update_velocities_from_gas_drag_cart(t_data &data, double dt)
{

	#pragma omp parallel for
    for (unsigned int i = 0; i < local_number_of_particles; ++i) {
	const double r = particles[i].get_distance_to_star();
	const double phi = particles[i].get_angle();

	// check if particle has left disc
	if ((r < RMIN) || (r > RMAX)) {
	    continue;
	}

	unsigned int n_radial_a_minus = 0, n_radial_a_plus = 1,
		 n_radial_b_minus = 0, n_radial_b_plus = 1;
    unsigned int n_azimuthal_a_minus = 0, n_azimuthal_a_plus = 0,
		 n_azimuthal_b_minus = 0, n_azimuthal_b_plus = 0;
    find_nearest(n_radial_a_minus, n_radial_a_plus, n_azimuthal_a_minus,
		 n_azimuthal_a_plus, n_radial_b_minus, n_radial_b_plus,
		 n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);

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
    
	const double vg_azimuthal_temp = interpolate_bilinear(
	data[t_data::V_AZIMUTHAL], false, true, n_radial_b_minus,
	n_radial_b_plus, n_azimuthal_a_minus, n_azimuthal_a_plus, r, phi);
    
	const double vg_azimuthal = corret_v_gas_azimuthal_omega_frame(vg_azimuthal_temp, r);

	// calculate gas velocities in cartesian coordinates
	const double vg_x = std::cos(phi) * vg_radial - std::sin(phi) * vg_azimuthal;
	const double vg_y = std::sin(phi) * vg_radial + std::cos(phi) * vg_azimuthal;

	// particles store cartesian data
	double &vx = particles[i].r_dot;
	double &vy = particles[i].phi_dot;

	// calculate relative velocities
	const double vrel_x = vx - vg_x;
	const double vrel_y = vy - vg_y;
	const double vrel = std::sqrt(std::pow(vrel_x, 2) + std::pow(vrel_y, 2));


	const double size = particles[i].radius;
	const double tstop = calc_tstop_only(size, rho, vrel, temperature);

	const double a_drag_x = - vrel_x/tstop;
	const double a_drag_y = - vrel_y/tstop;

	// update velocities
	vx += dt * a_drag_x;
	vy += dt * a_drag_y;

	// update stokes number
	const double omega_kepler = calculate_omega_kepler(r);
	particles[i].stokes = tstop*omega_kepler;

	if (parameters::particle_disk_gravity_enabled) {
	    update_velocity_from_disk_gravity_cart(
		n_radial_a_minus, n_radial_a_plus, n_azimuthal_b_minus,
		n_azimuthal_b_plus, r, phi, i, dt);
	}
    } // loop end

    /*
    // For testing purpose only, very slow
    if (parameters::particle_disk_gravity_enabled) {
	    update_velocity_from_disk_gravity_cart_old(data, dt);
    }
    */
}

void update_velocities_from_gas_drag(t_data &data, double dt)
{

	#pragma omp parallel for
    for (unsigned int i = 0; i < local_number_of_particles; ++i) {
	double r = particles[i].get_distance_to_star();
	double phi = particles[i].get_angle();

	// check if particle has left disc
	if ((r < RMIN) || (r > RMAX)) {
	    continue;
	}

	// find nearest cells
	unsigned int n_radial_a_minus = 0, n_radial_a_plus = 1,
		 n_radial_b_minus = 0, n_radial_b_plus = 1;
    unsigned int n_azimuthal_a_minus = 0, n_azimuthal_a_plus = 0,
		 n_azimuthal_b_minus = 0, n_azimuthal_b_plus = 0;
    find_nearest(n_radial_a_minus, n_radial_a_plus, n_azimuthal_a_minus,
		 n_azimuthal_a_plus, n_radial_b_minus, n_radial_b_plus,
		 n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);

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
    
	const double vg_azimuthal_temp = interpolate_bilinear(
	data[t_data::V_AZIMUTHAL], false, true, n_radial_b_minus,
	n_radial_b_plus, n_azimuthal_a_minus, n_azimuthal_a_plus, r, phi);
    
	const double vg_azimuthal = corret_v_gas_azimuthal_omega_frame(vg_azimuthal_temp, r);

	// calculate relative velocities
	const double vrel_r = particles[i].r_dot - vg_radial;
	const double vrel_phi = r * particles[i].phi_dot - vg_azimuthal;
	const double vrel =
	    std::sqrt(std::pow(vrel_r, 2) + std::pow(vrel_phi, 2));



	const double size = particles[i].radius;
	const double tstop = calc_tstop_only(size, rho, vrel, temperature);

	// const double fdrag_temp = -0.5 * Cd * M_PI *
	// 			  std::pow(particles[i].radius, 2) * rho *
	// 			  vrel / particles[i].mass;

	// const double fdrag_r = fdrag_temp * vrel_r;
	// const double fdrag_phi = fdrag_temp * vrel_phi / r;

	const double a_drag_r = - vrel_r / tstop;
	const double a_drag_phi = - vrel_phi / r / tstop;

	// update stokes number
	const double omega_kepler = calculate_omega_kepler(r);
	particles[i].stokes = tstop*omega_kepler;

	// update velocities
	particles[i].r_dot += dt*a_drag_r;
	particles[i].phi_dot += dt*a_drag_phi;

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

void integrate(t_data &data, const double current_time, const double dt)
{
	// const double q1 = particles[0].r;
	// const double q2 = particles[0].phi;
	// const double stokes = particles[0].stokes;
	// printf("t = %.3e, q1 = %.3e, q2 = %.3e, |q| = %.3e, St = %.3e\n", current_time, q1, q2, std::sqrt(q1*q1 + q2*q2), stokes);
	if (parameters::particle_gas_drag_enabled){
		compute_rho(data, current_time);
	}

	switch (parameters::particle_integrator) {
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
    case parameters::integrator_exponential_midpoint: {
	integrate_exponential_midpoint(data, dt);
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
	if (parameters::particle_dust_diffusion) {
		if (!parameters::particle_gas_drag_enabled) {
			check_tstop(data);
		}
		// TODO: should be before corrector step for implicit method: see Picogna+2018 App. B.2
		dust_diffusion::diffuse_dust(data, particles, dt, local_number_of_particles);
    }
	move();

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

    vx += dt * (std::cos(phi) * sg_radial - std::sin(phi) * sg_azimuthal);
    vy += dt * (std::sin(phi) * sg_radial + std::cos(phi) * sg_azimuthal);
}

/**
 * @brief compute_smoothing_isothermal
 * @param r
 * @return smoothing length for calculating gravitational force due to a
 * vertical gas column
 */
static double compute_smoothing_isothermal(double r)
{
    double smooth;
    const double scale_height =
	parameters::ASPECTRATIO_REF * std::pow(r, 1.0 + parameters::FLARINGINDEX); // = H
    smooth = parameters::thickness_smoothing * scale_height;
    return smooth;
}
/**
 * @brief update_velocity_from_disk_gravity_cart_old
 * @abstract horrible inefficient way to calculate gravitational force from gas
 * on dust particles, only use for testing
 * @param data
 * @param dt
 */
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
    MPI_Allgatherv(&particles[0], local_number_of_particles, mpi_particle,
		   all_particles, number_of_particles, particle_offsets,
		   mpi_particle, MPI_COMM_WORLD);

    std::memset(force_x, 0, sizeof(*force_x) * global_number_of_particles);
    std::memset(force_y, 0, sizeof(*force_y) * global_number_of_particles);

    double dphi = 2.0 * M_PI / (double)data[t_data::SIGMA].get_size_azimuthal();
    for (unsigned int n_radial = Zero_or_active; n_radial < Max_or_active;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::SIGMA].get_max_azimuthal();
	     ++n_azimuthal) {
	    double cell_angle = n_azimuthal * dphi;
	    double cell_x = Rmed[n_radial] * std::cos(cell_angle);
	    double cell_y = Rmed[n_radial] * std::sin(cell_angle);
	    double cell_mass =
		Surf[n_radial] * data[t_data::SIGMA](n_radial, n_azimuthal);
	    for (unsigned int i = 0; i < global_number_of_particles; ++i) {
		const double &x = all_particles[i].r;
		const double &y = all_particles[i].phi;
		double smoothing = compute_smoothing_isothermal(
		    all_particles[i].get_distance_to_star());
		double d_x = cell_x - x;
		double d_y = cell_y - y;
		double invdist3 = std::pow(std::pow(d_x, 2) + std::pow(d_y, 2) +
					       std::pow(smoothing, 2),
					   -3.0 / 2.0);

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
	#pragma omp parallel for
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

void integrate_exponential_midpoint(t_data &data, const double dt)
{
    // Semi implicit integrator in cylindrical coordinates (see Zhu et al. 2014,
    // eqs. A4-A12)
	#pragma omp parallel for
    for (unsigned int i = 0; i < local_number_of_particles; ++i) {

	// initialize
	const double r0 = particles[i].r;
	const double phi0 = particles[i].phi;

	const double r_dot0 = particles[i].r_dot;
	const double phi_dot0 = particles[i].phi_dot;

	const double l0 = r0 * r0 * phi_dot0;

	const double hfdt = 0.5 * dt;
	double tstop = dt;
	double vrel_r = 0.0;

	// Half-drift ////////////////////////////////////////////
	const double r_dot1 = r_dot0;
	const double l1 = l0;
	const double phi_dot1 = phi_dot0; // follows from  l1 = l0

	const double r1 = r0 + r_dot0 * hfdt;
	double phi1 =
	    phi0 + 0.5 * (l0 / std::pow(r0, 2) + l0 / std::pow(r1, 2)) * hfdt;
	check_angle(phi1);
	// END Half-drift
	// /////////////////////////////////////////////////////////

	// Kick ///////////////////////////////////////////////////
	const double r2 = r1;
	const double phi2 = phi1;
	double minus_l_rel = 0.0;

	if (parameters::particle_gas_drag_enabled) {
	    calculate_tstop2(r1, phi1, r_dot1, phi_dot1, data,
			     particles[i].radius, vrel_r, minus_l_rel, tstop,
			     r0, l0);
	} else {
	    tstop = 1e100;
	}

	const double rsmooth = calculate_dust_smoothing(r1, phi1, particles[i].stokes, data);

	double grav_r_ddot;
	double minus_grav_l_dot;
	if (parameters::ParticlesInCartesian) {
	    calculate_derivitives_from_star_and_planets_in_cart(
		grav_r_ddot, minus_grav_l_dot, r1, phi1, rsmooth, data);
	} else {

	    calculate_derivitives_from_star_and_planets(
		grav_r_ddot, minus_grav_l_dot, r1, phi1, rsmooth, data);
	}

	// exponential propagator  eq.33 (Mignone et al. 2019)
	const double exp_tstop = std::exp(-dt / tstop);
	const double h1 = tstop * (-expm1(-dt / tstop));

	// updating angular momentum
	double l2 = exp_tstop * l1 + h1 * minus_grav_l_dot;
	if (parameters::particle_gas_drag_enabled) {
	    const double l_gas = minus_l_rel + l0;
	    l2 += h1 * l_gas / tstop;
	}

	// Updating radial velocity
	double r_dot2 = exp_tstop * r_dot1;
	r_dot2 += h1 * 0.5 * (l1 * l1 + l2 * l2) / (r1 * r1 * r1);
	r_dot2 += h1 * grav_r_ddot;

	if (parameters::particle_gas_drag_enabled) {
	    const double v_r_g = vrel_r + r_dot1;
	    r_dot2 += h1 * v_r_g / tstop;
	}
	// END kick ///////////////////////////////////////

	// Half-drift
	const double r3 = r1 + r_dot2 * hfdt;
	const double phi3 =
	    phi2 + 0.5 * (l2 / std::pow(r2, 2) + l2 / std::pow(r3, 2)) * hfdt;

	// Update
	particles[i].r_dot = r_dot2;
	particles[i].r = r3;
	particles[i].phi = phi3;
	check_angle(particles[i].phi);
	particles[i].phi_dot = l2 / std::pow(particles[i].r, 2);
	particles[i].stokes = tstop*calculate_omega_kepler(r3);
    }
	// TODO: check for multiple calls of move
}

void integrate_semiimplicit(t_data &data, const double dt)
{
    // Semi implicit integrator in cylindrical coordinates (see Zhu et al. 2014,
    // eqs. A4-A12)
	#pragma omp parallel for
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
	double vrel_r = 0.0;

	// Half-drift
	const double r_dot1 = r_dot0;
	const double l1 = l0;
	const double phi_dot1 = phi_dot0; // follows from  l1 = l0

	const double r1 = r0 + r_dot0 * hfdt;
	double phi1 =
	    phi0 + 0.5 * (l0 / std::pow(r0, 2) + l0 / std::pow(r1, 2)) * hfdt;
	check_angle(phi1);

	// Kick
	const double r2 = r1;
	const double phi2 = phi1;
	double minus_l_rel = 0.0;

	if (parameters::particle_gas_drag_enabled) {
	    calculate_tstop2(r1, phi1, r_dot1, phi_dot1, data,
			     particles[i].radius, vrel_r, minus_l_rel, tstop,
			     r0, l0);
	    dt1 = dt / (1.0 + hfdt / tstop);
	}

	const double rsmooth = calculate_dust_smoothing(r1, phi1, particles[i].stokes, data);

	double grav_r_ddot;
	double minus_grav_l_dot;
	if (parameters::ParticlesInCartesian) {
	    calculate_derivitives_from_star_and_planets_in_cart(
		grav_r_ddot, minus_grav_l_dot, r1, phi1, rsmooth, data);
	} else {

	    calculate_derivitives_from_star_and_planets(
		grav_r_ddot, minus_grav_l_dot, r1, phi1, rsmooth, data);
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
	const double phi3 =
	    phi2 + 0.5 * (l2 / std::pow(r2, 2) + l2 / std::pow(r3, 2)) * hfdt;

	particles[i].r_dot = r_dot2;
	particles[i].r = r3;
	particles[i].phi = phi3;
	check_angle(particles[i].phi);
	particles[i].phi_dot = l2 / std::pow(particles[i].r, 2);
	particles[i].stokes = tstop*calculate_omega_kepler(r3);
    }
}

void integrate_implicit(t_data &data, const double dt)
{
    // Fully implicit integrator in cylindrical coordinates (see Zhu et al.
    // 2014, eqs. A15-A18)
    double r1, phi1, l1, r_dot1, tstop1, tstop0, dt0, dt1;
    double r_ddot0, minus_l_dot0, r_ddot1, minus_l_dot1, hfdt, minus_r_dot_rel0,
	minus_l_rel0, minus_r_dot_rel1;

    double minus_l_rel1 = 0.0;

    // initialize with failsave values to suppress compiler warning
    // "-Wmaybe-uninitialized"
    minus_l_rel0 = 0.0;
    minus_r_dot_rel0 = 0.0;
    tstop0 = 1e+300;

	#pragma omp parallel for
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
	phi1 = phi0 + 0.5 * (l0 / std::pow(r0, 2) + l0 / std::pow(r1, 2)) * dt;
	check_angle(phi1);

	// Kick
	// Current position
	if (parameters::particle_gas_drag_enabled) {
	    calculate_tstop(r0, phi0, r_dot0, phi_dot0, data,
			    particles[i].radius, minus_r_dot_rel0, minus_l_rel0,
			    tstop0);
	    dt0 = 1.0 + dt / tstop0;
	}

	const double rsmooth = calculate_dust_smoothing(r1, phi1, particles[i].stokes, data);

	if (parameters::ParticlesInCartesian) {
	    calculate_derivitives_from_star_and_planets_in_cart(
		r_ddot0, minus_l_dot0, r0, phi0, rsmooth, data);
	} else {
	    calculate_derivitives_from_star_and_planets(r_ddot0, minus_l_dot0,
							r0, phi0, rsmooth, data);
	}
	// Predicted position
	if (parameters::particle_gas_drag_enabled) {
	    calculate_tstop2(r1, phi1, r_dot0, phi_dot0, data,
			     particles[i].radius, minus_r_dot_rel1,
			     minus_l_rel1, tstop1, r0, l0);
	    dt1 = hfdt / (1.0 + hfdt / tstop0 + hfdt / tstop1 +
			  hfdt * dt / (tstop0 * tstop1));
	}

	if (parameters::ParticlesInCartesian) {
	    calculate_derivitives_from_star_and_planets_in_cart(
		r_ddot1, minus_l_dot1, r1, phi1, rsmooth, data);
	} else {
	    calculate_derivitives_from_star_and_planets(r_ddot1, minus_l_dot1,
							r1, phi1, rsmooth, data);
	}

	l1 = l0;
	l1 += (minus_l_dot0 + minus_l_dot1 * dt0) * dt1;
	if (parameters::particle_gas_drag_enabled) {
	    l1 += (minus_l_rel0 / tstop0 + minus_l_rel1 / tstop1 * dt0) * dt1;
	}

	r_dot1 = r_dot0;
	r_dot1 += (r_ddot0 + r_ddot1 * dt0) * dt1;
	r_dot1 += (std::pow(l0, 2) / std::pow(r0, 3) +
		   std::pow(l1, 2) / std::pow(r1, 3) * dt0) *
		  dt1;
	if (parameters::particle_gas_drag_enabled) {
	    r_dot1 +=
		(minus_r_dot_rel0 / tstop0 + minus_r_dot_rel1 / tstop1 * dt0) *
		dt1;
	}

	// Half-drift
	const double r2 = r0 + (r_dot0 + r_dot1) * hfdt;
	const double phi2 =
	    phi0 + (l0 / std::pow(r0, 2) + l1 / std::pow(r2, 2)) * hfdt;
	particles[i].r_dot = r_dot1;
	particles[i].r = r2;
	particles[i].phi = phi2;
	check_angle(particles[i].phi);
	particles[i].phi_dot = l1 / std::pow(particles[i].r, 2);
	particles[i].stokes = tstop1*calculate_omega_kepler(r2);
    }
}

void integrate_explicit_adaptive(t_data &data, const double dt)
{

    if (parameters::CartesianParticles) {
		die("Somethings wrong with integrating particles with the explicit method and cartesian calculation of forces, no clue what. Find it out!\n");
	// disk gravity on particles is inside gas_drag function
	if (parameters::particle_gas_drag_enabled)
	    update_velocities_from_gas_drag_cart(data, dt);
    } else {
	// disk gravity on particles is inside gas_drag function
	if (parameters::particle_gas_drag_enabled)
	    update_velocities_from_gas_drag(data, dt);
    }

    // as particles move independent of each other, we can integrate one after
    // one
	#pragma omp parallel for
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

	    const double rsmooth = calculate_dust_smoothing(
		particles[i].get_distance_to_star(), particles[i].get_angle(),
		particles[i].stokes, data);

	    // Cash–Karp method
	    // (http://en.wikipedia.org/wiki/Cash%E2%80%93Karp_method) cartesian
	    // coordinates are written inside the polar coordinates
	    if (parameters::CartesianParticles) {
		calculate_accelerations_from_star_and_planets_cart(
		    temp_ar, temp_aphi, temp_r, temp_phi, rsmooth, data);
	    } else {
		calculate_accelerations_from_star_and_planets(
		    temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi,
		    temp_phi_dot, rsmooth, data);
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

	    if (parameters::CartesianParticles) {
		calculate_accelerations_from_star_and_planets_cart(
		    temp_ar, temp_aphi, temp_r, temp_phi, rsmooth, data);
	    } else {
		calculate_accelerations_from_star_and_planets(
		    temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi,
		    temp_phi_dot, rsmooth, data);
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

	    if (parameters::CartesianParticles) {
		calculate_accelerations_from_star_and_planets_cart(
		    temp_ar, temp_aphi, temp_r, temp_phi, rsmooth, data);
	    } else {
		calculate_accelerations_from_star_and_planets(
		    temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi,
		    temp_phi_dot, rsmooth, data);
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

	    if (parameters::CartesianParticles) {
		calculate_accelerations_from_star_and_planets_cart(
		    temp_ar, temp_aphi, temp_r, temp_phi, rsmooth, data);
	    } else {
		calculate_accelerations_from_star_and_planets(
		    temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi,
		    temp_phi_dot, rsmooth, data);
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

	    if (parameters::CartesianParticles) {
		calculate_accelerations_from_star_and_planets_cart(
		    temp_ar, temp_aphi, temp_r, temp_phi, rsmooth, data);
	    } else {
		calculate_accelerations_from_star_and_planets(
		    temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi,
		    temp_phi_dot, rsmooth, data);
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

	    if (parameters::CartesianParticles) {
		calculate_accelerations_from_star_and_planets_cart(
		    temp_ar, temp_aphi, temp_r, temp_phi, rsmooth, data);
	    } else {
		calculate_accelerations_from_star_and_planets(
		    temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi,
		    temp_phi_dot, rsmooth, data);
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

	    if (!parameters::CartesianParticles) {
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

	    err = std::sqrt(err / (double)(4));
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

    if (parameters::CartesianParticles) {
		die("Somethings wrong with integrating particles with the explicit method and cartesian calculation of forces, no clue what. Find it out!\n");
	// disk gravity on particles is inside gas_drag function
	if (parameters::particle_gas_drag_enabled)
	    update_velocities_from_gas_drag_cart(data, dt);
    } else {
	// disk gravity on particles is inside gas_drag function
	if (parameters::particle_gas_drag_enabled)
	    update_velocities_from_gas_drag(data, dt);
    } // else {
    // disk gravity on particles is inside gas_drag function

    // as particles move independent of each other, we can integrate one after
    // one
	#pragma omp parallel for
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

	const double rsmooth = calculate_dust_smoothing(
	particles[i].get_distance_to_star(), particles[i].get_angle(),
	particles[i].stokes, data);

	// Cash–Karp method
	// (http://en.wikipedia.org/wiki/Cash%E2%80%93Karp_method) cartesian
	// coordinates are written inside the polar coordinates
	if (parameters::CartesianParticles) {
	    calculate_accelerations_from_star_and_planets_cart(
		temp_ar, temp_aphi, temp_r, temp_phi, rsmooth, data);
	} else {
	    calculate_accelerations_from_star_and_planets(
		temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi, temp_phi_dot,
		rsmooth, data);
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

	if (parameters::CartesianParticles) {
	    calculate_accelerations_from_star_and_planets_cart(
		temp_ar, temp_aphi, temp_r, temp_phi, rsmooth, data);
	} else {
	    calculate_accelerations_from_star_and_planets(
		temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi, temp_phi_dot,
		rsmooth, data);
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

	if (parameters::CartesianParticles) {
	    calculate_accelerations_from_star_and_planets_cart(
		temp_ar, temp_aphi, temp_r, temp_phi, rsmooth, data);
	} else {
	    calculate_accelerations_from_star_and_planets(
		temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi, temp_phi_dot,
		rsmooth, data);
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

	if (parameters::CartesianParticles) {
	    calculate_accelerations_from_star_and_planets_cart(
		temp_ar, temp_aphi, temp_r, temp_phi, rsmooth, data);
	} else {
	    calculate_accelerations_from_star_and_planets(
		temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi, temp_phi_dot,
		rsmooth, data);
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

	if (parameters::CartesianParticles) {
	    calculate_accelerations_from_star_and_planets_cart(
		temp_ar, temp_aphi, temp_r, temp_phi, rsmooth, data);
	} else {
	    calculate_accelerations_from_star_and_planets(
		temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi, temp_phi_dot,
		rsmooth, data);
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

	if (parameters::CartesianParticles) {
	    calculate_accelerations_from_star_and_planets_cart(
		temp_ar, temp_aphi, temp_r, temp_phi, rsmooth, data);
	} else {
	    calculate_accelerations_from_star_and_planets(
		temp_ar, temp_aphi, temp_r, temp_r_dot, temp_phi, temp_phi_dot,
		rsmooth, data);
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

	if (!parameters::CartesianParticles) {
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
    double local_r_min_squared = std::pow(local_r_min, 2);
    double local_r_max_squared = std::pow(local_r_max, 2);

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
	    particles.resize(particles_size);
	}

	MPI_Recv(&particles[0] + local_number_of_particles, number, mpi_particle,
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
	MPI_Send(&particles[0], 1, inward_type, CPU_Prev, 0, MPI_COMM_WORLD);
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
	    particles.resize(particles_size);
	}

	MPI_Recv(&particles[0] + local_number_of_particles, number, mpi_particle,
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
	MPI_Send(&particles[0], 1, outward_type, CPU_Next, 0, MPI_COMM_WORLD);
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

void write()
{
    MPI_File fh;
    MPI_Status status;

    std::string filename = output::snapshot_dir + "/particles.dat";

    // try to open file
    mpi_error_check_file_write(
	MPI_File_open(MPI_COMM_WORLD, filename.c_str(),
		      MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh),
	filename.c_str());

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
    MPI_File_write(fh, &particles[0], local_number_of_particles, mpi_particle,
		   &status);

    // close file
    MPI_File_close(&fh);
}

void rotate(const double angle)
{
	#pragma omp parallel for
    for (unsigned int i = 0; i < local_number_of_particles; ++i) {
	// rotate positions
	if (parameters::CartesianParticles) {
	    const double x_old = particles[i].r;
	    const double y_old = particles[i].phi;

	    const double vx_old = particles[i].r_dot;
	    const double vy_old = particles[i].phi_dot;

	    double &x = particles[i].r;
	    double &y = particles[i].phi;

	    double &vx = particles[i].r_dot;
	    double &vy = particles[i].phi_dot;

	    // rotate positions
	    x = x_old * std::cos(angle) + y_old * std::sin(angle);
	    y = -x_old * std::sin(angle) + y_old * std::cos(angle);
	    // rotate velocities
	    vx = vx_old * std::cos(angle) + vy_old * std::sin(angle);
	    vy = -vx_old * std::sin(angle) + vy_old * std::cos(angle);
	} else {
	    particles[i].phi -= angle;
	    check_angle(particles[i].phi);
	}
    }
}

} // namespace particles
