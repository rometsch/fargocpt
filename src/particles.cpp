/**
	\file particles.cpp
	\author Tobias Mueller <Tobias_Mueller@twam.info>
*/

#include <stdlib.h>
#include <random>
#include <math.h>
#include <mpi.h>
#include "LowTasks.h"
#include "util.h"
#include "particles.h"
#include "logging.h"
#include "parameters.h"
#include "constants.h"
#include "global.h"
#include "SourceEuler.h"
#include "Force.h"
#include "Pframeforce.h"

extern Pair IndirectTerm;

namespace particles {

/// local particle storage
t_particle* particles;

/// current size of local particle storage
unsigned int particles_size;

/// current number of particles on this node
unsigned int local_number_of_particles;

/// global number of particles;
unsigned int global_number_of_particles;

double local_r_min;
double local_r_max;

static MPI_Datatype mpi_particle;

static double interpolate_bilinear(t_polargrid &quantity, bool radial_a_grid, bool azimuthal_a_grid, unsigned int n_radial_minus, unsigned int n_radial_plus, unsigned int n_azimuthal_minus, unsigned int n_azimuthal_plus, double r, double phi);

// inverse of the Cumulative distribution function for a slope proportional to r^n
static double f(double x, double n) {
	double rmin = parameters::particle_minimum_radius;
	double rmax = parameters::particle_maximum_radius;

	if (n == -1) {
		return exp(log(rmax-rmin)*x)*rmin;
	} else {
		return exp(log(x*(-pow(rmin,n+1)+pow(rmax,n+1))+pow(rmin,n+1))/(n+1));
	}
}

void init() {
	// calculate number of initial local particles
	global_number_of_particles = parameters::number_of_particles;
	local_number_of_particles = global_number_of_particles/CPU_Number;

	if ((unsigned int)CPU_Rank < global_number_of_particles-CPU_Number*local_number_of_particles) {
		local_number_of_particles++;
	}

	local_r_min = Ra[CPU_Rank == 0 ? GHOSTCELLS_A-1 : CPUOVERLAP];
	local_r_max = Ra[CPU_Rank == CPU_Highest ? NRadial-GHOSTCELLS_A+1 : NRadial-CPUOVERLAP];

	// create storage
	particles_size = local_number_of_particles;
	particles = (t_particle*)malloc(sizeof(t_particle)*particles_size);

	// init particles with data
	unsigned int id_offset = global_number_of_particles/CPU_Number * CPU_Rank;
	if ((unsigned int)CPU_Rank < global_number_of_particles-CPU_Number*local_number_of_particles) {
		id_offset+=CPU_Rank;
	} else {
		id_offset+=global_number_of_particles-CPU_Number*local_number_of_particles;
	}

    logging::print(LOG_DEBUG "random generator seed: %u\n", id_offset);

    // random generator and distributions
    std::mt19937 generator(id_offset*parameters::random_seed);
    std::uniform_real_distribution<double> dis_one(0.0, 1.0);
    std::uniform_real_distribution<double> dis_twoPi(0.0, 2.0*PI);


	for (unsigned int i = 0; i < local_number_of_particles; ++i) {
        double semi_major_axis = f(dis_one(generator), parameters::particle_slope);
        double phi = dis_twoPi(generator);

		double eccentricity = 0.0;

		particles[i].radius = parameters::particle_radius;
		double volume = 4.0/3.0*PI*pow3(particles[i].radius);

		particles[i].mass = volume * parameters::particle_density;

		double r = semi_major_axis*(1.0+eccentricity);
		double v = sqrt(constants::G*(M+particles[i].mass)/semi_major_axis)*sqrt((1.0-eccentricity)/(1.0+eccentricity));

		particles[i].x = r*cos(phi);
		particles[i].y = r*sin(phi);

        particles[i].vx = -v*sin(phi);
        particles[i].vy = v*cos(phi);

		particles[i].id = id_offset+i;
	}

	// create MPI datatype
	int mpi_particle_count = 7;
	int mpi_particle_lengths[7] = {1, 1, 1, 1, 1, 1, 1};
	MPI_Aint mpi_particle_offsets[7];
	MPI_Datatype mpi_particle_types[7] = {MPI_UNSIGNED, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

	// calculate offsets
	MPI_Aint base;
    MPI_Get_address(particles, &base);
    MPI_Get_address(&particles[0].id, mpi_particle_offsets);
    MPI_Get_address(&particles[0].x, mpi_particle_offsets+1);
    MPI_Get_address(&particles[0].y, mpi_particle_offsets+2);
    MPI_Get_address(&particles[0].vx, mpi_particle_offsets+3);
    MPI_Get_address(&particles[0].vy, mpi_particle_offsets+4);
    MPI_Get_address(&particles[0].mass, mpi_particle_offsets+5);
    MPI_Get_address(&particles[0].radius, mpi_particle_offsets+6);

	for (int i = 0; i < 7; ++i) {
		mpi_particle_offsets[i] -= base;
	}

	MPI_Type_create_struct(mpi_particle_count, mpi_particle_lengths, mpi_particle_offsets, mpi_particle_types, &mpi_particle);
	MPI_Type_commit(&mpi_particle);

	// particles might be created on the wrong node, so move them to the correct one :)
	move();
}

void calculate_accelerations_from_star_and_planets(double &ax, double &ay, double x, double y, t_data& data) {
	ax = 0;
	ay = 0;

	// host star
	double r2 = pow2(x) + pow2(y);
	double factor = constants::G*M*pow(r2,-3.0/2.0);
	ax += factor*(-x);
	ay += factor*(-y);

	// planets
	for (unsigned int k = 0; k < data.get_planetary_system().get_number_of_planets(); ++k) {
		t_planet &planet = data.get_planetary_system().get_planet(k);
		double r2 = pow2(planet.get_x()-x) + pow2(planet.get_y()-y);
		double factor = constants::G*planet.get_mass()*pow(r2,-3.0/2.0);
		ax += factor*(planet.get_x()-x);
		ay += factor*(planet.get_y()-y);
	}
}

static double interpolate_bilinear(t_polargrid &quantity, bool radial_a_grid, bool azimuthal_a_grid, unsigned int n_radial_minus, unsigned int n_radial_plus, unsigned int n_azimuthal_minus, unsigned int n_azimuthal_plus, double r, double phi) {
	double dphi = 2.0*PI/(double)quantity.get_size_azimuthal();

	// values at corners
	double Qmm = quantity(n_radial_minus, n_azimuthal_minus);
	double Qpm = quantity(n_radial_plus, n_azimuthal_minus);
	double Qmp = quantity(n_radial_minus, n_azimuthal_plus);
	double Qpp = quantity(n_radial_plus, n_azimuthal_plus);

	double rm = radial_a_grid ? Ra[n_radial_minus] : Rb[n_radial_minus];
	double rp = radial_a_grid ? Ra[n_radial_plus] : Rb[n_radial_plus];

	double phim = azimuthal_a_grid ? n_azimuthal_minus*dphi : (n_azimuthal_minus+0.5)*dphi;
	double phip = azimuthal_a_grid ? n_azimuthal_plus*dphi : (n_azimuthal_plus+0.5)*dphi;

	if ((phim > 2*PI) || (phip > 2*PI)) {
		phim -= PI;
		phip -= PI;
		phi -= PI;
	}

	double Qm = ((rm*phip-rm*phi)*Qmm+(rm*phi-rm*phim)*Qmp)/(rm*phip-rm*phim);
	double Qp = ((rp*phip-rp*phi)*Qpm+(rp*phi-rp*phim)*Qpp)/(rp*phip-rp*phim);

	return ((rp-r)*Qm+(r-rm)*Qp)/(rp-rm);
}

void update_velocities_from_gas_drag(t_data &data, double dt) {
	double dphi = 2.0*PI/(double)data[t_data::DENSITY].get_size_azimuthal();

	compute_rho(data, true);

	for (unsigned int i = 0; i < local_number_of_particles; ++i) {
		double r = particles[i].get_distance_to_star();
		double phi = particles[i].get_angle();

		// check if particle has left disc
		if ((r<RMIN) || (r>RMAX)) {
			continue;
		}

		// find nearest n_radials
		unsigned int n_radial_a_minus = 0, n_radial_b_minus = 0, n_radial_a_plus = 1, n_radial_b_plus = 1;
		// we do not need to check if n_radial_a_plus < NRadial because it this would be the case, something in the move part went wrong :)
		while (Ra[n_radial_a_plus] < r)
			n_radial_a_plus++;
		n_radial_a_minus = n_radial_a_plus - 1;
		while (Rb[n_radial_b_plus] < r)
			n_radial_b_plus++;
		n_radial_b_minus = n_radial_b_plus - 1;

		// find nearest n_azimuthals
		unsigned int n_azimuthal_a_minus = 0, n_azimuthal_a_plus = 0;
		n_azimuthal_a_minus = trunc(phi/dphi);
		n_azimuthal_a_plus = n_azimuthal_a_minus+1;
		if (n_azimuthal_a_plus == data[t_data::DENSITY].get_size_azimuthal()) {
			n_azimuthal_a_plus = 0;
		}

		unsigned int n_azimuthal_b_minus = 0, n_azimuthal_b_plus = 0;
		n_azimuthal_b_minus = trunc((phi-0.5*dphi)/dphi);
		n_azimuthal_b_plus = n_azimuthal_b_minus+1;
		if (n_azimuthal_b_plus == data[t_data::DENSITY].get_size_azimuthal()) {
			n_azimuthal_b_plus = 0;
		}

		// calculate quantities need
		double rho = interpolate_bilinear(data[t_data::RHO], false, false, n_radial_b_minus, n_radial_b_plus, n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
		double vg_radial = interpolate_bilinear(data[t_data::V_RADIAL], true, false, n_radial_a_minus, n_radial_a_plus, n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);
		double vg_azimuthal = interpolate_bilinear(data[t_data::V_AZIMUTHAL], false, true, n_radial_b_minus, n_radial_b_plus, n_azimuthal_a_minus, n_azimuthal_a_plus, r, phi);
		double temperature = interpolate_bilinear(data[t_data::TEMPERATURE], false, false, n_radial_b_minus, n_radial_b_plus, n_azimuthal_b_minus, n_azimuthal_b_plus, r, phi);

		// calculate gas velocities in cartesian coordinates
		double vg_x = cos(phi) * vg_radial - sin(phi) * (vg_azimuthal + OmegaFrame*r );
		double vg_y = sin(phi) * vg_radial + cos(phi) * (vg_azimuthal + OmegaFrame*r );

		// calculate relative velocities
		double vrel_x = particles[i].vx - vg_x;
		double vrel_y = particles[i].vy - vg_y;
		double vrel = sqrt(pow2(vrel_x)+pow2(vrel_y));

		double m0 = parameters::MU * constants::m_u.get_code_value();
		double vthermal = sqrt(8.0*constants::k_B.get_code_value()*temperature/(PI*m0));
		// a0 = 1.5e-8 cm for molecular hydrogen
		double sigma = (PI*pow2(1.5e-8/units::length.get_cgs_factor()));
		double nu = 1.0/3.0*m0*vthermal/sigma;

		// calculate Reynolds number
		double reynolds = 2.0*rho*particles[i].radius*vrel/nu;

		// calculate coefficient
		double Cd;
		if (reynolds < 1.0) {
			Cd = 24.0/reynolds;
		} else if (reynolds < 800.0) {
			Cd = 24.0*pow(reynolds,-0.6);
		} else {
			Cd = 0.44;
		}

		double fdrag_temp = -0.5*Cd*PI*pow2(particles[i].radius)*rho*vrel;
		double fdrag_x = fdrag_temp * vrel_x;
		double fdrag_y = fdrag_temp * vrel_y;

		// update velocities
		particles[i].vx += dt * fdrag_x/particles[i].mass;
		particles[i].vy += dt * fdrag_y/particles[i].mass;
	}
}

void update_velocities_from_disk_gravity(t_data &data, double dt) {
	int *number_of_particles = (int*)malloc(sizeof(int)*CPU_Number);
	int *particle_offsets = (int*)malloc(sizeof(int)*CPU_Number);
	t_particle *all_particles = (t_particle*)malloc(sizeof(t_particle)*global_number_of_particles);
	double *force_x = (double*)malloc(sizeof(double)*global_number_of_particles);
	double *force_y = (double*)malloc(sizeof(double)*global_number_of_particles);

	// get particle amounts from all nodes
	MPI_Allgather(&local_number_of_particles, 1, MPI_INT, number_of_particles, 1, MPI_INT, MPI_COMM_WORLD);

	// calculate offsets
	for (int i = 0; i < CPU_Number; ++i) {
		particle_offsets[i] = 0;
		for (int j = 0; j < i; ++j) {
			particle_offsets[i] += number_of_particles[j];
		}
	}

	/*
	for (unsigned int i = 0; i < local_number_of_particles; ++i) {
		logging::print("Local Particle %u at %g,%g\n", i, particles[i].x, particles[i].y);
	} */

	// get all particles
	MPI_Allgatherv(particles, local_number_of_particles, mpi_particle, all_particles, number_of_particles, particle_offsets, mpi_particle, MPI_COMM_WORLD);

	/*
	for (unsigned int i = 0; i < global_number_of_particles; ++i) {
		logging::print("Global Particle %u at %g,%g\n", i, all_particles[i].x, all_particles[i].y);
	} */

	for (unsigned int i = 0; i < global_number_of_particles; ++i) {
		force_x[i] = 0;
		force_y[i] = 0;
	}

	double dphi = 2.0*PI/(double)data[t_data::DENSITY].get_size_azimuthal();
	for (unsigned int n_radial = Zero_or_active; n_radial < Max_or_active; ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
            double cell_angle = n_azimuthal*dphi;
			double cell_x =	Rmed[n_radial]*cos(cell_angle);
			double cell_y = Rmed[n_radial]*sin(cell_angle);
			double cell_mass = Surf[n_radial]*data[t_data::DENSITY](n_radial,n_azimuthal);
			for (unsigned int i = 0; i < global_number_of_particles; ++i) {
				double smoothing = compute_smoothing(all_particles[i].get_distance_to_star());
				double d_x = cell_x - all_particles[i].x;
				double d_y = cell_y - all_particles[i].y;
				double invdist3 = pow(pow2(d_x)+pow2(d_y)+pow2(smoothing),-3.0/2.0);

				force_x[i] += constants::G*cell_mass*d_x*invdist3;
				force_y[i] += constants::G*cell_mass*d_y*invdist3;
 			}
		}
	}

	// send force back
	MPI_Allreduce(MPI_IN_PLACE, force_x, global_number_of_particles, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, force_y, global_number_of_particles, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	Pair IndirectTerm = ComputeIndirectTerm();

	// update particles
	unsigned int offset = particle_offsets[CPU_Rank];
	for (unsigned int i = 0; i < local_number_of_particles; ++i) {
		particles[i].vx += dt * force_x[offset+i] + dt * IndirectTerm.x;
		particles[i].vy += dt * force_y[offset+i] + dt * IndirectTerm.y;
	}

	free(number_of_particles);
	free(particle_offsets);
	free(all_particles);
	free(force_x);
	free(force_y);
}

void integrate(t_data &data, double dt) {
	if ((parameters::calculate_disk) && (parameters::particle_gas_drag_enabled))
		update_velocities_from_gas_drag(data, dt);

	if (parameters::particle_disk_gravity_enabled) {
		update_velocities_from_disk_gravity(data, dt);
	}

	// as particles move independent of each other, we can integrate one after one
	for (unsigned int i = 0; i < local_number_of_particles; ++i) {
		double k1_x, k1_y, k1_vx, k1_vy;
		double k2_x, k2_y, k2_vx, k2_vy;
		double k3_x, k3_y, k3_vx, k3_vy;
		double k4_x, k4_y, k4_vx, k4_vy;
		double k5_x, k5_y, k5_vx, k5_vy;
		double k6_x, k6_y, k6_vx, k6_vy;

		double temp_ax, temp_ay, temp_x, temp_y;

		// Cash–Karp method (http://en.wikipedia.org/wiki/Cash%E2%80%93Karp_method)

		calculate_accelerations_from_star_and_planets(particles[i].ax, particles[i].ay, particles[i].x, particles[i].y, data);

		// calculate k1
		k1_vx = dt * particles[i].ax;
		k1_vy = dt * particles[i].ay;
		k1_x = dt * particles[i].vx;
		k1_y = dt * particles[i].vy;

		temp_x = particles[i].x + 0.2 * k1_x;
		temp_y = particles[i].y + 0.2 * k1_y;

		calculate_accelerations_from_star_and_planets(temp_ax, temp_ay, temp_x, temp_y, data);

		// calculate k2
		k2_vx = dt * temp_ax;
		k2_vy = dt * temp_ay;
		k2_x = dt * (particles[i].vx + 0.2 * k1_vx);
		k2_y = dt * (particles[i].vy + 0.2 * k1_vy);

		temp_x = particles[i].x + 0.075 * k1_x + 0.225 * k2_x;
		temp_y = particles[i].y + 0.075 * k1_y + 0.225 * k2_y;

		calculate_accelerations_from_star_and_planets(temp_ax, temp_ay, temp_x, temp_y, data);

		// calculate k3
		k3_vx = dt * temp_ax;
		k3_vy = dt * temp_ay;
		k3_x = dt * (particles[i].vx + 0.075 * k1_vx + 0.225 * k2_vx);
		k3_y = dt * (particles[i].vy + 0.075 * k1_vy + 0.225 * k2_vy);

		temp_x = particles[i].x + 0.3 * k1_x - 0.9 * k2_x + 1.2 * k3_x;
		temp_y = particles[i].y + 0.3 * k1_y - 0.9 * k2_y + 1.2 * k3_y;

		calculate_accelerations_from_star_and_planets(temp_ax, temp_ay, temp_x, temp_y, data);

		// calculate k4
		k4_vx = dt * temp_ax;
		k4_vy = dt * temp_ay;
		k4_x = dt * (particles[i].vx + 0.3 * k1_vx - 0.9 * k2_vx + 1.2 * k3_vx);
		k4_y = dt * (particles[i].vy + 0.3 * k1_vy - 0.9 * k2_vy + 1.2 * k3_vy);

		temp_x = particles[i].x - 11.0/54.0 * k1_x + 2.5 * k2_x - 70.0/27.0 * k3_x + 35.0/27.0 * k4_x;
		temp_y = particles[i].y - 11.0/54.0 * k1_y + 2.5 * k2_y - 70.0/27.0 * k3_y + 35.0/27.0 * k4_y;

		calculate_accelerations_from_star_and_planets(temp_ax, temp_ay, temp_x, temp_y, data);

		// calculate k5
		k5_vx = dt * temp_ax;
		k5_vy = dt * temp_ay;
		k5_x = dt * (particles[i].vx - 11.0/54.0 * k1_vx + 2.5 * k2_vx - 70.0/27.0 * k3_vx + 35.0/27.0 * k4_vx);
		k5_y = dt * (particles[i].vy - 11.0/54.0 * k1_vy + 2.5 * k2_vy - 70.0/27.0 * k3_vy + 35.0/27.0 * k4_vy);

		temp_x = particles[i].x + 1631.0/55296.0 * k1_x + 175.0/512.0 * k2_x + 575.0/13824.0 * k3_x + 44275.0/110592.0 * k4_x + 253.0/4096.0 * k5_x;
		temp_y = particles[i].y + 1631.0/55296.0 * k1_y + 175.0/512.0 * k2_y + 575.0/13824.0 * k3_y + 44275.0/110592.0 * k4_y + 253.0/4096.0 * k5_y;

		calculate_accelerations_from_star_and_planets(temp_ax, temp_ay, temp_x, temp_y, data);

		// calculate k6
		k6_vx = dt * temp_ax;
		k6_vy = dt * temp_ay;
		k6_x = dt * (particles[i].vx + 1631.0/55296.0 * k1_vx + 175.0/512.0 * k2_vx + 575.0/13824.0 * k3_vx + 44275.0/110592.0 * k4_vx + 253.0/4096.0 * k5_vx);
		k6_y = dt * (particles[i].vy + 1631.0/55296.0 * k1_vy + 175.0/512.0 * k2_vy + 575.0/13824.0 * k3_vy + 44275.0/110592.0 * k4_vy + 253.0/4096.0 * k5_vy);

		// update position & velocity
		particles[i].x += 37.0/378.0 * k1_x + 250.0/621.0 * k3_x + 125.0/594.0 * k4_x + 512.0/1771.0 * k6_x;
		particles[i].y += 37.0/378.0 * k1_y + 250.0/621.0 * k3_y + 125.0/594.0 * k4_y + 512.0/1771.0 * k6_y;

		particles[i].vx += 37.0/378.0 * k1_vx + 250.0/621.0 * k3_vx + 125.0/594.0 * k4_vx + 512.0/1771.0 * k6_vx;
		particles[i].vy += 37.0/378.0 * k1_vy + 250.0/621.0 * k3_vy + 125.0/594.0 * k4_vy + 512.0/1771.0 * k6_vy;
	}

	move();
}

void move(void) {
	// remove escaped particles
	for (unsigned int i = 0; i < local_number_of_particles; ++i) {
		if (particles[i].get_squared_distance_to_star() > pow2(parameters::particle_escape_radius)) {
			particles[i] = particles[local_number_of_particles-1];
			local_number_of_particles--;
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
		if (particles_size < local_number_of_particles+number) {
			particles_size = local_number_of_particles+number;
			particles = (t_particle*)realloc(particles, sizeof(t_particle)*particles_size);
		}

		MPI_Recv(particles+local_number_of_particles, number, mpi_particle, CPU_Next, 0, MPI_COMM_WORLD, &status);

		local_number_of_particles += number;
	}

	// check if we need to send particles to inner node
	int inward_offset[local_number_of_particles];
	int inward_size[local_number_of_particles];
	int inward_count = 0;
	for (unsigned int i = 0; i < local_number_of_particles; ++i) {
		if ((CPU_Rank > 0) && (particles[i].get_squared_distance_to_star() < local_r_min_squared)) {
			inward_offset[inward_count] = i;
			inward_size[inward_count] = 1;
			inward_count++;
		}
	}

	// send particles to inner node
	if (CPU_Rank > 0) {
		MPI_Datatype inward_type;
		MPI_Type_indexed(inward_count, inward_size, inward_offset, mpi_particle, &inward_type);
		MPI_Type_commit(&inward_type);
		MPI_Send(particles, 1, inward_type, CPU_Prev, 0, MPI_COMM_WORLD);
		MPI_Type_free(&inward_type);
	}

	// delete particles sent to inner node
	for (int i = 0; i < inward_count; ++i) {
		// swap particle with last particle
		particles[inward_offset[i]] = particles[local_number_of_particles-1];

		// check if the particle we switched is also in inward_offset
		for (int j = i+1; j < inward_count; ++j) {
			if (inward_offset[j] == (int)local_number_of_particles-1)
				inward_offset[j] = inward_offset[i];
		}

		local_number_of_particles--;
	}

	// move particles outwards starting from the inner most node

	// receive particles from outer node first
	if (CPU_Rank > 0) {
		MPI_Status status;
		MPI_Probe(CPU_Prev, 0, MPI_COMM_WORLD, &status);
		int number;
		MPI_Get_count(&status, mpi_particle, &number);

		// check if array is large enough for new particles
		if (particles_size < local_number_of_particles+number) {
			particles_size = local_number_of_particles+number;
			particles = (t_particle*)realloc(particles, sizeof(t_particle)*particles_size);
		}

		MPI_Recv(particles+local_number_of_particles, number, mpi_particle, CPU_Prev, 0, MPI_COMM_WORLD, &status);

		local_number_of_particles += number;
	}

	// check if we need to send particles to outer node
	int outward_offset[local_number_of_particles];
	int outward_size[local_number_of_particles];
	int outward_count = 0;
	for (unsigned int i = 0; i < local_number_of_particles; ++i) {
		if ((CPU_Rank < CPU_Highest) && (particles[i].get_squared_distance_to_star() > local_r_max_squared)) {
			outward_offset[outward_count] = i;
			outward_size[outward_count] = 1;
			outward_count++;
		}
	}

	// send particles to outer node
	if (CPU_Rank < CPU_Highest) {
		MPI_Datatype outward_type;
		MPI_Type_indexed(outward_count, outward_size, outward_offset, mpi_particle, &outward_type);
		MPI_Type_commit(&outward_type);
		MPI_Send(particles, 1, outward_type, CPU_Next, 0, MPI_COMM_WORLD);
		MPI_Type_free(&outward_type);
	}

	// delete particles sent to inner node
	for (int i = 0; i < outward_count; ++i) {
		// swap particle with last particle
		particles[outward_offset[i]] = particles[local_number_of_particles-1];

		// check if the particle we switched is also in outward_offset
		for (int j = i+1; j < outward_count; ++j) {
			if (outward_offset[j] == (int)local_number_of_particles-1)
				outward_offset[j] = outward_offset[i];
		}

		local_number_of_particles--;
	}

	// update global_number_of_particles
	MPI_Allreduce(&local_number_of_particles, &global_number_of_particles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

void write(unsigned int timestep) {
	MPI_File fh;
	MPI_Status status;
	int error, error_class, error_length;
	char *filename, error_string[MPI_MAX_ERROR_STRING+1];

	if (asprintf(&filename, "%s/particles%i.dat",OUTPUTDIR,timestep) <0) {
		die("Not enough memory!");
	}

	// try to open file
	error = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
	if (error != MPI_SUCCESS) {
		logging::print_master(LOG_ERROR "Error while writing to file '%s'. Check file permissions and IO support of MPI library\n", filename);

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
	unsigned int nodes_number_of_particles[CPU_Number];
	MPI_Allgather(&local_number_of_particles, 1, MPI_UNSIGNED, nodes_number_of_particles, 1, MPI_UNSIGNED, MPI_COMM_WORLD);

	// compute local offset
	unsigned int local_offset = 0;
	for (int cpu = 0; cpu < CPU_Rank; ++cpu) {
		local_offset += nodes_number_of_particles[cpu];
	}

	MPI_File_set_view(fh, 0, mpi_particle, mpi_particle, const_cast<char*>("native"), MPI_INFO_NULL);
	MPI_File_seek(fh, local_offset, MPI_SEEK_SET);
	MPI_File_write(fh, particles, local_number_of_particles, mpi_particle, &status);

	// close file
	MPI_File_close(&fh);
}

}
