/**
   \file selfgravity.cpp

   This file contains all functions for the calculation of SelfGravity.
   Previously they we're spreaded over several sg*.c files.

   All references to pages and equations are in "Toward predictive scenarios of
   planetary migration" by C. Baruteau
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef DISABLE_FFTW
#endif

#include <math.h>

#include "LowTasks.h"
#include "Theo.h"
#include "axilib.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "parameters.h"
#include "quantities.h"
#include "selfgravity.h"

#ifdef DISABLE_FFTW

namespace selfgravity {

	static void die_not_compiled() {
		logging::print_master(LOG_ERROR
			"Self-gravity is not compiled in. Please recompile with FFTW enabled.\n");
		PersonalExit(1);
	}

	void init([[maybe_unused]] t_data &data) { die_not_compiled(); }
	void compute([[maybe_unused]] t_data &data, [[maybe_unused]] double dt, [[maybe_unused]] bool update) { die_not_compiled(); }
	void init_azimuthal_velocity([[maybe_unused]] t_polargrid&) { die_not_compiled(); }
	void mpi_init(void) {}
	void mpi_finalize(void) {}

	double *g_radial;
	double *g_azimuthal;
}


#else // DISABLE_FFTW

#include <fftw3-mpi.h>

#ifndef NDEBUG
#undef FFTW_MEASURE
#define FFTW_MEASURE FFTW_ESTIMATE
#endif

namespace selfgravity
{
/// Parameters for selfgravity mode == Baruteau
/// smoothing parameter of potential (introduced as B on page 53)
double epsilon;
/// Parameters for selfgravity mode == Moldenhauer
/// smoothing parameters of potential as introduced by Tobi (to be published)
double lambda_sq;
double chi_sq;
/// Parameters for selfgravity mode == bessel kernel
double aspect_ratio;

/// mesh size in radial direction (introduced as Δu on page 53)
double r_step;
/// mesh size in azimuthal direction (introduced as Δphi on page 53)
double t_step;

/// plan for 2D FFT forward computation
fftw_plan fftplan_forward_K_radial;
fftw_plan fftplan_forward_K_azimuthal;
fftw_plan fftplan_forward_S_radial;
fftw_plan fftplan_forward_S_azimuthal;

/// plan for 2D FFT backward computation
fftw_plan fftplan_backward_acc_radial;
fftw_plan fftplan_backward_acc_azimuthal;

/* Arrays and fft-arrays, since (mpi)fft computes in-place transforms... */
double *K_radial;
double *K_azimuthal;
fftw_complex *FFT_K_radial;
fftw_complex *FFT_K_azimuthal;

double *S_radial;
double *S_azimuthal;
fftw_complex *FFT_S_radial;
fftw_complex *FFT_S_azimuthal;

double *acc_radial;
double *acc_azimuthal;
fftw_complex *FFT_acc_radial;
fftw_complex *FFT_acc_azimuthal;

// Arrays for use outside the sg module.
// Memory is allocated in the data construct.
// These variables only point to the respective Field array.
double *g_radial;
double *g_azimuthal;



void mpi_init(void)
{
    // init MPI
    fftw_mpi_init();
}

void mpi_finalize(void)
{
    // destroy plans
    fftw_destroy_plan(fftplan_forward_K_radial);
    fftw_destroy_plan(fftplan_forward_K_azimuthal);
    fftw_destroy_plan(fftplan_forward_S_radial);
    fftw_destroy_plan(fftplan_forward_S_azimuthal);
    fftw_destroy_plan(fftplan_backward_acc_radial);
    fftw_destroy_plan(fftplan_backward_acc_azimuthal);

    // free memory
    fftw_free(K_radial);
    fftw_free(K_azimuthal);
    fftw_free(FFT_K_radial);
    fftw_free(FFT_K_azimuthal);
    fftw_free(S_radial);
    fftw_free(S_azimuthal);
    fftw_free(FFT_S_radial);
    fftw_free(FFT_S_azimuthal);
    fftw_free(acc_radial);
    fftw_free(acc_azimuthal);
    fftw_free(FFT_acc_radial);
    fftw_free(FFT_acc_azimuthal);

    // cleanup MPI
#ifdef _OPENMP
	fftw_cleanup_threads();
#else
	fftw_mpi_cleanup();
#endif // _OPENMP
}

/**
 * @brief Obtain the aspect ratio of the disk and recalculate the average if needed.
 */
static double get_aspect_ratio(t_data &data) {
	double rv;
	if (parameters::Locally_Isothermal) {
		rv = parameters::aspectratio_ref;
	} else {
		rv = quantities::gas_allreduce_mass_average(data, data[t_data::ASPECTRATIO], RMAX);
	}
	// safety net
	if (rv == 0.0) {
		rv = parameters::aspectratio_ref;
	}
	return rv;
}

/**
 * @brief update_sg_constants, calculate mass weighted aspect ratio and use it
 * to update lambda_sq and chi_sq. see Moldenhauer 2018 for more info on the
 * constants or Baruteau, 2008 for info on the self gravity module.
 * @param data
 */
void update_sg_constants() {
	if (parameters::self_gravity_mode == parameters::t_sg::sg_S) {
		lambda_sq = std::pow(
			0.4571 * aspect_ratio + 0.6737 * std::sqrt(aspect_ratio), 2);
		chi_sq = std::pow((-0.7543 * aspect_ratio + 0.6472) * aspect_ratio, 2);
	} else {
		epsilon = parameters::thickness_smoothing_sg*aspect_ratio;
	}
}


/**
 * @brief upate_kernel, update the self-gravity kernel if needed.
 * @param data
 */
static void update_kernel(t_data &data)
{
	if (parameters::Locally_Isothermal) {
		return;
	}

	// only update every Nth timestep
    const int update_every_nth_step = parameters::self_gravity_steps_between_kernel_update;
    static int since_last_calculated = update_every_nth_step;

	if (since_last_calculated >= (update_every_nth_step - 1)) {
		since_last_calculated = 0;
	} else {
		since_last_calculated++;
		return;
	}

	// only update if the aspect ratio changed enough
	static double last_aspect_ratio = 0;
	aspect_ratio = get_aspect_ratio(data);

	if (std::abs(last_aspect_ratio - aspect_ratio) < parameters::self_gravity_aspectratio_change_threshold) {
		return;
	}
	last_aspect_ratio = aspect_ratio;

	update_sg_constants();
	compute_FFT_kernel();
}

/**
   Initializes self gravity.
*/
void init(t_data &data)
{
    // Self-gravity must be on a polar logarithmic grid
    if (parameters::radial_grid_type != parameters::logarithmic_spacing) {
	logging::print_master(
	    LOG_ERROR
	    "A logarithmic grid is needed to compute self-gravity with polar method. Try again!\n");
	PersonalExit(1);
    }



#ifdef _OPENMP
	int thread_success = fftw_init_threads();
	// Tell plans to use openmp
	// this must be done before creating plans
	if(thread_success != 0){
	int num_openmp_threads;
	#pragma omp parallel
	{
	num_openmp_threads = omp_get_num_threads();
	}

	fftw_plan_with_nthreads(num_openmp_threads);
	}
#endif

    r_step = std::log(Radii[GlobalNRadial] / Radii[0]) / (double)GlobalNRadial;
    t_step = 2.0 * M_PI / (double)NAzimuthal;

    // allocate memory
    K_radial = fftw_alloc_real(2 * total_local_size);
    K_azimuthal = fftw_alloc_real(2 * total_local_size);
    FFT_K_radial = fftw_alloc_complex(total_local_size);
    FFT_K_azimuthal = fftw_alloc_complex(total_local_size);
    S_radial = fftw_alloc_real(2 * total_local_size);
    S_azimuthal = fftw_alloc_real(2 * total_local_size);
    FFT_S_radial = fftw_alloc_complex(total_local_size);
    FFT_S_azimuthal = fftw_alloc_complex(total_local_size);
    acc_radial = fftw_alloc_real(2 * total_local_size);
    acc_azimuthal = fftw_alloc_real(2 * total_local_size);
    FFT_acc_radial = fftw_alloc_complex(total_local_size);
    FFT_acc_azimuthal = fftw_alloc_complex(total_local_size);

    g_radial = data[t_data::SG_ACCEL_RAD].Field;
    g_azimuthal = data[t_data::SG_ACCEL_AZI].Field;

    // create FFT plans
    fftplan_forward_K_radial = fftw_mpi_plan_dft_r2c_2d(
	2 * GlobalNRadial, NAzimuthal, K_radial, FFT_K_radial, MPI_COMM_WORLD,
	FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);
    fftplan_forward_K_azimuthal = fftw_mpi_plan_dft_r2c_2d(
	2 * GlobalNRadial, NAzimuthal, K_azimuthal, FFT_K_azimuthal,
	MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);
    fftplan_forward_S_radial = fftw_mpi_plan_dft_r2c_2d(
	2 * GlobalNRadial, NAzimuthal, S_radial, FFT_S_radial, MPI_COMM_WORLD,
	FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);
    fftplan_forward_S_azimuthal = fftw_mpi_plan_dft_r2c_2d(
	2 * GlobalNRadial, NAzimuthal, S_azimuthal, FFT_S_azimuthal,
	MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);

    fftplan_backward_acc_radial = fftw_mpi_plan_dft_c2r_2d(
	2 * GlobalNRadial, NAzimuthal, FFT_acc_radial, acc_radial,
	MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN);
    fftplan_backward_acc_azimuthal = fftw_mpi_plan_dft_c2r_2d(
	2 * GlobalNRadial, NAzimuthal, FFT_acc_azimuthal, acc_azimuthal,
	MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN);

	aspect_ratio = get_aspect_ratio(data);
	update_sg_constants();
	compute_FFT_kernel();
    compute(data, 0, false);

	// check that the ratio of outer and inner boundary radius
	// is compatible with the coefficients for the smoothing length
	if (parameters::self_gravity_mode == parameters::t_sg::sg_S) {
		const double radius_ratio = RMAX/RMIN;
		if (std::abs(radius_ratio - 12.5) > 2) {
			logging::print_master(LOG_WARNING "WARNING: The ratio of outer to inner boundary radius (%f) is substaintially different than 12.5. Consider recomputing the factors for the polynomials to calculate the smoothing length. Refer to section 2.2 page 24 of Tobias Moldenhauer's Master thesis.\n", radius_ratio);
		}
	}

}

void compute(t_data &data, double dt, bool update)
{
    update_kernel(data);

    auto &density = data[t_data::SIGMA];
    compute_FFT_density(density);
    compute_acceleration(density);

    if (update) {
	// Computes polar components of acceleration and updates values of vrad,
	// vazi at each step

	auto &v_radial = data[t_data::V_RADIAL];
	auto &v_azimuthal = data[t_data::V_AZIMUTHAL];
	update_velocities(v_radial, v_azimuthal, dt);
    }
}

void compute_FFT_density(t_polargrid &density)
{

    MPI_Request req;
	double *dens = density.Field;
    const unsigned int stride = 2 * (NAzimuthal / 2 + 1);

    // We communicate the hydro density field to the fftw domain decomposition
    // (d. d.)
    const int one_if_odd = (CPU_Number % 2 == 0 ? 0 : 1);

    // every cpu except that one with no 'friend' needs to interchange data
    if (CPU_Rank != CPU_NoFriend) {
	if (CPU_Rank >= (CPU_Number + one_if_odd) / 2) {
	    /* all upper cpus send data */
	    MPI_Isend(&dens[Zero_or_active * NAzimuthal],
		      active_hydro_totalsize, MPI_DOUBLE, CPU_Friend, 30,
		      MPI_COMM_WORLD, &req);
	} else {
	    /* all lower cpus recieve data */
	    MPI_Irecv(&dens_friend[0], active_hydro_totalsize_friend,
		      MPI_DOUBLE, CPU_Friend, 30, MPI_COMM_WORLD, &req);
	}
	MPI_Wait(&req, &global_MPI_Status);
    }

    // lower half of cpus
    if ((CPU_Rank < CPU_Number / 2) && (CPU_Rank != CPU_NoFriend)) {
	#pragma omp parallel for collapse(2)
	for (unsigned int i = 0; i < (unsigned int) ifront + 1; i++) {
	    for (unsigned int j = 0; j < NAzimuthal; j++) {
		const unsigned int l = i * stride + j;
		const unsigned int ih = i + Zero_or_active;
		const unsigned int lh = ih * NAzimuthal + j;
		/* S_r = sigma(u,phi) exp(u/2) */
		S_radial[l] = dens[lh] * std::sqrt(Rmed[ih] / GlobalRmed[0]);
		/* S_t = sigma(u,phi) exp(3*u/2) */
		S_azimuthal[l] = S_radial[l] * Rmed[ih] / GlobalRmed[0];
	    }
	}

	#pragma omp parallel for collapse(2)
	for (unsigned int i = (unsigned int) ifront + 1; i < (unsigned int)local_Nx; i++) {
	    for (unsigned int j = 0; j < NAzimuthal; j++) {
		const unsigned int l = i * stride + j;
		const unsigned int ih = i - ((unsigned int) ifront + 1);
		const unsigned int lh = ih * NAzimuthal + j;
		const unsigned int ir = i + IMIN + Zero_or_active;
		if ((i + local_i_start) < GlobalNRadial) {
		    S_radial[l] = dens_friend[lh] * std::sqrt(GlobalRmed[ir] / GlobalRmed[0]);
		    S_azimuthal[l] = S_radial[l] * GlobalRmed[ir] / GlobalRmed[0];
		} else {
		    S_radial[l] = 0.;
		    S_azimuthal[l] = 0.;
		}
	    }
	}
    }

    // cpu with no friend (upper most if odd total cpu number)
    if (CPU_Rank == CPU_NoFriend) {
	#pragma omp parallel for collapse(2)
	for (unsigned int i = 0; i < (unsigned int)local_Nx; i++) {
	    for (unsigned int j = 0; j < NAzimuthal; j++) {
		const unsigned int l = i * stride + j;
		const unsigned int ih = i + Zero_or_active;
		if ((i + local_i_start) < GlobalNRadial) {
		    unsigned int lh = ih * NAzimuthal + j;
		    S_radial[l] = dens[lh] * std::sqrt(Rmed[ih] / GlobalRmed[0]);
		    S_azimuthal[l] = S_radial[l] * Rmed[ih] / GlobalRmed[0];
		} else {
		    S_radial[l] = 0.;
		    S_azimuthal[l] = 0.;
		}
	    }
	}
    }

    // upper half of cpus
    if ((CPU_Rank >= CPU_Number / 2) && (CPU_Rank != CPU_NoFriend)) {
	#pragma omp parallel for collapse(2)
	for (unsigned int i = 0; i < (unsigned int)local_Nx; i++) {
	    for (unsigned int j = 0; j < NAzimuthal; j++) {
		const unsigned int l = i * stride + j;
		if ((i + local_i_start) >= GlobalNRadial) {
		    S_radial[l] = 0.;
		    S_azimuthal[l] = 0.;
		}
	    }
	}
    }

    // Now we can compute the in-place ffts of reduced density arrays.
    fftw_execute(fftplan_forward_S_radial);
    fftw_execute(fftplan_forward_S_azimuthal);
}

void compute_FFT_kernel()
{
    const unsigned int stride = 2 * (NAzimuthal / 2 + 1);

	#pragma omp parallel for
    for (unsigned int i = 0; i < (unsigned int)local_Nx; i++) {
    double u;
	if (i + local_i_start < GlobalNRadial) {
	    u = std::log(Radii[i + local_i_start] / Radii[0]);
	} else {
	    u = -std::log(Radii[2 * GlobalNRadial - (i + local_i_start)] /
			  Radii[0]);
	}

	for (unsigned int j = 0; j < NAzimuthal; j++) {
	    const unsigned int l = i * stride + j;
		const double theta = dphi * (double)j;

	if (parameters::self_gravity_mode == parameters::t_sg::sg_B) {
		const double denominator = std::pow(epsilon*epsilon*std::exp(u)
									+ 2.0 * (std::cosh(u) - std::cos(theta)),-1.5);

		K_radial[l] = 1.0 + epsilon*epsilon - std::cos(theta) * std::exp(-u);
		K_radial[l] *= denominator;
		K_azimuthal[l] = std::sin(theta);
		K_azimuthal[l] *= denominator;
	} else if (parameters::self_gravity_mode == parameters::t_sg::sg_BK) {

		/* This Kernel is an anlytical solution in the limit Q->oo
		* and can be faithfully used for Q>=20.
		* A further correction accounting for all Q values is coming 
		* at end 2023.
		*/
		if (u==0. && theta==0.) {
			/* At the singularity we must cancel the Kernels or 
			* use tapering functions (see Sect. 4.1 of 
			* https://doi.org/10.1051/0004-6361/202346178).
			*/
			K_radial[l] = 0.;
			K_azimuthal[l]  = 0.;

		} else {
			const double distance_squared = 2. * std::pow(aspect_ratio, -2.) 
												* (std::cosh(u) - std::cos(theta)) 
												/  std::cosh(u);

			/* The modified Bessel functions of the second kind are also 
			* known as "Irregular modified cylindrical Bessel functions".
			* This last naming is the one used in the cmath library
			*/

			/* L_sg would be defined in Rendon Restrepo et al. 2024 (not yet published) */
			/* The following computation of L_sg solve an exponential overflow/underflow for large distances
			* by using the exact Taylor expansion at infinite of the Bessel Kernel . 
			*/

			double L_sg;
            if (X_aux < 60) {
                L_sg = std::pow(M_PI, 0.5)
                     * X_aux
                     * std::exp(X_aux)
                     * ( std::cyl_bessel_kl(1., X_aux)
                     - std::cyl_bessel_kl(0., X_aux) );
            } else { // Taylor expansion at inifinity in order to avoid exp overflow
                L_sg = std::pow(M_PI, 0.5)
                     * X_aux
                     * 0.5 * std::pow(M_PI/2., 0.5)
                     * ( std::pow(X_aux, -1.5)
                       - 3./8.*std::pow(X_aux, -2.5)
                       + 45./128.*std::pow(X_aux, -3.5) );
            }
			
			K_radial[l] = (L_sg / 2. / M_PI / aspect_ratio)
						* std::pow(std::cosh(u), -0.5)
						* std::pow(std::cosh(u)-std::cos(theta), -1.)
						* (1.0-std::cos(theta)*std::exp(-u));

			K_azimuthal[l] =  (L_sg / 2. / M_PI / aspect_ratio)
							* std::pow(std::cosh(u), -0.5)
							* std::pow(std::cosh(u)-std::cos(theta), -1.)
							* std::sin(theta);
		}


	} else if (parameters::self_gravity_mode == parameters::t_sg::sg_S) {
		const double denominator = std::pow(
		2 * (std::cosh(u) - std::cos(theta)) +
		    lambda_sq * (std::exp(u) + std::exp(-u) - 2) + chi_sq, -1.5);

		K_radial[l] = 1.0 - std::cos(theta) * std::exp(-u);
		K_radial[l] *= denominator;
		K_azimuthal[l] = std::sin(theta);
		K_azimuthal[l] *= denominator;
	}
	}
    }

    fftw_execute(fftplan_forward_K_radial);
    fftw_execute(fftplan_forward_K_azimuthal);
}

void compute_acceleration(t_polargrid &density)
{
    MPI_Request req1, req3, req4, req5, req6;
    const unsigned int nr = density.Nrad;
    const unsigned int stride = 2 * (NAzimuthal / 2 + 1);
    const int one_if_odd = (CPU_Number % 2 == 0 ? 0 : 1);

    // First we compute sg_acc as a convolution product of reduced density and
    // kernel arrays. In fact, all bufffttabs are transposed arrays, since we
    // used the flag FFTW_TRANSPOSED_ORDER before. However, this is not a
    // problem: we compute double and imaginary parts of 1-dimension convolution
    // arrays buffer_FFT_acc_radial(t)

	#pragma omp parallel for
    for (unsigned int i = 0; i < total_local_size; i++) {
	FFT_acc_radial[i][0] =
	    -constants::G * (FFT_K_radial[i][0] * FFT_S_radial[i][0] -
			     FFT_K_radial[i][1] * FFT_S_radial[i][1]);
	FFT_acc_radial[i][1] =
	    -constants::G * (FFT_K_radial[i][0] * FFT_S_radial[i][1] +
			     FFT_K_radial[i][1] * FFT_S_radial[i][0]);

	FFT_acc_azimuthal[i][0] =
	    -constants::G * (FFT_K_azimuthal[i][0] * FFT_S_azimuthal[i][0] -
			     FFT_K_azimuthal[i][1] * FFT_S_azimuthal[i][1]);
	FFT_acc_azimuthal[i][1] =
	    -constants::G * (FFT_K_azimuthal[i][0] * FFT_S_azimuthal[i][1] +
			     FFT_K_azimuthal[i][1] * FFT_S_azimuthal[i][0]);
    }

    fftw_execute(fftplan_backward_acc_radial);
    fftw_execute(fftplan_backward_acc_azimuthal);

    // The use again of argument FFTW_TRANSPOSED_ORDER in these backward fourier
    // transforms ensures that arrays are not transposed anymore. Then we
    // transfer the exact necessary quantity of sg_acceleration arrays from the
    // fftw d.d. to the hydro mesh
    if (CPU_Rank != CPU_NoFriend) {
	if (CPU_Rank < CPU_Number / 2) {
		#pragma omp parallel for
	    for (unsigned int i = 0; i < transfer_size; i++) {
		if (i < transfer_size / 2)
		    ffttohydro_transfer[i] =
			acc_radial[((unsigned int) ifront + 1 - CPUOVERLAP) * stride + i];
		else
		    ffttohydro_transfer[i] =
			acc_azimuthal[((unsigned int) ifront + 1 - CPUOVERLAP) * stride + i -
				      transfer_size / 2];
	    }
	    MPI_Isend(ffttohydro_transfer, transfer_size, MPI_DOUBLE,
		      CPU_Friend, 40, MPI_COMM_WORLD, &req1);
	    MPI_Wait(&req1, &global_MPI_Status);
	} else {
	    MPI_Irecv(ffttohydro_transfer_friend, transfer_size_friend,
		      MPI_DOUBLE, CPU_Friend, 40, MPI_COMM_WORLD, &req1);
	    MPI_Wait(&req1, &global_MPI_Status);
	}
    }

    // We now compute sg_acceleration arrays on the hydro mesh
    if (CPU_Rank < (CPU_Number + one_if_odd) / 2) {
	if (CPU_Rank == 0) {
		#pragma omp parallel for collapse(2)
	    for (unsigned int i = 0; i < nr; i++) {
		for (unsigned int j = 0; j < NAzimuthal; j++) {
		    const unsigned int l = i * NAzimuthal + j;
		    g_radial[l] = acc_radial[i * stride + j];
		    g_azimuthal[l] = acc_azimuthal[i * stride + j];
		}
	    }
	} else {
		#pragma omp parallel for collapse(2)
	    for (unsigned int i = Zero_or_active; i < nr; i++) {
		for (unsigned int j = 0; j < NAzimuthal; j++) {
		    const unsigned int l = i * NAzimuthal + j;
		    g_radial[l] = acc_radial[(i - Zero_or_active) * stride + j];
		    g_azimuthal[l] =
			acc_azimuthal[(i - Zero_or_active) * stride + j];
		}
	    }
	}
    }

    if (CPU_Rank >= (CPU_Number + one_if_odd) / 2) {
	if (CPU_Rank == CPU_Highest) {
		#pragma omp parallel for collapse(2)
	    for (unsigned int i = 0; i < nr; i++) {
		for (unsigned int j = 0; j < NAzimuthal; j++) {
		    const unsigned int l = i * NAzimuthal + j;
		    g_radial[l] = ffttohydro_transfer_friend[i * stride + j];
		    g_azimuthal[l] =
			ffttohydro_transfer_friend[transfer_size_friend / 2 +
						   i * stride + j];
		}
	    }
	} else {
		#pragma omp parallel for collapse(2)
	    for (unsigned int i = 0; i < Max_or_active; i++) {
		for (unsigned int j = 0; j < NAzimuthal; j++) {
		    unsigned int l = i * NAzimuthal + j;
		    g_radial[l] = ffttohydro_transfer_friend[i * stride + j];
		    g_azimuthal[l] =
			ffttohydro_transfer_friend[transfer_size_friend / 2 +
						   i * stride + j];
		}
	    }
	}
    }

    // Now we exchange the correct amount of sg_acceleration between cpus to
    // fill ghosts.

    const unsigned int ghost_size = CPUOVERLAP * NAzimuthal;
    if (CPU_Number > 1) {
	if ((CPU_Rank > 0) && (CPU_Rank < (CPU_Number + one_if_odd) / 2)) {
	    MPI_Isend(&g_azimuthal[Zero_or_active * NAzimuthal], ghost_size,
		      MPI_DOUBLE, CPU_Prev, 60, MPI_COMM_WORLD, &req3);
	    MPI_Wait(&req3, &global_MPI_Status);
	    MPI_Irecv(&g_azimuthal[0], ghost_size, MPI_DOUBLE, CPU_Prev, 61,
		      MPI_COMM_WORLD, &req4);
	    MPI_Wait(&req4, &global_MPI_Status);
	}
	if ((CPU_Rank >= (CPU_Number + one_if_odd) / 2) &&
	    (CPU_Rank != CPU_Highest)) {
	    MPI_Irecv(&g_azimuthal[Max_or_active * NAzimuthal], ghost_size,
		      MPI_DOUBLE, CPU_Next, 60, MPI_COMM_WORLD, &req3);
	    MPI_Wait(&req3, &global_MPI_Status);
	    MPI_Isend(&g_azimuthal[(Max_or_active - CPUOVERLAP) * NAzimuthal],
		      ghost_size, MPI_DOUBLE, CPU_Next, 61, MPI_COMM_WORLD,
		      &req4);
	    MPI_Wait(&req4, &global_MPI_Status);
	}
	if ((CPU_Rank > 0) && (CPU_Rank < (CPU_Number + one_if_odd) / 2)) {
	    MPI_Isend(&g_radial[Zero_or_active * NAzimuthal], ghost_size,
		      MPI_DOUBLE, CPU_Prev, 50, MPI_COMM_WORLD, &req5);
	    MPI_Wait(&req5, &global_MPI_Status);
	    MPI_Irecv(&g_radial[0], ghost_size, MPI_DOUBLE, CPU_Prev, 51,
		      MPI_COMM_WORLD, &req6);
	    MPI_Wait(&req6, &global_MPI_Status);
	}
	if ((CPU_Rank >= (CPU_Number + one_if_odd) / 2) &&
	    (CPU_Rank != CPU_Highest)) {
	    MPI_Irecv(&g_radial[Max_or_active * NAzimuthal], ghost_size,
		      MPI_DOUBLE, CPU_Next, 50, MPI_COMM_WORLD, &req5);
	    MPI_Wait(&req5, &global_MPI_Status);
	    MPI_Isend(&g_radial[(Max_or_active - CPUOVERLAP) * NAzimuthal],
		      ghost_size, MPI_DOUBLE, CPU_Next, 51, MPI_COMM_WORLD,
		      &req6);
	    MPI_Wait(&req6, &global_MPI_Status);
	}
    }

    // We don't forget to renormalize acc arrays!
	#pragma omp parallel for
    for (unsigned int i = 0; i < nr; i++) {
	// g_r(u,phi) normalized with exp(-u/2)*Δu*Δphi/(2*N_r*N_phi) (3.43 page
	// 57) g_phi(u,phi) normalized with exp(-3*u/2)*Δu Δphi/(2*N_r*N_phi)
	// (3.44 page 57)
	double normaccr = r_step * t_step /
		   ((double)(2 * GlobalNRadial) * (double)NAzimuthal);
	double normacct = normaccr;
	normaccr /= std::sqrt(Rmed[i] / GlobalRmed[0]);
	normacct /=
	    (Rmed[i] / GlobalRmed[0] * std::sqrt(Rmed[i] / GlobalRmed[0]));
	for (unsigned int j = 0; j < NAzimuthal; j++) {
	    const unsigned int l = i * NAzimuthal + j;
	    g_radial[l] *= normaccr;
	    g_azimuthal[l] *= normacct;
	}
    }
#ifdef EPSILON_SMOOTHING_SG
	// Eventually, we take the compensation from selfforce into account
	// g_r(u,phi) is corrected by G*sigma(u,phi)*Δu*Δphi/B (3.43 page 57)
	#pragma omp parallel for
	for (unsigned int i = 0 ; i < nr; i++ ) {
		for (unsigned int j = 0; j < NAzimuthal; j++ ) {
			const unsigned int l = i*NAzimuthal + j;
			if ( (i+IMIN) < GlobalNRadial ) {
				g_radial[l] += constants::G*density.Field[l]*r_step*t_step / epsilon;
			}
		}
	}
#endif
}

/**
   Update the velocity fields to take into account self-gravity

   \param VRad
   \param Vazi
   \param Dt
*/
void update_velocities(t_polargrid &v_radial, t_polargrid &v_azimuthal,
		       double dt)
{
	const unsigned int Nr = v_radial.get_max_radial();

    // Here we update velocity fields to take into account self-gravity
	#pragma omp parallel for collapse(2)
	for (unsigned int n_radial = 0; n_radial < Nr;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0; n_azimuthal < NAzimuthal;
	     ++n_azimuthal) {
	    const unsigned int l = n_radial * NAzimuthal + n_azimuthal;
	    // We compute VRadial - half-centered in azimuth - from
	    // centered-in-cell radial sg acceleration
	    if (n_radial > 0) {
		v_radial(n_radial, n_azimuthal) +=
		    dt *
		    ((Rinf[n_radial] - Rmed[n_radial - 1]) * g_radial[l] +
		     (Rmed[n_radial] - Rinf[n_radial]) *
			 g_radial[l - NAzimuthal]) *
		    InvDiffRmed[n_radial];
	    }

	    // We compute VAzimuthal - half-centered in radius - from
	    // centered-in-cell azimutal sg acceleration
		unsigned int jm1;
	    if (n_azimuthal == 0) {
			jm1 = NAzimuthal - 1;
		} else {
			jm1 = n_azimuthal - 1;
		}
	    const unsigned int lm1 = n_radial * NAzimuthal + jm1;
	    v_azimuthal(n_radial, n_azimuthal) += 0.5 * dt * (g_azimuthal[l] + g_azimuthal[lm1]);
	}
    }
}

void init_azimuthal_velocity(t_polargrid &v_azimuthal)
{
    double *GLOBAL_AxiSGAccr = (double *)malloc(sizeof(double) * GlobalNRadial);
    mpi_make1Dprofile(g_radial, GLOBAL_AxiSGAccr);

	#pragma omp parallel for
    for (unsigned int n_radial = 0;
	 n_radial < v_azimuthal.get_size_radial() - GHOSTCELLS_B; ++n_radial) {

	// this corresponds to equation (3.42) in Baruteau, 2008
		// also considering derivative of smoothed potential
    double omega_cell;
    if(parameters::v_azimuthal_with_quadropole_support){
    omega_cell = initial_locally_isothermal_smoothed_v_az_with_quadropole_moment(Rmed[n_radial], hydro_center_mass) / Rmed[n_radial];
    } else {
    omega_cell = initial_locally_isothermal_smoothed_v_az(Rmed[n_radial], hydro_center_mass) / Rmed[n_radial];
    }
	double temp = std::pow(omega_cell, 2) - GLOBAL_AxiSGAccr[n_radial + IMIN] / Rmed[n_radial];
	if (temp < 0) {
	    logging::print(
		"Radicand %lg < 0 in init_azimuthal_velocity! Maybe ThicknessSmoothingSG (%lg) is too small!\n",
		temp, parameters::thickness_smoothing_sg);
	}
	const double omega = std::sqrt(temp);

	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < v_azimuthal.get_size_azimuthal(); ++n_azimuthal) {
	    v_azimuthal(n_radial, n_azimuthal) = Rmed[n_radial] * omega;
	}
    }

    free(GLOBAL_AxiSGAccr);
}

} // namespace selfgravity


#endif // DISABLE_FFTW