/**
	\file selfgravity.cpp

	This file contains all functions for the calculation of SelfGravity. Previously they we're spreaded over several sg*.c files.

	All references to pages and equations are in "Toward predictive scenarios of planetary migration" by C. Baruteau
*/

#include <math.h>
#include <fftw3-mpi.h>

#include "Theo.h"
#include "util.h"
#include "parameters.h"
#include "selfgravity.h"
#include "constants.h"
#include "LowTasks.h"
#include "axilib.h"
#include "time.h"
#include "global.h"
#include "logging.h"
#include "quantities.h"

#ifndef NDEBUG
#undef FFTW_MEASURE
#define FFTW_MEASURE FFTW_ESTIMATE
#endif


namespace selfgravity {

/// Initialization counter
static int init_counter = 0;

double lambda_sq;
double chi_sq;
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

double* acc_radial;
double* acc_azimuthal;
fftw_complex *FFT_acc_radial;
fftw_complex *FFT_acc_azimuthal;

double* g_radial;
double* g_azimuthal;

void mpi_init(void)
{
	// init MPI
	fftw_mpi_init();
}

void mpi_finalize(void)
{
	if (init_counter > 0) {
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

		free(g_radial);
		free(g_azimuthal);
	}

	// cleanup MPI
	fftw_mpi_cleanup();
}


/**
 * @brief update_tobi_constants, calculate mass weighted aspect ratio and use it to update lambda_sq and chi_sq.
 * see Tobi's master thesis for more info on the constants or Baruteau, 2008 for info on the self gravity module.
 * @param data
 * @param timestep
 * @param force_update
 */
void update_tobi_constants(t_data &data)
{
	// only update every 10th timestep
	constexpr int update_every_nth_step = 10;
	static int last_timestep_calculated = update_every_nth_step;
	last_timestep_calculated++;
	if(last_timestep_calculated < update_every_nth_step)
		return;
	last_timestep_calculated = 0;


	double aspect_ratio = quantities::gas_aspect_ratio(data);

	lambda_sq = pow2(0.4571*aspect_ratio + 0.6737*sqrt(aspect_ratio));
	chi_sq = pow2((-0.7543*aspect_ratio + 0.6472)*aspect_ratio);

	return;
}

/**
	Initializes self gravity.
*/
void init()
{
	// Self-gravity must be on a polar logarithmic grid
	if ( parameters::radial_grid_type != parameters::logarithmic_spacing ) {
		logging::print_master(LOG_ERROR "A logarithmic grid is needed to compute self-gravity with polar method. Try again!\n");
		PersonalExit(1);
	}

	lambda_sq = pow2(0.4571*ASPECTRATIO_REF + 0.6737*sqrt(ASPECTRATIO_REF));
	chi_sq = pow2((-0.7543*ASPECTRATIO_REF + 0.6472)*ASPECTRATIO_REF);
	r_step = log(Radii[GlobalNRadial]/Radii[0]) / (double)GlobalNRadial;
	t_step = 2.0*PI/(double)NAzimuthal;

	// allocate memory
	K_radial = fftw_alloc_real(2*total_local_size);
	K_azimuthal = fftw_alloc_real(2*total_local_size);
	FFT_K_radial = fftw_alloc_complex(total_local_size);
	FFT_K_azimuthal = fftw_alloc_complex(total_local_size);
	S_radial = fftw_alloc_real(2*total_local_size);
	S_azimuthal = fftw_alloc_real(2*total_local_size);
	FFT_S_radial = fftw_alloc_complex(total_local_size);
	FFT_S_azimuthal = fftw_alloc_complex(total_local_size);
	acc_radial = fftw_alloc_real(2*total_local_size);
	acc_azimuthal = fftw_alloc_real(2*total_local_size);
	FFT_acc_radial = fftw_alloc_complex(total_local_size);
	FFT_acc_azimuthal = fftw_alloc_complex(total_local_size);

	g_radial = (double*)malloc(sizeof(double)*hydro_totalsize);
	g_azimuthal = (double*)malloc(sizeof(double)*hydro_totalsize);

	// create FFT plans
	fftplan_forward_K_radial = fftw_mpi_plan_dft_r2c_2d(2*GlobalNRadial, NAzimuthal, K_radial, FFT_K_radial, MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);
	fftplan_forward_K_azimuthal = fftw_mpi_plan_dft_r2c_2d(2*GlobalNRadial, NAzimuthal, K_azimuthal, FFT_K_azimuthal, MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);
	fftplan_forward_S_radial = fftw_mpi_plan_dft_r2c_2d(2*GlobalNRadial, NAzimuthal, S_radial, FFT_S_radial, MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);
	fftplan_forward_S_azimuthal = fftw_mpi_plan_dft_r2c_2d(2*GlobalNRadial, NAzimuthal, S_azimuthal, FFT_S_azimuthal, MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);

	fftplan_backward_acc_radial = fftw_mpi_plan_dft_c2r_2d(2*GlobalNRadial, NAzimuthal, FFT_acc_radial, acc_radial, MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN);
	fftplan_backward_acc_azimuthal = fftw_mpi_plan_dft_c2r_2d(2*GlobalNRadial, NAzimuthal, FFT_acc_azimuthal, acc_azimuthal, MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN);
}

void compute(t_polargrid &density, t_polargrid &v_radial, t_polargrid &v_azimuthal, double dt, bool update)
{
	if (init_counter == 0) {
		init();
		compute_FFT_kernel();
	}

	compute_FFT_density(density);

	// Here we compute radial and azimutal components of sg acceleration as a convolution product of reduced density and kernel arrays
	compute_acceleration(density);

	if (update) {
		// Computes polar components of acceleration and updates values of vrad, vtheta at each step
		update_velocities(v_radial, v_azimuthal, dt);
	}

	init_counter++;
}

void compute_FFT_density(t_polargrid &density)
{

	MPI_Request req;
	unsigned int i, j;
	int l;
	int ih, lh, ir;
	int one_if_odd, stride;
	double *dens;
	dens = density.Field;
	stride = 2*(NAzimuthal/2+1);

	// We communicate the hydro density field to the fftw domain decomposition (d. d.)
	one_if_odd = (CPU_Number%2 == 0 ? 0 : 1);

	// every cpu except that one with no 'friend' needs to interchange data
	if ( CPU_Rank != CPU_NoFriend ) {
		if ( CPU_Rank >= (CPU_Number+one_if_odd)/2 ) {
			/* all upper cpus send data */
			MPI_Isend(&dens[Zero_or_active*NAzimuthal], active_hydro_totalsize, MPI_DOUBLE, CPU_Friend, 30, MPI_COMM_WORLD, &req);
		} else {
			/* all lower cpus recieve data */
			MPI_Irecv (&dens_friend[0], active_hydro_totalsize_friend, MPI_DOUBLE, CPU_Friend, 30, MPI_COMM_WORLD, &req);
		}
        MPI_Wait (&req, &global_MPI_Status);
	}

	// lower half of cpus
	if ( (CPU_Rank < CPU_Number/2) && (CPU_Rank!=CPU_NoFriend) ) {
		for ( i = 0; i < ifront+1; i++ ) {
			for ( j = 0; j < NAzimuthal; j++ ) {
				l = i*stride + j;
				ih = i+Zero_or_active;
				lh = ih*NAzimuthal + j;
				/* S_r = sigma(u,phi) exp(u/2) */
				S_radial[l] = dens[lh] * sqrt(Rmed[ih] / GlobalRmed[0]);
				/* S_t = sigma(u,phi) exp(3*u/2) */
				S_azimuthal[l] = S_radial[l] * Rmed[ih] / GlobalRmed[0];
			}
		}

		for (i = ifront+1; i < (unsigned int)local_Nx; i++ ) {
			for ( j = 0; j < NAzimuthal; j++ ) {
				l = i*stride + j;
				ih = i-(ifront+1);
				lh = ih*NAzimuthal + j;
				ir = i+IMIN+Zero_or_active;
				if ( (i+local_i_start) < GlobalNRadial ) {
					S_radial[l] = dens_friend[lh] * sqrt(GlobalRmed[ir] / GlobalRmed[0]);
					S_azimuthal[l] = S_radial[l] * GlobalRmed[ir] / GlobalRmed[0];
				} else {
					S_radial[l] = 0.;
					S_azimuthal[l] = 0.;
				}
			}
		}
	}

	// cpu with no friend (upper most if odd total cpu number)
	if ( CPU_Rank == CPU_NoFriend ) {
		for ( i = 0; i < (unsigned int)local_Nx; i++ ) {
			for ( j = 0; j < NAzimuthal; j++ ) {
				l = i*stride + j;
				ih = i+Zero_or_active;
				if ( (i+local_i_start) < GlobalNRadial ) {
					lh = ih*NAzimuthal + j;
					S_radial[l] = dens[lh] * sqrt(Rmed[ih] / GlobalRmed[0]);
					S_azimuthal[l] = S_radial[l] * Rmed[ih] / GlobalRmed[0];
				} else {
					S_radial[l] = 0.;
					S_azimuthal[l] = 0.;
				}
			}
		}
	}

	// upper half of cpus
	if ( (CPU_Rank >= CPU_Number/2) && (CPU_Rank!=CPU_NoFriend) ) {
		for ( i = 0; i < (unsigned int)local_Nx; i++ ) {
			for ( j = 0; j < NAzimuthal; j++ ) {
				l = i*stride + j;
				if ( (i+local_i_start) >= GlobalNRadial ) {
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
	double theta, u;
	int stride = 2*(NAzimuthal/2+1);

	for (unsigned int i = 0; i < (unsigned int)local_Nx; i++) {
		if (i+local_i_start < GlobalNRadial) {
			u = log(Radii[i+local_i_start]/Radii[0]);
		} else {
			u = -log(Radii[2*GlobalNRadial-(i+local_i_start)]/Radii[0]);
		}

		for (unsigned int j = 0; j < NAzimuthal; j++) {
            double denominator;
			int l = i*stride + j;
			theta = 2.0*PI*(double)j/(double)NAzimuthal;

      denominator = pow(2*(cosh(u) - cos(theta)) + lambda_sq*(exp(u)+exp(-u) -2) + chi_sq,-1.5);

			K_radial[l] = 1.0 - cos(theta)*exp(-u);
			K_radial[l] *= denominator;
			K_azimuthal[l] = sin(theta);
			K_azimuthal[l] *= denominator;
		}
	}

	fftw_execute(fftplan_forward_K_radial);
	fftw_execute(fftplan_forward_K_azimuthal);
}

void compute_acceleration(t_polargrid &density)
{
	MPI_Request req1, req3, req4, req5, req6;
	unsigned int i,j, nr;
	int l;
	int one_if_odd;
	int stride;
	int ghost_size;
	double normaccr, normacct;
	nr = density.Nrad;
	stride = 2*(NAzimuthal/2 + 1);
	one_if_odd = (CPU_Number%2 == 0 ? 0 : 1);

	// First we compute sg_acc as a convolution product of reduced density and kernel arrays.
	// In fact, all bufffttabs are transposed arrays, since we used the flag FFTW_TRANSPOSED_ORDER before.
	// However, this is not a problem: we compute double and imaginary parts of 1-dimension convolution arrays buffer_FFT_acc_radial(t)
	for (unsigned int i = 0; i < total_local_size; i++) {
		FFT_acc_radial[i][0] = -constants::G * (FFT_K_radial[i][0]*FFT_S_radial[i][0] - FFT_K_radial[i][1]*FFT_S_radial[i][1]);
		FFT_acc_radial[i][1] = -constants::G * (FFT_K_radial[i][0]*FFT_S_radial[i][1] + FFT_K_radial[i][1]*FFT_S_radial[i][0]);

		FFT_acc_azimuthal[i][0] = -constants::G * (FFT_K_azimuthal[i][0]*FFT_S_azimuthal[i][0] - FFT_K_azimuthal[i][1]*FFT_S_azimuthal[i][1]);
		FFT_acc_azimuthal[i][1] = -constants::G * (FFT_K_azimuthal[i][0]*FFT_S_azimuthal[i][1] + FFT_K_azimuthal[i][1]*FFT_S_azimuthal[i][0]);
	}

	fftw_execute(fftplan_backward_acc_radial);
	fftw_execute(fftplan_backward_acc_azimuthal);

	// The use again of argument FFTW_TRANSPOSED_ORDER in these backward fourier transforms ensures that arrays are not transposed anymore.
	// Then we transfer the exact necessary quantity of sg_acceleration arrays from the fftw d.d. to the hydro mesh
	if ( CPU_Rank != CPU_NoFriend ) {
		if ( CPU_Rank < CPU_Number/2 ) {
			for ( i = 0; i < transfer_size; i++ ) {
				if ( i < transfer_size/2 )
					ffttohydro_transfer[i] = acc_radial[(ifront+1-CPUOVERLAP)*stride + i];
				else
					ffttohydro_transfer[i] = acc_azimuthal[(ifront+1-CPUOVERLAP)*stride + i - transfer_size/2];
			}
			MPI_Isend(ffttohydro_transfer, transfer_size, MPI_DOUBLE, CPU_Friend, 40, MPI_COMM_WORLD, &req1);
            MPI_Wait(&req1, &global_MPI_Status);
		} else {
			MPI_Irecv(ffttohydro_transfer_friend, transfer_size_friend, MPI_DOUBLE, CPU_Friend, 40, MPI_COMM_WORLD, &req1);
            MPI_Wait(&req1, &global_MPI_Status);
		}
	}

	// We now compute sg_acceleration arrays on the hydro mesh
	if ( CPU_Rank < (CPU_Number+one_if_odd)/2 )  {
		if ( CPU_Rank == 0 ) {
			for ( i = 0 ; i < nr; i++ ) {
				for ( j = 0; j < NAzimuthal; j++ ) {
 					l = i*NAzimuthal + j;
					g_radial[l] = acc_radial[i*stride + j];
					g_azimuthal[l] = acc_azimuthal[i*stride + j];
				}
			}
		} else {
			for ( i = Zero_or_active ; i < nr; i++ ) {
				for ( j = 0; j < NAzimuthal; j++ ) {
					l = i*NAzimuthal + j;
					g_radial[l] = acc_radial[(i-Zero_or_active)*stride + j];
					g_azimuthal[l] = acc_azimuthal[(i-Zero_or_active)*stride + j];
				}
			}
		}
	}

	if ( CPU_Rank >= (CPU_Number+one_if_odd)/2 ) {
		if ( CPU_Rank == CPU_Highest ) {
			for ( i = 0 ; i < nr; i++ ) {
				for ( j = 0; j < NAzimuthal; j++ ) {
					l = i*NAzimuthal + j;
					g_radial[l] = ffttohydro_transfer_friend[i*stride + j];
					g_azimuthal[l] = ffttohydro_transfer_friend[transfer_size_friend/2 + i*stride + j];
				}
			}
		} else {
			for ( i = 0 ; i < Max_or_active; i++ ) {
				for ( j = 0; j < NAzimuthal; j++ ) {
					l = i*NAzimuthal + j;
					g_radial[l] = ffttohydro_transfer_friend[i*stride + j];
					g_azimuthal[l] = ffttohydro_transfer_friend[transfer_size_friend/2 + i*stride + j];
				}
			}
		}
	}

	// Now we exchange the correct amount of sg_acceleration between cpus to fill ghosts.

	ghost_size = CPUOVERLAP * NAzimuthal;
	if ( CPU_Number > 1 ) {
		if ( (CPU_Rank > 0) && (CPU_Rank < (CPU_Number+one_if_odd)/2) ) {
			MPI_Isend(&g_azimuthal[Zero_or_active*NAzimuthal], ghost_size, MPI_DOUBLE, CPU_Prev, 60, MPI_COMM_WORLD, &req3);
            MPI_Wait(&req3, &global_MPI_Status);
			MPI_Irecv(&g_azimuthal[0], ghost_size, MPI_DOUBLE, CPU_Prev, 61, MPI_COMM_WORLD, &req4);
            MPI_Wait(&req4, &global_MPI_Status);
		}
		if ( (CPU_Rank >= (CPU_Number+one_if_odd)/2) && (CPU_Rank != CPU_Highest) ) {
			MPI_Irecv(&g_azimuthal[Max_or_active*NAzimuthal], ghost_size, MPI_DOUBLE, CPU_Next, 60, MPI_COMM_WORLD, &req3);
            MPI_Wait(&req3, &global_MPI_Status);
			MPI_Isend(&g_azimuthal[(Max_or_active-CPUOVERLAP)*NAzimuthal], ghost_size, MPI_DOUBLE, CPU_Next, 61, MPI_COMM_WORLD, &req4);
            MPI_Wait(&req4, &global_MPI_Status);
		}
		if ( (CPU_Rank > 0) && (CPU_Rank < (CPU_Number+one_if_odd)/2) ) {
			MPI_Isend(&g_radial[Zero_or_active*NAzimuthal], ghost_size, MPI_DOUBLE, CPU_Prev, 50, MPI_COMM_WORLD, &req5);
            MPI_Wait(&req5, &global_MPI_Status);
			MPI_Irecv(&g_radial[0], ghost_size, MPI_DOUBLE, CPU_Prev, 51, MPI_COMM_WORLD, &req6);
            MPI_Wait(&req6, &global_MPI_Status);
		}
		if ( (CPU_Rank >= (CPU_Number+one_if_odd)/2) && (CPU_Rank != CPU_Highest) ) {
			MPI_Irecv(&g_radial[Max_or_active*NAzimuthal], ghost_size, MPI_DOUBLE, CPU_Next, 50, MPI_COMM_WORLD, &req5);
            MPI_Wait(&req5, &global_MPI_Status);
			MPI_Isend(&g_radial[(Max_or_active-CPUOVERLAP)*NAzimuthal], ghost_size, MPI_DOUBLE, CPU_Next, 51, MPI_COMM_WORLD, &req6);
            MPI_Wait(&req6, &global_MPI_Status);
		}
	}

	// We don't forget to renormalize acc arrays!
	for ( i = 0 ; i < nr; i++ ) {
		// g_r(u,phi) normalized with exp(-u/2)*Δu*Δphi/(2*N_r*N_phi) (3.43 page 57)
		// g_phi(u,phi) normalized with exp(-3*u/2)*Δu Δphi/(2*N_r*N_phi) (3.44 page 57)
		normaccr = r_step * t_step / ( (double)(2*GlobalNRadial) * (double)NAzimuthal );
		normacct = normaccr;
		normaccr /= sqrt(Rmed[i] / GlobalRmed[0]);
		normacct /= ( Rmed[i] / GlobalRmed[0] * sqrt(Rmed[i] / GlobalRmed[0]) );
		for ( j = 0; j < NAzimuthal; j++ ) {
			l = i*NAzimuthal + j;
			g_radial[l] *= normaccr;
			g_azimuthal[l] *= normacct;
		}
	}
}

/**
	Update the velocity fields to take into account self-gravity

	\param VRad
	\param VTheta
	\param Dt
*/
void update_velocities(t_polargrid &v_radial, t_polargrid &v_azimuthal, double dt)
{
	int l;
	int jm1, lm1;

	// Here we update velocity fields to take into account self-gravity
	for (unsigned int n_radial = 0 ; n_radial < v_radial.get_max_radial(); ++n_radial ) {
		for (unsigned int n_azimuthal = 0; n_azimuthal < NAzimuthal; ++n_azimuthal) {
			l = n_radial*NAzimuthal + n_azimuthal;
			// We compute VRadial - half-centered in azimuth - from centered-in-cell radial sg acceleration
			if (n_radial > 0) {
				v_radial(n_radial, n_azimuthal) += dt*( (Rinf[n_radial] - Rmed[n_radial-1]) * g_radial[l] + (Rmed[n_radial] - Rinf[n_radial]) * g_radial[l-NAzimuthal] ) * InvDiffRmed[n_radial];
			}

			// We compute VAzimuthal - half-centered in radius - from centered-in-cell azimutal sg acceleration
			if (n_azimuthal == 0)
				jm1 = NAzimuthal-1;
			else
				jm1 = n_azimuthal-1;
			lm1 = n_radial*NAzimuthal + jm1;
			v_azimuthal(n_radial, n_azimuthal) += 0.5 * dt * (g_azimuthal[l] + g_azimuthal[lm1]);
		}
	}
}

void init_azimuthal_velocity(t_polargrid &v_azimuthal)
{
	double omega;

	double *GLOBAL_AxiSGAccr = (double*)malloc(sizeof(double) * GlobalNRadial);
	mpi_make1Dprofile(g_radial, GLOBAL_AxiSGAccr);

	for (unsigned int n_radial = 0; n_radial <= v_azimuthal.get_max_radial() - GHOSTCELLS_B; ++n_radial) {
		// this corresponds to equation (3.42) in Baruteau, 2008
		double temp = pow2(omega_kepler(Rmed[n_radial])) *( 1.0 - (1.+SIGMASLOPE-2.0*FLARINGINDEX) * pow2(ASPECTRATIO_REF)*pow(Rmed[n_radial],2.0*FLARINGINDEX) ) - GLOBAL_AxiSGAccr[n_radial+IMIN]/Rmed[n_radial];
		if (temp < 0) {
			logging::print("Radicand %lg < 0 in init_azimuthal_velocity! Maybe ThicknessSmoothingSG (%lg) is too small!\n", temp, parameters::thickness_smoothing_sg);
		}
		omega = sqrt(temp);

		for (unsigned int n_azimuthal = 0; n_azimuthal <= v_azimuthal.get_max_azimuthal(); ++n_azimuthal) {
			v_azimuthal(n_radial, n_azimuthal) = Rmed[n_radial]*omega;
		}
	}

	free(GLOBAL_AxiSGAccr);
}

void init_planetary_system(t_data &data)
{
	int ipl;
	double dist, ri, rip1, dr, sgacc;

	double *GLOBAL_AxiSGAccr = (double*)malloc(sizeof(double) * GlobalNRadial);

	mpi_make1Dprofile(g_radial, GLOBAL_AxiSGAccr);

	// Planetary system initialization in self-gravity cases: planets are put in a fixed circular orbit, we need to know radial sg acceleration felt by planets

	for (unsigned int k = 0; k < data.get_planetary_system().get_number_of_planets(); ++k) {
		/* dist denotes the planet's semi-major axis */
		if(data.get_planetary_system().get_planet(k).get_feeldisk())
		{
			dist = data.get_planetary_system().get_planet(k).get_semi_major_axis();
			// TODO remove quick&dirty solution and think about it

			if (dist<RMAX) {
				ipl = 0;
				while (GlobalRmed[ipl] <= dist)
					ipl++;
				ri = GlobalRmed[ipl];
				rip1 = GlobalRmed[ipl+1];
				dr = rip1 - ri;
				sgacc = (dist - ri)*GLOBAL_AxiSGAccr[ipl+1] + (rip1 - dist)*GLOBAL_AxiSGAccr[ipl];
				sgacc /= dr;
				// TODO: fix that only y velocity is updated
				/* sgacc is the radial sg acc. at the planet's semi-major axis */
				double new_vy = data.get_planetary_system().get_planet(k).get_vy() * (double)sqrt (1. - dist*dist*sgacc);
				data.get_planetary_system().get_planet(k).set_vy(new_vy);
			}
		}
	}

	free(GLOBAL_AxiSGAccr);
}

}
