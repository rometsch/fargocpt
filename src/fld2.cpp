#include "fld2.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "parameters.h"
#include "constants.h"
#include "units.h"
#include "opacity.h"
#include "pvte_law.h"
#include "SourceEuler.h"
#include "boundary_conditions.h"
#include "logging.h"
#include "compute.h"


namespace fld2 {

// Radiative transport using Flux-Limited-Diffusion
// This was first implemented by Tobias Müller according to his Ph.D. thesis 2013. 
// See Appendix A therein for the equations of discretization.
// The document can be found at http://nbn-resolving.de/urn:nbn:de:bsz:21-opus-72189
// In this implementation we changed the formulation back to FLD in 3 dimensions
// but with the radiation transport restricted to the midplane.
// Instead of using the surface density and the scale height at various locations,
// we compute the midplane density once and then use the precomputed values.
// This makes the implementation much cleaner.

t_polargrid Ka; // diffusion coefficient on radial boundaries
t_polargrid Kb; // diffusion coefficient on azimuthal boundaries
t_polargrid A, B, C, D, E; // Matrix elements of linear equation
t_polargrid Erad, Trad, Xold, Xtmp; // Intermediate store for old values.

double Erad_max;
double Erad_min;

bool use_temperature;
bool constant_fluxlimiter;

// The full hydro algorithm needs multiple overlapping cells between
// two adjacent cores because many substeps in which information can flow are involved.
// This number is store in the global variable CPUOVERLAP.
// For radiative transport, only one overlap cell is needed
// because in every iteration of the matrix solver information flows only one step.
// Hence we can safe time during communication by reducing the overlap to one inside this module.
// We need to consider that the arrays are still allocated with the global CPUOVERLAP,
// so we define an nstart and nstop for the iteration of the solver iteration loop.
const unsigned int cpuoverlap = 1;
unsigned int nstart;
unsigned int nstop;
unsigned int nskip;

double *SendInnerBoundary;
double *SendOuterBoundary;
double *RecvInnerBoundary;
double *RecvOuterBoundary;

static inline double flux_limiter(const double R)
{
	if (constant_fluxlimiter) {
		return 1.0/3;
	}
    if (R <= 2) {
	return 2.0 / (3 + std::sqrt(9 + 10 * std::pow(R, 2)));
    } else {
	return 10.0 / (10 * R + 9 + std::sqrt(180 * R + 81));
    }
}

void init(const unsigned int Nrad, const unsigned int Naz) { 
	Ka.set_vector(true);
	Ka.set_size(Nrad, Naz);
	Kb.set_scalar(true);
	Kb.set_size(Nrad, Naz);

	A.set_scalar(true);
	A.set_size(Nrad, Naz);
	B.set_scalar(true);
	B.set_size(Nrad, Naz);
	C.set_scalar(true);
	C.set_size(Nrad, Naz);
	D.set_scalar(true);
	D.set_size(Nrad, Naz);
	E.set_scalar(true);
	E.set_size(Nrad, Naz);

	Trad.set_scalar(true);
	Trad.set_size(Nrad, Naz);

	Erad.set_scalar(true);
	Erad.set_size(Nrad, Naz);

	Xold.set_scalar(true);
	Xold.set_size(Nrad, Naz);

	// if (parameters::radiative_diffusion_solver == parameters::t_radiative_diffusion_solver::Jacobi) {
		Xtmp.set_scalar(true);
		Xtmp.set_size(Nrad, Naz);
	// }

    SendInnerBoundary = (double *)malloc(Naz * cpuoverlap * sizeof(double));
    SendOuterBoundary = (double *)malloc(Naz * cpuoverlap * sizeof(double));
    RecvInnerBoundary = (double *)malloc(Naz * cpuoverlap * sizeof(double));
    RecvOuterBoundary = (double *)malloc(Naz * cpuoverlap * sizeof(double));

	constant_fluxlimiter = parameters::radiative_diffusion_test_2d || parameters::radiative_diffusion_test_1d;
	use_temperature = parameters::radiative_diffusion_variable == parameters::t_radiative_diffusion_variable::temperature;

	const double aR = 4*constants::sigma/constants::c;
	Erad_min = aR*std::pow(parameters::minimum_temperature, 4);
	Erad_max = aR*std::pow(parameters::maximum_temperature, 4);

	// define start and end radial index for the solver iteration loop
	if (CPU_Rank == 0) {
		nstart = 1;
	} else {
		nstart = CPUOVERLAP;
	}
	nskip = nstart;

	if (CPU_Rank == CPU_Highest) {
		nstop = NRadial - 1;
	} else {
		nstop = NRadial - CPUOVERLAP;
	}
}

void finalize() {
    free(SendInnerBoundary);
    free(SendOuterBoundary);
    free(RecvInnerBoundary);
    free(RecvOuterBoundary);
}


static void set_minimum_energy(t_data &data) {
    
    auto &Sigma = data[t_data::SIGMA];
    auto &Energy = data[t_data::ENERGY];

    const unsigned int Nphi = Sigma.get_size_azimuthal();

	// We set minimum Temperature here such that H (thru Cs) is also computed with the minimum Temperature
	#pragma omp parallel for
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const unsigned int nr_max = Energy.get_max_radial();

	if (CPU_Rank == 0 &&
	    parameters::boundary_inner == parameters::boundary_condition_open) {
	    const double Tmin = parameters::minimum_temperature;

	    const double mu = pvte::get_mu(data, 1, naz);
	    const double gamma_eff = pvte::get_gamma_eff(data, 1, naz);
	    Sigma(0, naz) = Sigma(1, naz);

	    const double minimum_energy = Tmin * Sigma(1, naz) / mu * constants::R / (gamma_eff - 1.0);

	    Energy(0, naz) = minimum_energy;
	}

	if (CPU_Rank == CPU_Highest &&
	    parameters::boundary_outer == parameters::boundary_condition_open) {
	    const double Tmin = parameters::minimum_temperature; 

	    const double mu = pvte::get_mu(data, nr_max - 1, naz);
	    const double gamma_eff = pvte::get_gamma_eff(data, nr_max - 1, naz);
	    Sigma(nr_max, naz) = Sigma(nr_max - 1, naz);

	    const double minimum_energy = Tmin * Sigma(nr_max - 1, naz) / mu * constants::R / (gamma_eff - 1.0);

	    Energy(nr_max, naz) = minimum_energy;
	}
    }
}

static void radial_boundary_diffusion_coeff() {

    const unsigned int Naz = Ka.get_size_azimuthal();
    #pragma omp parallel for
	for (unsigned int naz = 0; naz < Naz; ++naz) {
		const unsigned int nr_max = Ka.get_max_radial();
		if(CPU_Rank == CPU_Highest && parameters::boundary_outer == parameters::boundary_condition_reflecting){
		Ka(nr_max-1, naz) = 0.0;
		}


		if(CPU_Rank == 0 && parameters::boundary_inner == parameters::boundary_condition_reflecting){
		Ka(1, naz) = 0.0;
		}
	}

	#pragma omp parallel for
	// Similar to Tobi's original implementation
	for (unsigned int naz = 0; naz < Naz; ++naz) {
		const unsigned int nr_max = Ka.get_max_radial();
		if(CPU_Rank == CPU_Highest && !(parameters::boundary_outer == parameters::boundary_condition_reflecting || parameters::boundary_outer == parameters::boundary_condition_open)){
		Ka(nr_max-1, naz) = Ka(nr_max-2, naz);
		}


		if(CPU_Rank == 0 && !(parameters::boundary_inner == parameters::boundary_condition_reflecting || parameters::boundary_inner == parameters::boundary_condition_open)){
		Ka(1, naz) = Ka(2, naz);
		}
	}


    // Similar to Tobi's original implementation
    for (unsigned int naz = 0; naz < Naz; ++naz) {
	const unsigned int nr_max = Ka.get_max_radial();
	if (CPU_Rank == CPU_Highest &&
	    !(parameters::boundary_outer ==
		  parameters::boundary_condition_reflecting ||
	      parameters::boundary_outer ==
		  parameters::boundary_condition_open)) {
	    Ka(nr_max - 1, naz) = Ka(nr_max - 2, naz);
	}

	if (CPU_Rank == 0 && !(parameters::boundary_inner ==
				   parameters::boundary_condition_reflecting ||
			       parameters::boundary_inner ==
				   parameters::boundary_condition_open)) {
	    Ka(1, naz) = Ka(2, naz);
	}
    }
}

static inline double get_opacity(const double temperature, const double density) {
	// Convert to cgs units and call the opacity routine.

	const double temperatureCGS = temperature * units::temperature;
	const double densityCGS = density * units::density;

	const double kappaCGS = opacity::opacity(densityCGS, temperatureCGS);
	const double kappa = parameters::kappa_factor * kappaCGS * units::opacity.get_inverse_cgs_factor();

	return kappa;
}


static inline double diffusion_coeff_energy(const double lrad, const double E, const double nabla_E) {
	// FLD diffusion coefficient for diffusion equation of energy

	// Levermore & Pomraning 1981
	// R = |nabla E|/E * lrad
	// lrad = 1/(rho * kappa), photon mean free path

	const double R = nabla_E / E * lrad;
	const double lambda = flux_limiter(R);

	const double c = constants::c;
	const double K = lambda * c * lrad;
	return K;
}

static inline double diffusion_coeff_temperature(const double lrad, const double T, const double nabla_T) {
	// FLD diffusion coefficient for diffusion equation of temperature

	// Levermore & Pomraning 1981
	// R = |nabla E|/E * lrad = 4 |nabla T|/T * lrad
	// lrad = 1/(rho * kappa), photon mean free path
	const double R = 4.0 * nabla_T / T * lrad;
	const double lambda = flux_limiter(R);

	const double sigRad = constants::sigma;
	const double K = lambda * 16 * sigRad * lrad * std::pow(T, 3);
	return K;
}

static inline double Erad2Trad(const double E) {
	const double aR = 4*constants::sigma/constants::c;
	return std::pow(E/aR, 0.25);
}

static inline double diffusion_coeff(const double rho, const double x, const double nabla_x) {
	if (use_temperature) {
		const double T = x;
		const double kappa = get_opacity(T, rho);
		const double lrad = 1 / (rho * kappa);
		return diffusion_coeff_temperature(lrad, x, nabla_x);
	} else {
		const double T = Erad2Trad(x);
		const double kappa = get_opacity(T, rho);
		const double lrad = 1 / (rho * kappa);
		return diffusion_coeff_energy(lrad, x, nabla_x);
	}
}

static void compute_diffusion_coeff_radial(t_data &data, t_polargrid &X) {
	// Calculate the diffusion coefficient on at the radial interface locations.

    auto &Density = data[t_data::RHO];

    const unsigned int Nrad = Ka.get_size_radial();
    const unsigned int Naz = Ka.get_size_azimuthal();
    const unsigned int Naz_last = Ka.get_max_azimuthal();

    // calcuate Ka for K(i/2,j)
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nrad - 1; ++nr) {
	for (unsigned int naz = 0; naz < Naz; ++naz) {
	    const unsigned int naz_next = (naz == Naz_last ? 0 : naz + 1);
	    const unsigned int naz_last = (naz == 0 ? Naz_last : naz - 1);

	    // average radially
	    const double x = 0.5 * (X(nr - 1, naz) + X(nr, naz));
	    const double rho = 0.5 * (Density(nr - 1, naz) + Density(nr, naz));

		// calculate gradient length
	    const double dX_dr = (X(nr, naz) - X(nr - 1, naz)) * InvDiffRmed[nr];
        const double xnext = 0.5 * (X(nr - 1, naz_next) + X(nr, naz_next));
        const double xlast = 0.5 * (X(nr - 1, naz_last) + X(nr, naz_last));
	    const double dX_dphi = InvRinf[nr] * (xnext-xlast) / (2.0 * dphi);

	    const double nabla_X = std::hypot(dX_dr, dX_dphi);

		Ka(nr, naz) = diffusion_coeff(rho, x, nabla_X);
	}
    }
	
}


static void compute_diffusion_coeff_azimuthal(t_data &data, t_polargrid &X) {
	// Calculate the diffusion coefficient on at the azimuthal interface locations.

    // calcuate Kb for K(i,j/2)
    auto &Density = data[t_data::RHO];

	const unsigned int Nrad = Kb.get_size_radial();
    const unsigned int Naz = Kb.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nrad-1; ++nr) {
	for (unsigned int naz = 0; naz < Naz; ++naz) {
		const unsigned int naz_prev = (naz == 0 ? Naz-1 : naz - 1);

	    // average azimuthally
	    const double x = 0.5 * (X(nr, naz_prev) + X(nr, naz));
		const double rho = 0.5 * (Density(nr, naz_prev) + Density(nr, naz));

		// calculate gradient length
		// values are located on the azimuthal interface
		const double Router = Ra[nr + 1];
		const double Rinner = Ra[nr - 1];
		const double xouter = 0.5 * (X(nr + 1, naz_prev) + X(nr + 1, naz));
		const double xinner = 0.5 * (X(nr - 1, naz_prev) + X(nr - 1, naz));
	    const double dX_dr = (xouter - xinner) / (Router - Rinner);
        const double dX_dphi = InvRmed[nr] * (X(nr, naz) - X(nr, naz_prev)) / dphi;

	    const double nabla_X = std::hypot(dX_dr, dX_dphi);

	    Kb(nr, naz) = diffusion_coeff(rho, x, nabla_X);
	}
    }
}

static void save_old_values(t_polargrid &X) {
	const unsigned int Nrad = X.get_size_radial();
    const unsigned int Naz = X.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nrad; ++nr) {
		for (unsigned int naz = 0; naz < Naz; ++naz) {
			Xold(nr, naz) = X(nr, naz);
		}
	}
}

static void calculate_matrix_elements(t_data &data, const double dt) {
	// calculate the matrix elements of the linear equation 

	auto &Density = data[t_data::RHO];

	// TODO: check whether we need gamma_eff here
    const double c_v = constants::R / (parameters::MU * (parameters::ADIABATICINDEX - 1.0));
	// printf("c_v = %e\n", c_v);
	// const double c_v = 1;
	const unsigned int Nrad = A.get_size_radial();
    const unsigned int Naz = A.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nrad - 1; ++nr) {
	for (unsigned int naz = 0; naz < Naz; ++naz) {
		double common_factor;
		if (use_temperature) {
			const double rho = Density(nr, naz);
			// this refers to the common factor in Müller 2013 Ph.D. thesis (A.1.10) adjusted for 3D FLD in midplane
			common_factor = -dt / (rho * c_v);
		} else {
			common_factor = -dt;
		}
		if (nr == 20 && naz == 20) {
			printf("common_factor = %e\n", common_factor);
		}

	    // 2/(dR^2)
	    const double common_AC = common_factor * 2.0 / (std::pow(Ra[nr + 1], 2) - std::pow(Ra[nr], 2));
	    A(nr, naz) = common_AC * Ka(nr, naz) * Ra[nr] * InvDiffRmed[nr];
	    C(nr, naz) = common_AC * Ka(nr + 1, naz) * Ra[nr + 1] * InvDiffRmed[nr + 1];

	    // 1/(r^2 dphi^2)
		const unsigned int naz_next = naz == Kb.get_max_azimuthal() ? 0 : naz + 1;
	    const double common_DE = common_factor / (std::pow(Rb[nr], 2) * std::pow(dphi, 2));
	    D(nr, naz) = common_DE * Kb(nr, naz);
	    E(nr, naz) = common_DE * Kb(nr, naz_next);

	    B(nr, naz) = -A(nr, naz) - C(nr, naz) - D(nr, naz) - E(nr, naz) + 1.0;
	}
    }

}

// Communicate the temperature values among adjacent nodes
static void communicate_parallelization_boundaries(t_polargrid &Temperature) {

	const unsigned int Nrad = Temperature.Nrad;
	const unsigned int Naz = NAzimuthal;

    const int l = cpuoverlap * NAzimuthal;
    const int oo = (Nrad - cpuoverlap) * NAzimuthal;
    const int o = (Nrad - 2 * cpuoverlap) * NAzimuthal;

	const unsigned int start_active = CPUOVERLAP*Naz;
	const unsigned int end_active = (Nrad - CPUOVERLAP)*Naz;
	const unsigned int last_active = end_active - l;

	const unsigned int start_inner_boundary = (CPUOVERLAP-cpuoverlap)*Naz;
	const unsigned int start_outer_boundary = end_active;

	// communicate with other nodes
	memcpy(SendInnerBoundary, Temperature.Field + start_active, l * sizeof(double));
	memcpy(SendOuterBoundary, Temperature.Field + last_active, l * sizeof(double));

	MPI_Request req1, req2, req3, req4;

	if (CPU_Rank % 2 == 0) {
	    if (CPU_Rank != 0) {
		MPI_Isend(SendInnerBoundary, l,
			  MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
		MPI_Irecv(RecvInnerBoundary, l,
			  MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req2);
	    }
	    if (CPU_Rank != CPU_Highest) {
		MPI_Isend(SendOuterBoundary, l,
			  MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req3);
		MPI_Irecv(RecvOuterBoundary, l,
			  MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req4);
	    }
	} else {
	    if (CPU_Rank != CPU_Highest) {
		MPI_Irecv(RecvOuterBoundary, l,
			  MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req3);
		MPI_Isend(SendOuterBoundary, l,
			  MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req4);
	    }
	    if (CPU_Rank != 0) {
		MPI_Irecv(RecvInnerBoundary, l,
			  MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
		MPI_Isend(SendInnerBoundary, l,
			  MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req2);
	    }
	}

	// copy values revceived from inner core to the inner local boundary
	if (CPU_Rank != 0) {
	    MPI_Wait(&req1, &global_MPI_Status);
	    MPI_Wait(&req2, &global_MPI_Status);
	    memcpy(Temperature.Field + start_inner_boundary, RecvInnerBoundary, l * sizeof(double));
	}

	// copy values received from outer core to the outer local boundary
	if (CPU_Rank != CPU_Highest) {
	    MPI_Wait(&req3, &global_MPI_Status);
	    MPI_Wait(&req4, &global_MPI_Status);
	    memcpy(Temperature.Field + start_outer_boundary, RecvOuterBoundary,
		   l * sizeof(double));
	}
}

static void boundary_T_SOR(t_polargrid &X) {
	double minval;
	if (use_temperature) {
		minval = parameters::minimum_temperature;
	} else {
		minval = Erad_min;
	}


    const unsigned int Naz = X.get_size_azimuthal();
    const unsigned int nrad_last = X.get_max_radial();

    const bool at_outer_boundary = CPU_Rank == CPU_Highest;
    const bool outer_boundary_open = parameters::boundary_outer == parameters::boundary_condition_open;
	if (at_outer_boundary && outer_boundary_open) {
		// set temperature to T_min in outermost cells
		for (unsigned int naz = 0; naz < Naz; ++naz) {
			X(nrad_last, naz) = minval;
		}
	}

    const bool at_inner_boundary = CPU_Rank == 0;
    const bool inner_boundary_open = parameters::boundary_inner == parameters::boundary_condition_open;
	if (at_inner_boundary && inner_boundary_open) {
        // set temperature to T_min in innermost cells
		for (unsigned int naz = 0; naz < Naz; ++naz) {
			X(0, naz) = minval;
		}
	}
	// boundary_conditions::apply_boundary_condition(data, current_time, 0.0, false);
}

static inline bool is_active(const unsigned int nr) {
	// Check whether a radial cell index belongs to an active cell.
	// An active cell is neither a ghost cell nor an overlap cell.
    const unsigned int nrad_last = NRadial - 1;
    const unsigned int nghost_right = (CPU_Rank == CPU_Highest) ? GHOSTCELLS_B : CPUOVERLAP;
    const unsigned int nghost_left = (CPU_Rank == 0) ? GHOSTCELLS_B : CPUOVERLAP;
    
    const bool isnot_ghostcell_rank_0 = nr > nghost_left;
    const bool isnot_ghostcell_rank_highest = nr < (nrad_last - nghost_right);

    return isnot_ghostcell_rank_0 && isnot_ghostcell_rank_highest;
}

static inline void floorceil_T_single(double &T) {
	// Set temperature floor and ceiling
	if(T < parameters::minimum_temperature){
		T = parameters::minimum_temperature;
	}

	if(T > parameters::maximum_temperature){
		logging::print(LOG_INFO "max temp inside FLD SOR loop %e\n", T);
		T = parameters::maximum_temperature;
	}
}

static inline void floorceil_E_single(double &E) {
	// Set temperature floor and ceiling
	if(E < Erad_min){
		E = Erad_min;
	}

	if(E > Erad_max){
		logging::print(LOG_INFO "max temp inside FLD SOR loop E = %e, T = %e\n", E, Erad2Trad(E));
		E = Erad_max;
	}
}

static inline void floorceil_x_single(double &x) {
	if (use_temperature) {
		floorceil_T_single(x);
	} else {
		floorceil_E_single(x);
	}
}


static void Jacobi(t_polargrid &X) {
	// Solve the linear system using the Jacobi method
	// Based on the book "Iterative Methods for Sparse Linear Systems" by Yousef Saad
	// DOI: 10.1137/1.9780898718003
	// downloaded 2023-06-28 from the author's academic webpage: 
	// https://www-users.cse.umn.edu/~saad/IterMethBook_2ndEd.pdf
	// The method is described in Chapter 4.1 and the relevant equation is Eq. (4.4)

    const unsigned int Nrad = X.get_size_radial();
    const unsigned int Naz = X.get_size_azimuthal();

    unsigned int iterations = 0;
	double absolute_norm = std::numeric_limits<double>::max();
	double norm_change = std::numeric_limits<double>::max();

    // do SOR
	const double tolerance = parameters::radiative_diffusion_tolerance;
	const unsigned int maxiter = parameters::radiative_diffusion_max_iterations;

	// initialize solution vector with zero
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = nstart; nr < nstop; ++nr) {
		for (unsigned int naz = 0; naz < Naz; ++naz) {
			X(nr, naz) = 0.0;
		}
	}

    while ((norm_change > tolerance) && (maxiter > iterations)) {

    boundary_T_SOR(X);

	copy_polargrid(Xtmp, X);

	norm_change = absolute_norm;
	absolute_norm = 0.0;
	#pragma omp parallel for collapse(2) reduction(+ : absolute_norm)
	for (unsigned int nr = nstart; nr < nstop; ++nr) {
		for (unsigned int naz = 0; naz < Naz; ++naz) {

        const unsigned int naz_last = X.get_max_azimuthal();
		const unsigned int naz_next = (naz == naz_last ? 0 : naz + 1);
		const unsigned int naz_prev = (naz == 0 ? naz_last : naz - 1);
		const double old_value = X(nr, naz);

		const double aii = B(nr, naz);
		const double xn = Xold(nr, naz); // value at iteration n
		const double sum =
			A(nr, naz) * X(nr - 1, naz) +
			C(nr, naz) * X(nr + 1, naz) +
			D(nr, naz) * X(nr, naz_prev) +
			E(nr, naz) * X(nr, naz_next);

		// perform the update step given in Eq. (4.4) in the reference given above
		const double xup = (xn - sum)/aii; // update value at iteration n+1

		Xtmp(nr, naz) = xup;
		// ensure minimum and maximum temperature
		floorceil_x_single(Xtmp(nr, naz));

		// only non ghostcells to norm and don't count overlap cell's
		// twice
		if (is_active(nr)) {
			absolute_norm += std::pow(old_value - Xtmp(nr, naz), 2);
		}
		} // azimuthal loop
	} // radial loop

	copy_polargrid(X, Xtmp);

    const double tmp = absolute_norm;
	MPI_Allreduce(&tmp, &absolute_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	absolute_norm = std::sqrt(absolute_norm) / (GlobalNRadial * NAzimuthal);

	norm_change = fabs(absolute_norm - norm_change);
	// printf("iter %u, norm_change = %e\n", iterations, norm_change);
	iterations++;

    communicate_parallelization_boundaries(X);

    } // END Loop


    if (iterations == parameters::radiative_diffusion_max_iterations) {
	logging::print_master(
	    LOG_WARNING
	    "Maximum iterations (%u) reached in radiative_diffusion with Jacobi solver. Norm is %lg with a last change of %lg.\n",
	    parameters::radiative_diffusion_max_iterations, absolute_norm, norm_change);
    }

    logging::print_master(LOG_VERBOSE "%u iterations\n", iterations);
}


static void SOR(t_polargrid &X) {

    const unsigned int Nrad = X.get_size_radial();
    const unsigned int Naz = X.get_size_azimuthal();

    static unsigned int old_iterations = parameters::radiative_diffusion_max_iterations;
    static int direction = 1;
    static double omega = parameters::radiative_diffusion_omega;

    unsigned int iterations = 0;
	double absolute_norm = std::numeric_limits<double>::max();
	double norm_change = std::numeric_limits<double>::max();

    // do SOR
	const double tolerance = parameters::radiative_diffusion_tolerance;
	const unsigned int maxiter = parameters::radiative_diffusion_max_iterations;

    while ((norm_change > tolerance) && (maxiter > iterations)) {

    boundary_T_SOR(X);

	norm_change = absolute_norm;
	absolute_norm = 0.0;

	// static const unsigned int chunk_size = std::ceil((float)(nstop-nstart+1)/Thread_Number);
	// if (iterations==0) {
	// 	printf("[%i] Thread number = %i, chunk size = %u, size = %i\n", CPU_Rank, Thread_Number, chunk_size, nstop - nstart+1);
	// }
// #pragma omp parallel for collapse(2) shared(X) reduction(+ : absolute_norm) schedule(static,1)
// #pragma omp parallel for reduction(+ : absolute_norm) schedule(dynamic)
// #pragma omp parallel for reduction(+ : absolute_norm) schedule(dynamic, chunk_size)
	// int id = omp_get_thread_num();
    // int b = id * chunk_size;
    // int e = id == n_threads - 1 ? n : b + chunk_size;
    // printf("thread %d: %d items\n", id, e - b);
    // for (int i = b; i < e; i++) {
    //   // process item i
    // }

#pragma omp parallel for collapse(2) shared(X) reduction(+ : absolute_norm) schedule(static,1)
	for (unsigned int nr = nstart; nr < nstop; ++nr) {
		for (unsigned int naz = 0; naz < Naz; ++naz) {

		const double old_value = X(nr, naz);
        const unsigned int naz_last = X.get_max_azimuthal();
		const unsigned int naz_next = (naz == naz_last ? 0 : naz + 1);
		const unsigned int naz_prev = (naz == 0 ? naz_last : naz - 1);

		X(nr, naz) =
		    (1.0 - omega) * X(nr, naz) -
		    omega / B(nr, naz) *
			(A(nr, naz) * X(nr - 1, naz) +
			 C(nr, naz) * X(nr + 1, naz) +
			 D(nr, naz) * X(nr, naz_prev) +
			 E(nr, naz) * X(nr, naz_next) - Xold(nr, naz));

		// ensure minimum and maximum temperature
		floorceil_x_single(X(nr, naz));

		// only non ghostcells to norm and don't count overlap cell's
		// twice
        if (is_active(nr)) {
		    absolute_norm += std::pow(old_value - X(nr, naz), 2);
		}
	    }
	}

    const double tmp = absolute_norm;
	MPI_Allreduce(&tmp, &absolute_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	absolute_norm = std::sqrt(absolute_norm) / (GlobalNRadial * NAzimuthal);

	norm_change = fabs(absolute_norm - norm_change);
	// printf("iter %u, norm_change = %e\n", iterations, norm_change);
	iterations++;

    communicate_parallelization_boundaries(X);

    } // END SOR

    if (iterations == parameters::radiative_diffusion_max_iterations) {
	logging::print_master(
	    LOG_WARNING
	    "Maximum iterations (%u) reached in radiative_diffusion (omega = %lg). Norm is %lg with a last change of %lg.\n",
	    parameters::radiative_diffusion_max_iterations, omega,
	    absolute_norm, norm_change);
    }

    // adapt omega
    if (old_iterations < iterations) {
	direction *= -1;
    }

    if (parameters::radiative_diffusion_omega_auto_enabled) {
	omega += direction * 0.01;
    }

    if (omega >= 2.0) {
	omega = 1.99;
	direction = -1;
    }

    if (omega <= 1.0) {
	omega = 1.0;
	direction = 1;
    }

    old_iterations = iterations;

    logging::print_master(LOG_VERBOSE "%u iterations, omega=%lf\n", iterations,
			  omega);
}

static void update_energy(t_data &data, t_polargrid &T) {

    auto &Sigma = data[t_data::SIGMA];
    auto &Energy = data[t_data::ENERGY];

    const unsigned int Nrad = T.get_size_radial();
    const unsigned int Naz = T.get_size_azimuthal();

    const double c_v = constants::R / (parameters::ADIABATICINDEX - 1.0) / parameters::MU;
    // compute energy from temperature
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nrad - 1; ++nr) {
        for (unsigned int naz = 0; naz < Naz; ++naz) {
            Energy(nr, naz) = c_v * T(nr, naz) * Sigma(nr, naz);
        }
    }
}

static void set_constant_midplane_density(t_data &data){
	// set the midplane density to a constant 1 g/cm3 for testing purposes
	auto &Rho = data[t_data::RHO];
	
	const unsigned int Nrad = Rho.get_size_radial();
	const unsigned int Naz = Rho.get_size_azimuthal();

	const double density = parameters::radiative_diffusion_test_2d_density;

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nrad; ++nr) {
        for (unsigned int naz = 0; naz < Naz; ++naz) {
            Rho(nr, naz) = density;
        }
    }

}


static void print_temperature_over_one(t_data &data){
	auto &T = data[t_data::TEMPERATURE];
	
	const unsigned int Nrad = T.get_size_radial();
	const unsigned int Naz = T.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nrad; ++nr) {
        for (unsigned int naz = 0; naz < Naz; ++naz) {
			if (T(nr, naz) > 1) {
				printf("!!!!!!!!!!!!!!!!! T(%u,%u) = %e\n", nr, naz, T(nr, naz));
				const double aR = 4*constants::sigma/constants::c;
				printf("!!!!!!!!!!!!!!!!! aR T^4 = %e\n", aR*std::pow(T(nr, naz), 4));
			}
        }
    }

}

static void compute_radiation_energy_density(t_polargrid &E, t_polargrid &T){
	
	const unsigned int Nrad = E.get_size_radial();
	const unsigned int Naz = E.get_size_azimuthal();

	const double aR = 4*constants::sigma/constants::c;

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nrad; ++nr) {
        for (unsigned int naz = 0; naz < Naz; ++naz) {
			E(nr, naz) =  aR*std::pow(T(nr, naz), 4);
        }
    }
}

static void compute_radiation_temperature(t_polargrid &T, t_polargrid &E){
	
	const unsigned int Nrad = E.get_size_radial();
	const unsigned int Naz = E.get_size_azimuthal();

	const double aR = 4*constants::sigma/constants::c;

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nrad; ++nr) {
        for (unsigned int naz = 0; naz < Naz; ++naz) {
			T(nr, naz) =  std::pow(E(nr, naz)/aR, 0.25);
        }
    }
}

static void floor_T(t_polargrid &T){
	
	const unsigned int Nrad = T.get_size_radial();
	const unsigned int Naz = T.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nrad; ++nr) {
        for (unsigned int naz = 0; naz < Naz; ++naz) {
			if (T(nr, naz) < parameters::minimum_temperature) {
				T(nr, naz) = parameters::minimum_temperature;
			}
        }
    }
}


void check_solution(t_polargrid &X) {
	
	double diff = 0;
	double diff_abs = 0;

	const unsigned int Nrad = X.get_size_radial();
	const unsigned int Naz = X.get_size_azimuthal();

	#pragma omp parallel for collapse(2) reduction(+ : diff, diff_abs)
	for (unsigned int nr = nstart; nr < nstop; ++nr) {
		for (unsigned int naz = 0; naz < Naz; ++naz) {

        const unsigned int naz_last = X.get_max_azimuthal();
		const unsigned int naz_next = (naz == naz_last ? 0 : naz + 1);
		const unsigned int naz_prev = (naz == 0 ? naz_last : naz - 1);

		const double xold = Xold(nr, naz); // value at iteration n
		const double xnew =
			A(nr, naz) * X(nr - 1, naz) +
			B(nr, naz) * X(nr, naz) +
			C(nr, naz) * X(nr + 1, naz) +
			D(nr, naz) * X(nr, naz_prev) +
			E(nr, naz) * X(nr, naz_next);

		// only non ghostcells to norm and don't count overlap cell's
		// twice
		if (is_active(nr)) {
			diff += xnew - xold;
			diff_abs += std::abs(xnew - xold);
		}
		} // azimuthal loop
	} // radial loop


    double tmp = diff;
	MPI_Allreduce(&tmp, &diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	printf("Sum of differences in RD solver : %e\n", diff);

	tmp = diff_abs;
	MPI_Allreduce(&tmp, &diff_abs, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	printf("Sum of |differences| in RD solver : %e\n", diff_abs);
}


void check_matrix() {
	
	double diff = 0;
	double diff_abs = 0;

	const unsigned int Nrad = A.get_size_radial();
	const unsigned int Naz = A.get_size_azimuthal();

	#pragma omp parallel for collapse(2) reduction(+ : diff, diff_abs)
	for (unsigned int nr = nstart; nr < nstop; ++nr) {
		for (unsigned int naz = 0; naz < Naz; ++naz) {

		if (!is_active(nr)) {
			continue;
		}

		const double d = B(nr, naz);
		const double a1 = A(nr, naz);
		const double a2 = C(nr, naz);
		const double a3 = D(nr, naz);
		const double a4 = E(nr, naz);

		if (std::abs(a1) > std::abs(d)) {
			printf("Off diagonal element A = %e > diagonal %e at (%u, %u)\n", a1, d, nr, naz);
		}

		if (std::abs(a2) > std::abs(d)) {
			printf("Off diagonal element C = %e > diagonal %e at (%u, %u)\n", a2, d, nr, naz);
		}

		if (std::abs(a3) > std::abs(d)) {
			printf("Off diagonal element D = %e > diagonal %e at (%u, %u)\n", a3, d, nr, naz);
		}

		if (std::abs(a4) > std::abs(d)) {
			printf("Off diagonal element E = %e > diagonal %e at (%u, %u)\n", a4, d, nr, naz);
		}


		} // azimuthal loop
	} // radial loop


    // double tmp = diff;
	// MPI_Allreduce(&tmp, &diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	// printf("Sum of differences in RD solver : %e\n", diff);

}


static void radiative_diffusion_energy(t_data &data, const double dt) {

	// instantaneously equilibrate photon gas temperature to ideal gas temperature
	compute_radiation_energy_density(Erad, data[t_data::TEMPERATURE]);

	// calculate diffusion coefficient
	compute_diffusion_coeff_radial(data, Erad);
    radial_boundary_diffusion_coeff();
    compute_diffusion_coeff_azimuthal(data, Erad);

	// setup linear equation and solve the diffusion equation
    calculate_matrix_elements(data, dt);
	save_old_values(Erad);

	if (parameters::radiative_diffusion_solver == parameters::t_radiative_diffusion_solver::SOR) {
		SOR(Erad);
	} else {
		Jacobi(Erad);
	}

	check_matrix();
	check_solution(Erad);

	// instantaneously equilibrate ideal gas temperature to photon gas temperature
	compute_radiation_temperature(Trad, Erad);
	copy_polargrid(data[t_data::TEMPERATURE], Trad);
}

static void radiative_diffusion_temperature(t_data &data, const double dt) {

	// instantaneously equilibrate photon gas temperature to ideal gas temperature
	copy_polargrid(Trad, data[t_data::TEMPERATURE]);

	// calculate diffusion coefficient
	compute_diffusion_coeff_radial(data, Trad);
    radial_boundary_diffusion_coeff();
    compute_diffusion_coeff_azimuthal(data, Trad);

	// setup linear equation and solve the diffusion equation
    calculate_matrix_elements(data, dt);
	save_old_values(Trad);

	if (parameters::radiative_diffusion_solver == parameters::t_radiative_diffusion_solver::SOR) {
		SOR(Trad);
	} else {
		Jacobi(Trad);
	}

	check_matrix();
	check_solution(Trad);

	// instantaneously equilibrate ideal gas temperature to photon gas temperature
	copy_polargrid(data[t_data::TEMPERATURE], Trad);

	compute_radiation_energy_density(Erad, Trad);
}


void radiative_diffusion(t_data &data, const double current_time, const double dt)
{
    set_minimum_energy(data);
	
    // update temperature, soundspeed and aspect ratio
	compute_temperature(data);
	// floor_T(data[t_data::TEMPERATURE]);
	compute_sound_speed(data, current_time);
	// compute_scale_height(data, current_time); done in compute::midplane_density
	compute::midplane_density(data, current_time);

	if (parameters::radiative_diffusion_test_2d) {
		set_constant_midplane_density(data);
	}

	if (use_temperature) {
		radiative_diffusion_temperature(data, dt);
	} else {
		radiative_diffusion_energy(data, dt);
	}

	update_energy(data, Trad);

    // ensure energy is consistent with min/max temperature
	SetTemperatureFloorCeilValues(data, __FILE__, __LINE__);
}


}