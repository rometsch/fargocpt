#include "fld.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "parameters.h"
#include "constants.h"
#include "config.h"
#include "units.h"
#include "opacity.h"
#include "pvte_law.h"
#include "SourceEuler.h"
#include "boundary_conditions/boundary_conditions.h"
#include "logging.h"
#include "compute.h"
#include "output.h"


namespace fld {

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
t_polargrid Trad, Told; // Intermediate store for old values.

bool constant_fluxlimiter;

// set tolrance value relative to temperature floor
double tolerance;

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


// Parameters
/// enable radiative diffusion
bool radiative_diffusion_enabled;
/// omega for SOR in radiative diffusion
double radiative_diffusion_omega;
/// enable automatic omega in SOR in radiative diffusion
bool radiative_diffusion_omega_auto_enabled;
/// maximum iterations in SOR in radiative diffusion
unsigned int radiative_diffusion_max_iterations;
/// tolerance for the radiative diffusion
double radiative_diffusion_tolerance;
// enable 2d test
bool radiative_diffusion_test_2d;
// constant density for 2d fld test
double radiative_diffusion_test_2d_density;
// constant K for 2d fld test
double radiative_diffusion_test_2d_K;
// number of diffusion steps to perform for the 2d test
unsigned int radiative_diffusion_test_2d_steps;
// enable 1d test
bool radiative_diffusion_test_1d;
// save internal grids of the fld module at each snapshot
bool radiative_diffusion_dump_data;
// check solution of linear system
bool radiative_diffusion_check_solution;


void config() {

	radiative_diffusion_enabled = config::cfg.get_flag("RadiativeDiffusion", "No");

    radiative_diffusion_omega = config::cfg.get<double>("RadiativeDiffusionOmega", 1.5);
    radiative_diffusion_omega_auto_enabled = config::cfg.get_flag("RadiativeDiffusionAutoOmega", "no");
    radiative_diffusion_max_iterations = config::cfg.get<unsigned int>("RadiativeDiffusionMaxIterations", 50000);
    radiative_diffusion_tolerance = config::cfg.get<double>("RadiativeDiffusionTolerance", 1.5);

	// test parameters
    radiative_diffusion_test_2d_K = config::cfg.get<double>("RadiativeDiffusionTest2DK", 1.0);
    radiative_diffusion_test_2d_steps = config::cfg.get<unsigned int>("RadiativeDiffusionTest2DSteps", 1);
	units::precise_unit L0 = units::L0;
    units::precise_unit M0 = units::M0;
    radiative_diffusion_test_2d_density = config::cfg.get<double>("RadiativeDiffusionTest2DDensity", "1.0 g/cm3", M0/(L0*L0*L0));
    radiative_diffusion_test_2d = config::cfg.get_flag("RadiativeDiffusionTest2D", "no");
    radiative_diffusion_test_1d = config::cfg.get_flag("RadiativeDiffusionTest1D", "no");
    radiative_diffusion_dump_data = config::cfg.get_flag("RadiativeDiffusionDumpData", "no");
    radiative_diffusion_check_solution = config::cfg.get_flag("RadiativeDiffusionCheckSolution", "no");


	logging::print_master(LOG_INFO
	"Radiative diffusion is %s. Using %s omega = %lf with a maximum %u interations.\n",
	radiative_diffusion_enabled ? "enabled" : "disabled",
	radiative_diffusion_omega_auto_enabled ? "auto" : "fixed",
	radiative_diffusion_omega, radiative_diffusion_max_iterations);
}


/*
Flux limiter after Kley (1989) https://ui.adsabs.harvard.edu/abs/1989A&A...208...98K
*/
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

/*
Setup variables and arrays for the fld module.
*/
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

	Told.set_scalar(true);
	Told.set_size(Nrad, Naz);

	Ka.set_name("Ka");
	Kb.set_name("Kb");
	A.set_name("A");
	B.set_name("B");
	C.set_name("C");
	D.set_name("D");
	E.set_name("E");
	Trad.set_name("Trad");
	Told.set_name("Told");

    SendInnerBoundary = (double *)malloc(Naz * cpuoverlap * sizeof(double));
    SendOuterBoundary = (double *)malloc(Naz * cpuoverlap * sizeof(double));
    RecvInnerBoundary = (double *)malloc(Naz * cpuoverlap * sizeof(double));
    RecvOuterBoundary = (double *)malloc(Naz * cpuoverlap * sizeof(double));

	constant_fluxlimiter = radiative_diffusion_test_2d || radiative_diffusion_test_1d;

	const double reltol = radiative_diffusion_tolerance;
	tolerance = reltol*parameters::minimum_temperature;
	logging::print_master(LOG_INFO "FLD solver has tolerance %e in code units.\n", tolerance);

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


/*
If requested, write out all 2d arrays to the current snapshot directory.
*/
void handle_output() {
	if (!radiative_diffusion_dump_data) {
		return;
	}
	Ka.write2D();
	Kb.write2D();
	A.write2D();
	B.write2D();
	C.write2D();
	D.write2D();
	E.write2D();
	Trad.write2D();
	Told.write2D();
}

void finalize() {
    free(SendInnerBoundary);
    free(SendOuterBoundary);
    free(RecvInnerBoundary);
    free(RecvOuterBoundary);
}


/*
Set floor value for internal energy.
These values are different from the floor values of the radiative energy!
*/
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

/*
Boundary conditions for the diffusion coefficient.
*/
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

/*
Get opacity in code units.

First convert temperature and density to cgs and then convert opacity from cgs to code units.
Convert to cgs units and call the opacity routine.
*/
static inline double get_opacity(const double temperature, const double density) {

	const double temperatureCGS = temperature * units::temperature;
	const double densityCGS = density * units::density;

	const double kappaCGS = opacity::opacity(densityCGS, temperatureCGS);
	const double kappa = parameters::kappa_factor * kappaCGS * units::opacity.get_inverse_cgs_factor();

	return kappa;
}

/*
FLD diffusion coefficient for diffusion equation of temperature
*/
static inline double diffusion_coeff(const double rho, const double T, const double nabla_T) {

	// Levermore & Pomraning 1981
	// R = |nabla E|/E * lrad = 4 |nabla T|/T * lrad
	// lrad = 1/(rho * kappa), photon mean free path
	const double kappa = get_opacity(T, rho);
	const double lrad = 1 / (rho * kappa);

	const double R = 4.0 * nabla_T / T * lrad;
	const double lambda = flux_limiter(R);

	const double sigRad = constants::sigma;
	const double K = lambda * 16 * sigRad * lrad * std::pow(T, 3);
	return K;
}


/*
Calculate the diffusion coefficient on at the radial interface locations.
*/
static void compute_diffusion_coeff_radial(t_polargrid &Density, t_polargrid &T) {

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
	    const double x = 0.5 * (T(nr - 1, naz) + T(nr, naz));
	    const double rho = 0.5 * (Density(nr - 1, naz) + Density(nr, naz));

		// calculate gradient length
	    const double dX_dr = (T(nr, naz) - T(nr - 1, naz)) * InvDiffRmed[nr];
        const double xnext = 0.5 * (T(nr - 1, naz_next) + T(nr, naz_next));
        const double xlast = 0.5 * (T(nr - 1, naz_last) + T(nr, naz_last));
	    const double dX_dphi = InvRinf[nr] * (xnext-xlast) / (2.0 * dphi);

	    const double nabla_X = std::hypot(dX_dr, dX_dphi);

		Ka(nr, naz) = diffusion_coeff(rho, x, nabla_X);
	}
    }
	
}


/*
Calculate the diffusion coefficient on at the azimuthal interface locations.
*/
static void compute_diffusion_coeff_azimuthal(t_polargrid &Density, t_polargrid &T) {

    // calcuate Kb for K(i,j/2)

	const unsigned int Nrad = Kb.get_size_radial();
    const unsigned int Naz = Kb.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nrad-1; ++nr) {
	for (unsigned int naz = 0; naz < Naz; ++naz) {
		const unsigned int naz_prev = (naz == 0 ? Naz-1 : naz - 1);

	    // average azimuthally
	    const double x = 0.5 * (T(nr, naz_prev) + T(nr, naz));
		const double rho = 0.5 * (Density(nr, naz_prev) + Density(nr, naz));

		// calculate gradient length
		// values are located on the azimuthal interface
		const double Router = Ra[nr + 1];
		const double Rinner = Ra[nr - 1];
		const double xouter = 0.5 * (T(nr + 1, naz_prev) + T(nr + 1, naz));
		const double xinner = 0.5 * (T(nr - 1, naz_prev) + T(nr - 1, naz));
	    const double dX_dr = (xouter - xinner) / (Router - Rinner);
        const double dX_dphi = InvRmed[nr] * (T(nr, naz) - T(nr, naz_prev)) / dphi;

	    const double nabla_X = std::hypot(dX_dr, dX_dphi);

	    Kb(nr, naz) = diffusion_coeff(rho, x, nabla_X);
	}
    }
}


/*
Calculate the matrix elements resulting from the diffusion equation discretization.
This sets up the linear system that is solved later.
*/
static void calculate_matrix_elements(t_polargrid &Density, const double dt) {

	// TODO: check whether we need gamma_eff here
    const double c_v = constants::R / (parameters::MU * (parameters::ADIABATICINDEX - 1.0));

	const unsigned int Nrad = A.get_size_radial();
    const unsigned int Naz = A.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nrad - 1; ++nr) {
	for (unsigned int naz = 0; naz < Naz; ++naz) {
		
		const double rho = Density(nr, naz);
		// this refers to the common factor in Müller 2013 Ph.D. thesis (A.1.10) adjusted for 3D FLD in midplane
		double common_factor;
		if (radiative_diffusion_test_2d) {
			common_factor = -dt;
		} else {
			common_factor = -dt / (rho * c_v);
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

/*
Communicate values among adjacent nodes during the matrix solver loop.
This communication function only uses one overlap cell 
which makes it different from the standard communication function used
in the hydro solver.
*/
static void communicate_parallelization_boundaries(t_polargrid &Temperature) {

	const unsigned int Nrad = Temperature.Nrad;
	const unsigned int Naz = NAzimuthal;

    const unsigned int Noc = cpuoverlap * NAzimuthal; // Number of overlap cells

	const unsigned int start_active = CPUOVERLAP*Naz;
	const unsigned int end_active = (Nrad - CPUOVERLAP)*Naz;
	const unsigned int last_active = end_active - Noc;

	const unsigned int start_inner_boundary = (CPUOVERLAP-cpuoverlap)*Naz;
	const unsigned int start_outer_boundary = end_active;

	// communicate with other nodes
	memcpy(SendInnerBoundary, Temperature.Field + start_active, Noc * sizeof(double));
	memcpy(SendOuterBoundary, Temperature.Field + last_active, Noc * sizeof(double));

	MPI_Request req1, req2, req3, req4;

	if (CPU_Rank % 2 == 0) {
	    if (CPU_Rank != 0) {
		MPI_Isend(SendInnerBoundary, Noc,
			  MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
		MPI_Irecv(RecvInnerBoundary, Noc,
			  MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req2);
	    }
	    if (CPU_Rank != CPU_Highest) {
		MPI_Isend(SendOuterBoundary, Noc,
			  MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req3);
		MPI_Irecv(RecvOuterBoundary, Noc,
			  MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req4);
	    }
	} else {
	    if (CPU_Rank != CPU_Highest) {
		MPI_Irecv(RecvOuterBoundary, Noc,
			  MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req3);
		MPI_Isend(SendOuterBoundary, Noc,
			  MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req4);
	    }
	    if (CPU_Rank != 0) {
		MPI_Irecv(RecvInnerBoundary, Noc,
			  MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
		MPI_Isend(SendInnerBoundary, Noc,
			  MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req2);
	    }
	}

	// copy values revceived from inner core to the inner local boundary
	if (CPU_Rank != 0) {
	    MPI_Wait(&req1, &global_MPI_Status);
	    MPI_Wait(&req2, &global_MPI_Status);
	    memcpy(Temperature.Field + start_inner_boundary, RecvInnerBoundary, Noc * sizeof(double));
	}

	// copy values received from outer core to the outer local boundary
	if (CPU_Rank != CPU_Highest) {
	    MPI_Wait(&req3, &global_MPI_Status);
	    MPI_Wait(&req4, &global_MPI_Status);
	    memcpy(Temperature.Field + start_outer_boundary, RecvOuterBoundary, Noc * sizeof(double));
	}
}

/*
Boundary condition to be applied inside the matrix solver loop.
*/
static void boundary_inside_SOR(t_polargrid &T) {

	const double minval = parameters::minimum_temperature;

    const unsigned int Naz = T.get_size_azimuthal();
    const unsigned int nrad_last = T.get_max_radial();

    const bool at_outer_boundary = CPU_Rank == CPU_Highest;
    const bool outer_boundary_open = parameters::boundary_outer == parameters::boundary_condition_open;
	if (at_outer_boundary && outer_boundary_open) {
		// set temperature to T_min in outermost cells
		for (unsigned int naz = 0; naz < Naz; ++naz) {
			T(nrad_last, naz) = minval;
		}
	}

    const bool at_inner_boundary = CPU_Rank == 0;
    const bool inner_boundary_open = parameters::boundary_inner == parameters::boundary_condition_open;
	if (at_inner_boundary && inner_boundary_open) {
        // set temperature to T_min in innermost cells
		for (unsigned int naz = 0; naz < Naz; ++naz) {
			T(0, naz) = minval;
		}
	}
}

/* 
Check whether a radial cell index belongs to an active cell.
An active cell is neither a ghost cell nor an overlap cell.
*/
static inline bool is_active_cell(const unsigned int nr) {
    const unsigned int nrad_last = NRadial - 1;
    const unsigned int nghost_right = (CPU_Rank == CPU_Highest) ? GHOSTCELLS_B : CPUOVERLAP;
    const unsigned int nghost_left = (CPU_Rank == 0) ? GHOSTCELLS_B : CPUOVERLAP;
    
    const bool isnot_ghostcell_rank_0 = nr > nghost_left;
    const bool isnot_ghostcell_rank_highest = nr < (nrad_last - nghost_right);

    return isnot_ghostcell_rank_0 && isnot_ghostcell_rank_highest;
}


/*
Set temperature floor and ceiling.
*/
static inline void floorceil_single(double &T) {
	if(T < parameters::minimum_temperature){
		T = parameters::minimum_temperature;
	}

	if(T > parameters::maximum_temperature){
		logging::print(LOG_INFO "max temp inside FLD SOR loop %e\n", T);
		T = parameters::maximum_temperature;
	}
}


/*
Solve the linear system using succesive over-relaxation.
*/
static void SOR(t_polargrid &T) {

    const unsigned int Naz = T.get_size_azimuthal();

    static unsigned int old_iterations = radiative_diffusion_max_iterations;
    static int direction = 1;
    static double omega = radiative_diffusion_omega;

    unsigned int iterations = 0;
	// track change of differences
	double absolute_norm = std::numeric_limits<double>::max();
	double avg_absolute_change = std::numeric_limits<double>::max();
	double last_avg_absolute_norm = 0.0;

    // do SOR
	const unsigned int maxiter = radiative_diffusion_max_iterations;

    while ((avg_absolute_change > tolerance) && (maxiter > iterations)) {

    boundary_inside_SOR(T);

	absolute_norm = 0.0;

	// Each thread should get a single consecutive range of rings just as in the MPI parallelization.
	// We first calculate the size of these consecutive ring ranges (chunk_size) and then tell openmp
	// to give every thread chunk_size consecutive rings with the schedule(dynamic, chunk_size) directive.
	static const unsigned int chunk_size = std::ceil((float)(nstop-nstart+1)/Thread_Number);
	#pragma omp parallel for reduction(+ : absolute_norm) schedule(dynamic, chunk_size)
	for (unsigned int nr = nstart; nr < nstop; ++nr) {
		for (unsigned int naz = 0; naz < Naz; ++naz) {

		const double old_value = T(nr, naz);
        const unsigned int naz_last = T.get_max_azimuthal();
		const unsigned int naz_next = (naz == naz_last ? 0 : naz + 1);
		const unsigned int naz_prev = (naz == 0 ? naz_last : naz - 1);

		T(nr, naz) =
		    (1.0 - omega) * T(nr, naz) -
		    omega / B(nr, naz) *
			(A(nr, naz) * T(nr - 1, naz) +
			 C(nr, naz) * T(nr + 1, naz) +
			 D(nr, naz) * T(nr, naz_prev) +
			 E(nr, naz) * T(nr, naz_next) - Told(nr, naz));

		// ensure minimum and maximum temperature
		floorceil_single(T(nr, naz));

		// only non ghostcells to norm and don't count overlap cell's
		// twice
        if (is_active_cell(nr)) {
		    absolute_norm += std::pow(old_value - T(nr, naz), 2);
		}
	    }
	}

    const double tmp = absolute_norm;
	MPI_Allreduce(&tmp, &absolute_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	// calculate the absolute change averaged over all cells
	const double avg_absolute_norm = std::sqrt(absolute_norm) / (GlobalNRadial * NAzimuthal);

	avg_absolute_change = fabs(avg_absolute_norm - last_avg_absolute_norm);
	last_avg_absolute_norm = avg_absolute_norm;
	// printf("iter %u, norm_change = %e\n", iterations, norm_change);
	iterations++;

    communicate_parallelization_boundaries(T);

    } // END SOR

    if (iterations == radiative_diffusion_max_iterations) {
	logging::print_master(
	    LOG_WARNING
	    "Maximum iterations (%u) reached in radiative_diffusion (omega = %lg). Norm is %lg with a last change of %lg.\n",
	    radiative_diffusion_max_iterations, omega,
	    absolute_norm, avg_absolute_change);
    }

    // adapt omega
    if (old_iterations < iterations) {
	direction *= -1;
    }

    if (radiative_diffusion_omega_auto_enabled) {
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

/*
Use radiation temperature to update internal energy.
The underlying assumption is that gas temperature 
and radiation temperature equilibrate instantaneously.
Thus, we can simply use the radiation temperature in place of the gas temperature.
*/
static void update_energy(t_polargrid &Energy, t_polargrid &Sigma, t_polargrid &T) {
	
    const unsigned int Nrad = T.get_size_radial();
    const unsigned int Naz = T.get_size_azimuthal();

    const double c_v = constants::R / (parameters::ADIABATICINDEX - 1.0) / parameters::MU;

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nrad - 1; ++nr) {
        for (unsigned int naz = 0; naz < Naz; ++naz) {
            Energy(nr, naz) = c_v * T(nr, naz) * Sigma(nr, naz);
        }
    }
}

/*
Set every value in the grid to a constant.
*/
static void set_constant_value(t_polargrid &X, const double constant){
	
	logging::print(LOG_INFO "Setting %s to constant value %e\n", X.get_name(), constant);

	const unsigned int Nrad = X.get_size_radial();
	const unsigned int Naz = X.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nrad; ++nr) {
        for (unsigned int naz = 0; naz < Naz; ++naz) {
            X(nr, naz) = constant;
        }
    }

}

/*
Perform the matrix multiplication with the solution vector to check agains old values
and print the difference.

This function can be used to check whether the tolerance is low enough
such that the solution converges.
*/
static void check_solution(t_polargrid &T) {
	
	double diff = 0;
	double diff_abs = 0;
	double xmax = 0;

	const unsigned int Naz = T.get_size_azimuthal();

	#pragma omp parallel for collapse(2) reduction(+ : diff, diff_abs, xmax)
	for (unsigned int nr = nstart; nr < nstop; ++nr) {
		for (unsigned int naz = 0; naz < Naz; ++naz) {

        const unsigned int naz_last = T.get_max_azimuthal();
		const unsigned int naz_next = (naz == naz_last ? 0 : naz + 1);
		const unsigned int naz_prev = (naz == 0 ? naz_last : naz - 1);

		const double xold = Told(nr, naz); // constant vector in matrix equation
		const double xnew =
			A(nr, naz) * T(nr - 1, naz) +
			B(nr, naz) * T(nr, naz) +
			C(nr, naz) * T(nr + 1, naz) +
			D(nr, naz) * T(nr, naz_prev) +
			E(nr, naz) * T(nr, naz_next);

		if (is_active_cell(nr)) {
			diff += xnew - xold;
			diff_abs += std::abs(xnew - xold);
			if (xnew > xmax) {
				xmax = xnew;
			}
		}
		} // azimuthal loop
	} // radial loop


    double tmp = diff;
	MPI_Allreduce(&tmp, &diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	printf("Sum of differences in RD solver : %e\n", diff);

	tmp = diff_abs;
	MPI_Allreduce(&tmp, &diff_abs, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	printf("Sum of |differences| in RD solver : %e\n", diff_abs);

	tmp = xmax;
	MPI_Allreduce(&tmp, &xmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	printf("Max value in solution vector : %e\n", xmax);
}



/*
This test checks the implementation of the diffusion process
and the linear system solver.

Instead of the radiative temperature, Trad is the radiative energy in this case
and this quantity is diffused with a constant diffusion coefficient.

The equation solved is de/dt = nabla * K nabla e

The soluction can then be used to compare against an analytical solution.
*/
static void run_2d_diffusion_test(t_data &data, const double dt) {
	const unsigned int Nrad = A.get_size_radial();
	const unsigned int Naz = A.get_size_azimuthal();
	t_polargrid x;
	x.set_size(Nrad, Naz);
	x.set_name("f_FLD2Dtest");

	// load initial condition for diffusion test
	const std::string filename_in = output::outdir + "f_FLD2Dtest_input.dat";
	x.read2D(filename_in.c_str());

	// set constant midplane density
	// const double rho = radiative_diffusion_test_2d_density;
	const double rho = 1.0;
	auto &Density = data[t_data::RHO];
	set_constant_value(Density, rho);

	const double K = radiative_diffusion_test_2d_K;
	const unsigned int steps = radiative_diffusion_test_2d_steps;

	logging::print(LOG_INFO "Running 2D FLD test with K=%e and %u steps at dt=%e\n", K, steps, dt);

	set_constant_value(Ka, K);
	set_constant_value(Kb, K);

	calculate_matrix_elements(Density, dt);
	copy_polargrid(Told, x);

	for (unsigned int i = 0; i < steps; ++i) {
		printf("Step %u\n", i);
		SOR(x);
		copy_polargrid(Told, x);
		check_solution(x);
	}

	const std::string filename_out = output::outdir + "f_FLD2Dtest_output.dat";
	x.write2D(filename_out);

	handle_output();

	die("Everything is done for the 2D FLD test.\n");
}





/*
Perform radiation transport in the one-fluid FLD approximation.

Refer to Kley & Lin (1996) doi:10.1086/177115 for a description of the approximation.
Instead of evolving a system of coupled eqations for the radiative energy and the internal energy,
it is assumed that the contribution from the radiation energy is negligible and the diffusion term
is added to the equation for the internal energy directly.
Then, the variable that is evolved is the temperature directly.

It sould be noted that this is only a valid assumption in the optically thick case because 
in the optically thin case the radiation energy is not negligible anymore.
*/
static void one_fluid_fld(t_data &data, const double dt) {
	auto &Density = data[t_data::RHO];

	// calculate diffusion coefficient
	compute_diffusion_coeff_radial(Density, Trad);
    radial_boundary_diffusion_coeff();
    compute_diffusion_coeff_azimuthal(Density, Trad);

	// setup linear equation and solve the diffusion equation
    calculate_matrix_elements(Density, dt);
	copy_polargrid(Told, Trad);

	SOR(Trad);
}


/*
Perform radiative diffusion using the one-fluid FLD method.

This method solves a linear equation system resulting 
from a diffusion equation for the temperature.
The system is solved using the SOR method.
*/
void radiative_diffusion(t_data &data, const double current_time, const double dt)
{
    auto &Sigma = data[t_data::SIGMA];
    auto &Energy = data[t_data::ENERGY];
	auto &Tgas = data[t_data::TEMPERATURE];

    set_minimum_energy(data);
	
    // update temperature, soundspeed and aspect ratio
	compute_temperature(data);
	// floor_T(data[t_data::TEMPERATURE]);
	compute_sound_speed(data, current_time);
	// compute_scale_height(data, current_time); done in compute::midplane_density
	compute::midplane_density(data, current_time);

	if (radiative_diffusion_test_2d) {
		run_2d_diffusion_test(data, dt);
		return;
	}

	// instantaneously equilibrate photon gas temperature to ideal gas temperature
	copy_polargrid(Trad, Tgas);

	one_fluid_fld(data, dt);

	if (radiative_diffusion_check_solution) {
		check_solution(Trad);
	}

	// instantaneously equilibrate ideal gas temperature to photon gas temperature
	copy_polargrid(Tgas, Trad);

	update_energy(Energy, Sigma, Tgas);
	// ensure energy is consistent with min/max temperature
	SetTemperatureFloorCeilValues(data, __FILE__, __LINE__);
}


}