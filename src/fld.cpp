#include "fld.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "parameters.h"
#include "nongnu.h"
#include "constants.h"
#include "units.h"
#include "opacity.h"
#include "pvte_law.h"
#include "SourceEuler.h"
#include "boundary_conditions.h"
#include "logging.h"


namespace fld {

// Radiative transport using Flux-Limited-Diffusion
// See Ph.D. thesis of Tobias W.A. Müller 2013, Appendix A for the equations.
// The document can be found at http://nbn-resolving.de/urn:nbn:de:bsz:21-opus-72189

t_polargrid Ka; // diffusion coefficient on radial boundaries
t_polargrid Kb; // diffusion coefficient on azimuthal boundaries
t_polargrid A, B, C, D, E; // Matrix elements of linear equation
t_polargrid Told; // Intermediate store for old values.

double *SendInnerBoundary;
double *SendOuterBoundary;
double *RecvInnerBoundary;
double *RecvOuterBoundary;

static inline double flux_limiter(const double R)
{
	#ifdef CONSTANT_FLD_FLUXLIMITER
	return 1.0/3;
	#else
    if (R <= 2) {
	return 2.0 / (3 + std::sqrt(9 + 10 * std::pow(R, 2)));
    } else {
	return 10.0 / (10 * R + 9 + std::sqrt(180 * R + 81));
    }
	#endif
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

	Told.set_scalar(true);
	Told.set_size(Nrad, Naz);

    SendInnerBoundary = (double *)malloc(Naz * CPUOVERLAP * sizeof(double));
    SendOuterBoundary = (double *)malloc(Naz * CPUOVERLAP * sizeof(double));
    RecvInnerBoundary = (double *)malloc(Naz * CPUOVERLAP * sizeof(double));
    RecvOuterBoundary = (double *)malloc(Naz * CPUOVERLAP * sizeof(double));
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

	    const double minimum_energy =
		Tmin * Sigma(1, naz) / mu * constants::R / (gamma_eff - 1.0);

	    Energy(0, naz) = minimum_energy;
	}

	if (CPU_Rank == CPU_Highest &&
	    parameters::boundary_outer == parameters::boundary_condition_open) {
	    const double Tmin = parameters::minimum_temperature; 

	    const double mu = pvte::get_mu(data, nr_max - 1, naz);
	    const double gamma_eff = pvte::get_gamma_eff(data, nr_max - 1, naz);
	    Sigma(nr_max, naz) = Sigma(nr_max - 1, naz);

	    const double minimum_energy = Tmin * Sigma(nr_max - 1, naz) / mu *
					  constants::R / (gamma_eff - 1.0);

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

static void calculate_diffusion_coeff_radial(t_data &data) {
	// Calculate the diffusion coefficient on at the radial interface locations.

    auto &Temp = data[t_data::TEMPERATURE];
    auto &Surface_density = data[t_data::SIGMA];
    auto &Scale_height = data[t_data::SCALE_HEIGHT];

    const unsigned int Nrad = Ka.get_size_radial();
    const unsigned int Naz = Ka.get_size_azimuthal();
    const unsigned int Naz_last = Ka.get_max_azimuthal();

    // calcuate Ka for K(i/2,j)
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nrad - 1; ++nr) {
	for (unsigned int naz = 0; naz < Naz; ++naz) {
	    const unsigned int naz_next = (naz == Naz_last ? 0 : naz + 1);
	    const unsigned int naz_last = (naz == 0 ? Naz_last : naz - 1);

	    // average temperature radially
	    const double T = 0.5 * (Temp(nr - 1, naz) + Temp(nr, naz));
	    const double Sigma = 0.5 * (Surface_density(nr - 1, naz) + Surface_density(nr, naz));
	    const double H = 0.5 * (Scale_height(nr - 1, naz) + Scale_height(nr, naz));

		const double density = Sigma / (parameters::density_factor * H);
		const double kappa = get_opacity(T, density);
	    const double denom = 1.0 / (Sigma * kappa);

		// calculate temperature gradient length
	    const double dT_dr = (Temp(nr, naz) - Temp(nr - 1, naz)) * InvDiffRmed[nr];
        const double Tnext = 0.5 * (Temp(nr - 1, naz_next) + Temp(nr, naz_next));
        const double Tlast = 0.5 * (Temp(nr - 1, naz_last) +	Temp(nr, naz_last));
	    const double dT_dphi = InvRinf[nr] * (Tnext-Tlast) / (2.0 * dphi);

	    const double nabla_T = std::sqrt(std::pow(dT_dr, 2) + std::pow(dT_dphi, 2));

		// Levermore & Pomraning 1981
	    // R = 4 |nabla T\/T * lrad
		// lrad = 1/(rho kappa), photon mean free path
		const double lrad = 1  / (density * kappa);
	    const double R = 4.0 * nabla_T / T * lrad;

	    const double lambda = flux_limiter(R);

        const double sigRad = constants::sigma.get_code_value();
		const double fsq = std::pow(parameters::density_factor,2);

	    Ka(nr, naz) = lambda * fsq * 16 * sigRad * denom * H * std::pow(T, 3);
	}
    }

	
}


static void calculate_diffusion_coeff_azimuthal(t_data &data) {
	// Calculate the diffusion coefficient on at the azimuthal interface locations.

    // calcuate Kb for K(i,j/2)
    auto &Temp = data[t_data::TEMPERATURE];
    auto &Surface_density = data[t_data::SIGMA];
    auto &Scale_height = data[t_data::SCALE_HEIGHT];

	const unsigned int Nrad = Kb.get_size_radial();
    const unsigned int Naz = Kb.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nrad-1; ++nr) {
	for (unsigned int naz = 0; naz < Naz; ++naz) {
	    // unsigned int n_azimuthal_plus = (n_azimuthal ==
	    // Kb.get_max_azimuthal() ? 0 : n_azimuthal + 1);
		const unsigned int naz_prev = (naz == 0 ? Naz-1 : naz - 1);

	    // average azimuthally
	    const double T = 0.5 * (Temp(nr, naz_prev) + Temp(nr, naz));
		const double Sigma = 0.5 * (Surface_density(nr, naz_prev) + Surface_density(nr, naz));
	    const double H = 0.5 * (Scale_height(nr, naz_prev) + Scale_height(nr, naz));

		const double density = Sigma / (parameters::density_factor * H);
		const double kappa = get_opacity(T, density);
	    const double denom = 1.0 / (Sigma * kappa);

		// calculate temperature gradient length
		// values are located on the azimuthal interface
		const double Router = Ra[nr + 1];
		const double Rinner = Ra[nr - 1];
		const double Touter = 0.5 * (Temp(nr + 1, naz_prev) + Temp(nr + 1, naz));
		const double Tinner = 0.5 * (Temp(nr - 1, naz_prev) + Temp(nr - 1, naz));
	    const double dT_dr = (Touter - Tinner) / (Router - Rinner);
        const double dT_dphi = InvRmed[nr] * (Temp(nr, naz) - Temp(nr, naz_prev)) / dphi;

	    const double nabla_T = std::sqrt(std::pow(dT_dr, 2) + std::pow(dT_dphi, 2));

	    // Levermore & Pomraning 1981
	    // R = 4 |nabla T\/T * lrad
		// lrad = 1/(rho kappa), photon mean free path
		const double lrad = 1  / (density * kappa);
	    const double R = 4.0 * nabla_T / T * lrad;

	    const double lambda = flux_limiter(R);

        const double sigRad = constants::sigma.get_code_value();
		const double fsq = std::pow(parameters::density_factor,2);

	    Kb(nr, naz) = lambda * fsq * 16 * sigRad * denom * H * std::pow(T, 3);
	}
    }
}

static void calculate_matrix_elements(t_data &data, const double dt) {
	// calculate the matrix elements of the linear equation 

    auto &Temperature = data[t_data::TEMPERATURE];
    auto &Surface_density = data[t_data::SIGMA];
	auto &Scale_height = data[t_data::SCALE_HEIGHT];

	// TODO: check whether we need gamma_eff here
    const double c_v = constants::R / (parameters::MU * (parameters::ADIABATICINDEX - 1.0));
	const unsigned int Nrad = Temperature.get_size_radial();
    const unsigned int Naz = Temperature.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nrad - 1; ++nr) {
	for (unsigned int naz = 0; naz < Naz; ++naz) {
		const double Sigma = Surface_density(nr, naz);
		const double H = Scale_height(nr, naz);
		// this refers to the common factor in Müller 2013 Ph.D. thesis (A.1.10)
	    const double common_factor = -dt * H / (Sigma * c_v);

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

	    Told(nr, naz) = Temperature(nr, naz);

	    /*double energy_change = dt*data[t_data::QPLUS](n_radial,
	    n_azimuthal)
		- dt*data[t_data::QMINUS](n_radial, n_azimuthal)
		- dt*data[t_data::P_DIVV](n_radial, n_azimuthal);

	    double temperature_change =
		MU/R*(parameters::ADIABATICINDEX-1.0)*energy_change/Sigma(n_radial,n_azimuthal);
	    Told(n_radial, n_azimuthal) += temperature_change;

	    if (Told(n_radial, n_azimuthal) <
	    parameters::minimum_temperature*units::temperature.get_inverse_cgs_factor())
		{ Temperature(n_radial, n_azimuthal) =
	    parameters::minimum_temperature*units::temperature.get_inverse_cgs_factor();
	    }
	    */
	}
    }

}

// Communicate the temperature values among adjacent nodes
static void communicate_parallelization_boundaries(t_polargrid &Temperature) {

    const int l = CPUOVERLAP * NAzimuthal;
    const int oo = (Temperature.Nrad - CPUOVERLAP) * NAzimuthal;
    const int o = (Temperature.Nrad - 2 * CPUOVERLAP) * NAzimuthal;

	// communicate with other nodes
	memcpy(SendInnerBoundary, Temperature.Field + l, l * sizeof(double));
	memcpy(SendOuterBoundary, Temperature.Field + o, l * sizeof(double));

	MPI_Request req1, req2, req3, req4;

	if (CPU_Rank % 2 == 0) {
	    if (CPU_Rank != 0) {
		MPI_Isend(SendInnerBoundary, NAzimuthal * CPUOVERLAP,
			  MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
		MPI_Irecv(RecvInnerBoundary, NAzimuthal * CPUOVERLAP,
			  MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req2);
	    }
	    if (CPU_Rank != CPU_Highest) {
		MPI_Isend(SendOuterBoundary, NAzimuthal * CPUOVERLAP,
			  MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req3);
		MPI_Irecv(RecvOuterBoundary, NAzimuthal * CPUOVERLAP,
			  MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req4);
	    }
	} else {
	    if (CPU_Rank != CPU_Highest) {
		MPI_Irecv(RecvOuterBoundary, NAzimuthal * CPUOVERLAP,
			  MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req3);
		MPI_Isend(SendOuterBoundary, NAzimuthal * CPUOVERLAP,
			  MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req4);
	    }
	    if (CPU_Rank != 0) {
		MPI_Irecv(RecvInnerBoundary, NAzimuthal * CPUOVERLAP,
			  MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
		MPI_Isend(SendInnerBoundary, NAzimuthal * CPUOVERLAP,
			  MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req2);
	    }
	}

	if (CPU_Rank != 0) {
	    MPI_Wait(&req1, &global_MPI_Status);
	    MPI_Wait(&req2, &global_MPI_Status);
	    memcpy(Temperature.Field, RecvInnerBoundary, l * sizeof(double));
	}

	if (CPU_Rank != CPU_Highest) {
	    MPI_Wait(&req3, &global_MPI_Status);
	    MPI_Wait(&req4, &global_MPI_Status);
	    memcpy(Temperature.Field + oo, RecvOuterBoundary,
		   l * sizeof(double));
	}
}

static void boundary_T_SOR(t_polargrid &Temperature) {

    const unsigned int Naz = Temperature.get_size_azimuthal();
    const unsigned int nrad_last = Temperature.get_max_radial();

    const bool at_outer_boundary = CPU_Rank == CPU_Highest;
    const bool outer_boundary_open = parameters::boundary_outer == parameters::boundary_condition_open;
	if (at_outer_boundary && outer_boundary_open) {
		// set temperature to T_min in outermost cells
		for (unsigned int naz = 0; naz < Naz; ++naz) {
			Temperature(nrad_last, naz) = parameters::minimum_temperature;
		}
	}

    const bool at_inner_boundary = CPU_Rank == 0;
    const bool inner_boundary_open = parameters::boundary_inner == parameters::boundary_condition_open;
	if (at_inner_boundary && inner_boundary_open) {
        // set temperature to T_min in innermost cells
		for (unsigned int naz = 0; naz < Naz; ++naz) {
			Temperature(0, naz) = parameters::minimum_temperature;
		}
	}
	// boundary_conditions::apply_boundary_condition(data, current_time, 0.0, false);
}

static bool is_active(const unsigned int nr) {
	// Check whether a radial cell index belongs to an active cell.
	// An active cell is neither a ghost cell nor an overlap cell.
    const unsigned int nrad_last = NRadial - 1;
    const unsigned int nghost_right = (CPU_Rank == CPU_Highest) ? GHOSTCELLS_B : CPUOVERLAP;
    const unsigned int nghost_left = (CPU_Rank == 0) ? GHOSTCELLS_B : CPUOVERLAP;
    
    const bool isnot_ghostcell_rank_0 = nr > nghost_left;
    const bool isnot_ghostcell_rank_highest = nr < (nrad_last - nghost_right);

    return isnot_ghostcell_rank_0 && isnot_ghostcell_rank_highest;
}

static void SOR(t_data &data) {

    auto &Temperature = data[t_data::TEMPERATURE];
    
    const unsigned int Nrad = Temperature.get_size_radial();
    const unsigned int Naz = Temperature.get_size_azimuthal();

    static unsigned int old_iterations = parameters::radiative_diffusion_max_iterations;
    static int direction = 1;
    static double omega = parameters::radiative_diffusion_omega;

    unsigned int iterations = 0;
	double absolute_norm = std::numeric_limits<double>::max();
	double norm_change = std::numeric_limits<double>::max();

    // do SOR
	const double tolerance = parameters::radiative_diffusion_tolerance;

    boundary_T_SOR(Temperature);

	norm_change = absolute_norm;
	absolute_norm = 0.0;
#pragma omp parallel for collapse(2) reduction(+ : absolute_norm)
	for (unsigned int nr = 1; nr < Nrad - 1; ++nr) {
		for (unsigned int naz = 0; naz < Naz; ++naz) {

		const double old_value = Temperature(nr, naz);
        const unsigned int naz_last = Temperature.get_max_azimuthal();
		const unsigned int naz_next = (naz == naz_last ? 0 : naz + 1);
		const unsigned int naz_prev = (naz == 0 ? naz_last : naz - 1);

		Temperature(nr, naz) =
		    (1.0 - omega) * Temperature(nr, naz) -
		    omega / B(nr, naz) *
			(A(nr, naz) * Temperature(nr - 1, naz) +
			 C(nr, naz) * Temperature(nr + 1, naz) +
			 D(nr, naz) * Temperature(nr, naz_prev) +
			 E(nr, naz) * Temperature(nr, naz_next) - Told(nr, naz));

        // Set temperature floor and ceiling
        if(Temperature(nr, naz) < parameters::minimum_temperature){
            Temperature(nr, naz) = parameters::minimum_temperature;
        }

        if(Temperature(nr, naz) > parameters::maximum_temperature){
            logging::print(LOG_INFO "max temp inside FLD SOR loop at %d, %d, %e\n", nr, naz, Temperature(nr, naz));
            Temperature(nr, naz) = parameters::maximum_temperature;
        }

		// only non ghostcells to norm and don't count overlap cell's
		// twice
        if (is_active(nr)) {
		    absolute_norm += std::pow(old_value - Temperature(nr, naz), 2);
		}
	    }
	}

    const double tmp = absolute_norm;
	MPI_Allreduce(&tmp, &absolute_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	absolute_norm = std::sqrt(absolute_norm) / (GlobalNRadial * NAzimuthal);

	norm_change = fabs(absolute_norm - norm_change);
	iterations++;

    communicate_parallelization_boundaries(Temperature);

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

static void update_energy(t_data &data) {

    auto &Temperature = data[t_data::TEMPERATURE];
    auto &Sigma = data[t_data::SIGMA];
    auto &Energy = data[t_data::ENERGY];

    const unsigned int Nrad = Temperature.get_size_radial();
    const unsigned int Naz = Temperature.get_size_azimuthal();

    const double factor = constants::R / (parameters::ADIABATICINDEX - 1.0) / parameters::MU;
    // compute energy from temperature
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nrad - 1; ++nr) {
        for (unsigned int naz = 0; naz < Naz; ++naz) {
            Energy(nr, naz) = factor * Temperature(nr, naz) * Sigma(nr, naz);
        }
    }
}

void radiative_diffusion(t_data &data, const double current_time, const double dt)
{
    set_minimum_energy(data);

    // update temperature, soundspeed and aspect ratio
	compute_temperature(data);
	compute_sound_speed(data, current_time);
	compute_scale_height(data, current_time);

	// calculate diffusion coefficient
	calculate_diffusion_coeff_radial(data);
    radial_boundary_diffusion_coeff();
    calculate_diffusion_coeff_azimuthal(data);

	// setup linear equation and solve the diffusion equation
    calculate_matrix_elements(data, dt);
    SOR(data);

    update_energy(data);

    // ensure energy is consistent with min/max temperature
	SetTemperatureFloorCeilValues(data, __FILE__, __LINE__);
}


}