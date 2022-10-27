#include "sts.h"

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <vector>

#include "LowTasks.h"
#include "Pframeforce.h"
#include "SideEuler.h"
#include "SourceEuler.h"
#include "Theo.h"
#include "TransportEuler.h"
#include "boundary_conditions.h"
#include "commbound.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "nongnu.h"
#include "parameters.h"
#include "quantities.h"
#include "stress.h"
#include "units.h"
#include "util.h"
#include "viscosity/viscosity.h"
#include <cstring>

static void StsStep(t_data &data, double dt)
{
    // Now we can update energy with source terms
    for (int n_radial = 0; n_radial <= data[t_data::ENERGY].get_max_radial();
	 ++n_radial) {
	for (int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal();
	     ++n_azimuthal) {

	    data[t_data::ENERGY](n_radial, n_azimuthal) +=
		dt * data[t_data::QPLUS](n_radial, n_azimuthal);
	}
    }
}

/**
	In this substep we add the articifial viscous pressure source terms.
	Shocks are spread over CVNR zones: von Neumann-Richtmyer viscosity
   constant; Beware of misprint in Stone and Norman's paper : use C2^2 instead
   of C2
*/
static void StsStep2(t_data &data, double dt)
{

    if (parameters::Adiabatic) {
	// clear up all Qplus terms
	data[t_data::QPLUS].clear();
    }

	/// Do not Apply sub keplerian boundary for boundary conditions that set
	/// Vphi themselves
	const bool add_kep_inner =
	(parameters::boundary_inner !=
	parameters::boundary_condition_center_of_mass_initial) &&
	(parameters::boundary_inner !=
	 parameters::boundary_condition_evanescent) &&
	(parameters::boundary_inner !=
	 parameters::boundary_condition_boundary_layer) &&
	(parameters::boundary_inner !=
	 parameters::boundary_condition_precribed_time_variable) &&
	(!parameters::domegadr_zero);

	if (add_kep_inner) {
	ApplySubKeplerianBoundaryInner(data[t_data::V_AZIMUTHAL]);
	}

	if ((parameters::boundary_outer !=
	 parameters::boundary_condition_center_of_mass_initial) &&
	(parameters::boundary_outer !=
	 parameters::boundary_condition_zero_gradient) &&
	(parameters::boundary_outer !=
	 parameters::boundary_condition_evanescent) &&
	(parameters::boundary_outer !=
	 parameters::boundary_condition_boundary_layer) &&
	(parameters::boundary_outer !=
	 parameters::boundary_condition_precribed_time_variable) &&
	(!parameters::massoverflow) && (!parameters::domegadr_zero)) {
	ApplySubKeplerianBoundaryOuter(data[t_data::V_AZIMUTHAL],
					   add_kep_inner);
	}

    if (parameters::artificial_viscosity ==
	parameters::artificial_viscosity_SN) {

	// calculate q_r and q_phi
	for (int n_radial = 0; n_radial <= data[t_data::Q_R].get_max_radial();
	     ++n_radial) {
	    for (int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::Q_R].get_max_azimuthal();
		 ++n_azimuthal) {
		double dv_r =
		    data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) -
		    data[t_data::V_RADIAL](n_radial, n_azimuthal);
		if (dv_r < 0.0) {
		    data[t_data::Q_R](n_radial, n_azimuthal) =
			std::pow(parameters::artificial_viscosity_factor, 2) *
			data[t_data::SIGMA](n_radial, n_azimuthal) *
			std::pow(dv_r, 2);
		} else {
		    data[t_data::Q_R](n_radial, n_azimuthal) = 0.0;
		}

		double dv_phi =
		    data[t_data::V_AZIMUTHAL](
			n_radial,
			n_azimuthal ==
				data[t_data::V_AZIMUTHAL].get_max_azimuthal()
			    ? 0
			    : n_azimuthal + 1) -
		    data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal);
		if (dv_phi < 0.0) {
		    data[t_data::Q_PHI](n_radial, n_azimuthal) =
			std::pow(parameters::artificial_viscosity_factor, 2) *
			data[t_data::SIGMA](n_radial, n_azimuthal) *
			std::pow(dv_phi, 2);
		} else {
		    data[t_data::Q_PHI](n_radial, n_azimuthal) = 0.0;
		}
	    }
	}

	// If gas disk is adiabatic, we add artificial viscosity as a source
	// term for advection of thermal energy polargrid
	if (parameters::Adiabatic) {

		if (parameters::artificial_viscosity_dissipation) {
		for (int n_radial = 0;
			 n_radial <= data[t_data::QPLUS].get_max_radial();
			 ++n_radial) {
			const double dxtheta = dphi * Rmed[n_radial];
			const double invdxtheta = 1.0 / dxtheta;
			for (int n_azimuthal = 0;
			 n_azimuthal <= data[t_data::QPLUS].get_max_azimuthal();
			 ++n_azimuthal) {
			data[t_data::QPLUS](n_radial, n_azimuthal) =
				-data[t_data::Q_R](n_radial, n_azimuthal) *
				(data[t_data::V_RADIAL](n_radial + 1,
							n_azimuthal) -
				 data[t_data::V_RADIAL](n_radial,
							n_azimuthal)) *
				InvDiffRsup[n_radial] -
				data[t_data::Q_PHI](n_radial, n_azimuthal) *
				(data[t_data::V_AZIMUTHAL](
					 n_radial,
					 n_azimuthal == data[t_data::V_AZIMUTHAL]
							.get_max_azimuthal()
					 ? 0
					 : n_azimuthal + 1) -
				 data[t_data::V_AZIMUTHAL](n_radial,
							   n_azimuthal)) *
				invdxtheta;
			}
		}
		}
	}

	// add artificial viscous pressure source term to v_radial
	for (int n_radial = 1;
	     n_radial <= data[t_data::V_RADIAL].get_max_radial() - 1;
	     ++n_radial) {
	    for (int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::V_RADIAL].get_max_azimuthal();
		 ++n_azimuthal) {
		// 1/Sigma dq_r/dr : Sigma is calculated as a mean value between
		// the neightbour cells
		data[t_data::V_RADIAL](n_radial, n_azimuthal) =
		    data[t_data::V_RADIAL](n_radial, n_azimuthal) -
		    dt * 2.0 /
			(data[t_data::SIGMA](n_radial, n_azimuthal) +
			 data[t_data::SIGMA](n_radial - 1, n_azimuthal)) *
			(data[t_data::Q_R](n_radial, n_azimuthal) -
			 data[t_data::Q_R](n_radial - 1, n_azimuthal)) *
			InvDiffRmed[n_radial];
	    }
	}

	// add artificial viscous pressure source term to v_azimuthal
	for (int n_radial = 0;
	     n_radial <= data[t_data::V_AZIMUTHAL].get_max_radial();
	     ++n_radial) {
	    const double dxtheta = dphi * Rmed[n_radial];
	    const double invdxtheta = 1.0 / dxtheta;
	    for (int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::V_AZIMUTHAL].get_max_azimuthal();
		 ++n_azimuthal) {
		// 1/Sigma 1/r dq_phi/dphi : Sigma is calculated as a mean value
		// between the neightbour cells
		data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) =
		    data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) -
		    dt * 2.0 /
			(data[t_data::SIGMA](n_radial, n_azimuthal) +
			 data[t_data::SIGMA](
			     n_radial,
			     n_azimuthal == 0
				 ? data[t_data::SIGMA].get_max_azimuthal()
				 : n_azimuthal - 1)) *
			(data[t_data::Q_PHI](n_radial, n_azimuthal) -
			 data[t_data::Q_PHI](
			     n_radial,
			     n_azimuthal == 0
				 ? data[t_data::Q_PHI].get_max_azimuthal()
				 : n_azimuthal - 1)) *
			invdxtheta;
	    }
	}
    }
}

static void calculateQvis(t_data &data)
{

    if (parameters::heating_viscous_enabled) {
	/* We calculate the heating source term Qplus from i=1 to max-1 */
	for (int n_radial = 1;
	     n_radial <= data[t_data::QPLUS].get_max_radial() - 1; ++n_radial) {
	    for (int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::QPLUS].get_max_azimuthal();
		 ++n_azimuthal) {
		if (data[t_data::VISCOSITY](n_radial, n_azimuthal) != 0.0) {
		    // average tau_r_phi over 4 cells
		    double tau_r_phi =
			0.25 *
			(data[t_data::TAU_R_PHI](n_radial, n_azimuthal) +
			 data[t_data::TAU_R_PHI](n_radial + 1, n_azimuthal) +
			 data[t_data::TAU_R_PHI](
			     n_radial,
			     n_azimuthal ==
				     data[t_data::TAU_R_PHI].get_max_azimuthal()
				 ? 0
				 : n_azimuthal + 1) +
			 data[t_data::TAU_R_PHI](
			     n_radial + 1,
			     n_azimuthal ==
				     data[t_data::TAU_R_PHI].get_max_azimuthal()
				 ? 0
				 : n_azimuthal + 1));

		    double qplus =
			1.0 /
			(2.0 * data[t_data::VISCOSITY](n_radial, n_azimuthal) *
			 data[t_data::SIGMA](n_radial, n_azimuthal)) *
			(std::pow(data[t_data::TAU_R_R](n_radial, n_azimuthal),
				  2) +
			 2 * std::pow(tau_r_phi, 2) +
			 std::pow(
			     data[t_data::TAU_PHI_PHI](n_radial, n_azimuthal),
			     2));
		    qplus +=
			(2.0 / 9.0) *
			data[t_data::VISCOSITY](n_radial, n_azimuthal) *
			data[t_data::SIGMA](n_radial, n_azimuthal) *
			std::pow(data[t_data::DIV_V](n_radial, n_azimuthal), 2);

		    qplus *= parameters::heating_viscous_factor;
		    data[t_data::QPLUS](n_radial, n_azimuthal) += qplus;
		}
	    }
	}

	/// We calculate the heating source term Qplus for i=max
	/*for (int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::QPLUS].get_max_azimuthal();
	     ++n_azimuthal) {
	    if (data[t_data::VISCOSITY](data[t_data::QPLUS].get_max_radial(),
					n_azimuthal) != 0.0) {
		// power-law extrapolation
		double qplus =
		    data[t_data::QPLUS](
			data[t_data::QPLUS].get_max_radial() - 1, n_azimuthal) *
		    std::exp(
			std::log(data[t_data::QPLUS](
				     data[t_data::QPLUS].get_max_radial() - 1,
				     n_azimuthal) /
				 data[t_data::QPLUS](
				     data[t_data::QPLUS].get_max_radial() - 2,
				     n_azimuthal)) *
			std::log(
			    Rmed[data[t_data::QPLUS].get_max_radial()] /
			    Rmed[data[t_data::QPLUS].get_max_radial() - 1]) /
			std::log(
			    Rmed[data[t_data::QPLUS].get_max_radial() - 1] /
			    Rmed[data[t_data::QPLUS].get_max_radial() - 2]));

		data[t_data::QPLUS](data[t_data::QPLUS].get_max_radial(),
				    n_azimuthal) += qplus;
	    }
	}*/
	/// We calculate the heating source term Qplus for i=0
	/*for (int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::QPLUS].get_max_azimuthal();
	     ++n_azimuthal) {
	    if (data[t_data::VISCOSITY](0, n_azimuthal) != 0.0) {
		// power-law extrapolation
		double qplus =
		    data[t_data::QPLUS](1, n_azimuthal) *
		    std::exp(std::log(data[t_data::QPLUS](1, n_azimuthal) /
				      data[t_data::QPLUS](2, n_azimuthal)) *
			     std::log(Rmed[0] / Rmed[1]) /
			     std::log(Rmed[1] / Rmed[2]));

		data[t_data::QPLUS](0, n_azimuthal) += qplus;
		}
	}*/
    }
}

/*
namespace Detail
{
static double constexpr sqrtNewtonRaphson(double x, double curr, double prev)
{
	return curr == prev ? curr
			: sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
}
} // namespace Detail


 * Constexpr version of the square root
 * Return value:
 *   - For a finite and non-negative value of "x", returns an approximation for
 * the square root of "x"
 *   - Otherwise, returns NaN

static double constexpr const_sqrt(double x)
{
	return x >= 0 && x < std::numeric_limits<double>::infinity()
		   ? Detail::sqrtNewtonRaphson(x, x, 0)
		   : std::numeric_limits<double>::quiet_NaN();
} */

static constexpr int STS_MAX_STEPS = 1024;
double ts[STS_MAX_STEPS];

/* ********************************************************************* */
static void STS_ComputeSubSteps(double dtex, double *tau, int ssorder)
/*
 *
 *********************************************************************** */
{
    int i;
    for (i = 0; i < ssorder; i++) {
	tau[i] = dtex / ((-1.0 + parameters::STS_NU) * std::cos(((2.0 * i + 1.0) * M_PI) /
						    (2.0 * ssorder)) +
			 1.0 + parameters::STS_NU);
    }
}

/* ********************************************************************* */
static double STS_FindRoot(double x0, double dtr, double dta)
/*
 *
 *********************************************************************** */
{
    double b, c;
    double Ns, Ns1;
    int n;

    n = 0;

    Ns = x0 + 1.0;
    Ns1 = x0;

    while (fabs(Ns - Ns1) >= 1.0e-5) {
	Ns = Ns1;
	const static double a =
	    (1.0 - std::sqrt(parameters::STS_NU)) / (1.0 + std::sqrt(parameters::STS_NU));
	const static double sqrt_sts_nu = std::sqrt(parameters::STS_NU);

	b = std::pow(a, 2.0 * Ns);
	c = (1.0 - b) / (1.0 + b);
	Ns1 =
	    Ns + (dta - dtr * Ns / (2.0 * sqrt_sts_nu) * c) /
		     (dtr / (2.0 * sqrt_sts_nu) *
		      (c - 2.0 * Ns * b * std::log(a) * (1.0 + c) / (1.0 + b)));
	n += 1;
	if (n == 128) {
	    printf("! STS_FindRoot: max number of iterations exceeded");
	    PersonalExit(1);
	}
    }
    return (Ns);
}

/* ********************************************************************* */
static double STS_CorrectTimeStep(int n0, double dta)
/*
 *
 *********************************************************************** */
{
    double dtr;

    const static double a =
	(1.0 - std::sqrt(parameters::STS_NU)) / (1.0 + std::sqrt(parameters::STS_NU));
    const static double sqrt_sts_nu = std::sqrt(parameters::STS_NU);
    const double b = std::pow(a, 2.0 * n0);
    const double c = (1.0 - b) / (1.0 + b);

    dtr = dta * 2.0 * sqrt_sts_nu / (n0 * c);
    return (dtr);
}

void Sts(t_data &data, const double current_time, const double dt)
{

    if (parameters::heating_viscous_enabled) {
	// clear up all Qplus terms
	data[t_data::QPLUS].clear();
    }

	boundary_conditions::apply_boundary_condition(data, current_time, 0.0, false);
    double dt_par;

    MPI_Allreduce(&dt_parabolic_local, &dt_par, 1, MPI_DOUBLE, MPI_MIN,
		  MPI_COMM_WORLD);

    double tau;

    int m = 0;
    int n = STS_MAX_STEPS;

    // get Number of steps
    double N = ceil(STS_FindRoot(1.0, dt_par, dt));
    n = (int)N;

    // printf("\nN = %d %.5e	%.5e	%.5e	%.5e	N = %.5e\n",n, N,
    // dt_par, dt, dt_parabolic_local, STS_FindRoot(1.0, dt_par, dt));
    if (n > STS_MAX_STEPS) {
	logging::print(LOG_ERROR "STS: the number of substeps (%d) is > %d\n",
		       N, STS_MAX_STEPS);
	PersonalExit(1);
    }

    if (n > 1) {
	dt_par = STS_CorrectTimeStep(n, dt);
	STS_ComputeSubSteps(dt_par, ts, n);
    }
    if (n == 1)
	ts[0] = dt;

    // begin Super time stepping
    while (m < n) {

	tau = ts[n - m - 1];

	StsStep2(data, tau);
	recalculate_viscosity(data, current_time);
	viscosity::compute_viscous_stress_tensor(data);
	viscosity::update_velocities_with_viscosity(
		data, tau);

	boundary_conditions::apply_boundary_condition(data, current_time, 0.0, false);

	if (parameters::Adiabatic) {

	    calculateQvis(data);
	    StsStep(data, tau);
	    SetTemperatureFloorCeilValues(data, __FILE__, __LINE__);
	}

	m++;
    }

    if (parameters::Adiabatic) {
	// clear up all Qplus terms
	data[t_data::QPLUS].clear();
    }
}
