/**
\file TransportEuler.c

Functions that handle the transport substep of a hydrodynamical time step. The
FARGO algorithm is implemented here. The transport is performed in a manner
similar to what is done for the ZEUS code (Stone & Norman, 1992), except for the
momenta transport (we define a left and right momentum for each zone, which we
declare zone centered; we then transport then normally, and deduce the new
velocity in each zone by a proper averaging).
*/

#include <cmath>
#include <cstring>

#include "polargrid.h"
#include "radialgrid.h"

#include "LowTasks.h"
#include "SideEuler.h"
#include "SourceEuler.h"
#include "TransportEuler.h"
#include "boundary_conditions.h"
#include "constants.h"
#include "global.h"
#include "parameters.h"
#include "util.h"
#include "frame_of_reference.h"

// radial momentum is split to be cell centered
static t_polargrid radial_momentum_plus;
static t_polargrid radial_momentum_minus;
static t_polargrid angular_momentum_plus;
static t_polargrid angular_momentum_minus;

// residual azimuthal velocity (needed for FARGO)
static t_polargrid v_azimuthal_res;

// average azimuthal velocity (needed for FARGO)
static t_radialgrid v_azimuthal_mean;

static PolarGrid *Work, *QRStar, *DensityStar;

static double *TempShift;
static double *dq;

// static double VMed[MAX1D];
static int Nshift[MAX1D];
static boolean NoSplitAdvection[MAX1D];
static boolean UniformTransport;
extern boolean OpenInner;

/**
	Initializes (allocates) all variables needed.
*/
void InitTransport()
{
    radial_momentum_plus.set_scalar(true);
    radial_momentum_plus.set_size(NRadial, NAzimuthal);
    radial_momentum_plus.set_name("radial_momentum_plus");
    radial_momentum_minus.set_scalar(true);
    radial_momentum_minus.set_size(NRadial, NAzimuthal);
    radial_momentum_minus.set_name("radial_momentum_minus");
    angular_momentum_plus.set_scalar(true);
    angular_momentum_plus.set_size(NRadial, NAzimuthal);
    angular_momentum_plus.set_name("angular_momentum_plus");
    angular_momentum_minus.set_scalar(true);
    angular_momentum_minus.set_size(NRadial, NAzimuthal);
    angular_momentum_minus.set_name("angular_momentum_minus");

    v_azimuthal_res.set_vector(false);
    v_azimuthal_res.set_size(NRadial, NAzimuthal);
    v_azimuthal_res.set_name("v_azimuthal_res");

    v_azimuthal_mean.set_vector(false);
    v_azimuthal_mean.set_size(NRadial);
    v_azimuthal_mean.set_name("v_azimuthal_mean");

    Work = CreatePolarGrid(NRadial, NAzimuthal, "WorkGrid");

    QRStar = CreatePolarGrid(NRadial, NAzimuthal, "QRStar");
    // TODO: This has to be checked!!!
    QRStar->set_vector(true);
    QRStar->set_size(NRadial, NAzimuthal);

    DensityStar = CreatePolarGrid(NRadial, NAzimuthal, "DensityStar");
    // TODO: This has to be checked! set_vector(true) is wrong, since there
    // cannot be flow through outermost ghost cell
    DensityStar->set_vector(true);
    DensityStar->set_size(NRadial, NAzimuthal);

    TempShift = (double *)malloc(NRadial * NAzimuthal * sizeof(double));

    dq = (double *)malloc(NRadial * NAzimuthal * sizeof(double));
}

/**
	Frees all variables allocate before by InitTransport.
*/
void FreeTransport()
{

    delete DensityStar;
    delete Work;
    delete QRStar;

    free(TempShift);
    free(dq);
}

void Transport(t_data &data, PolarGrid *Density, PolarGrid *VRadial,
		   PolarGrid *VAzimuthal, PolarGrid *Energy, const double dt)
{
    compute_momenta_from_velocities(*Density, *VRadial, *VAzimuthal);

    /* No-Alternate Directionnal Splitting */
    OneWindRad(data, Density, VRadial, Energy, dt);
    OneWindTheta(data, Density, VAzimuthal, Energy, dt);

    compute_velocities_from_momenta(*Density, *VRadial, *VAzimuthal);

	// assure minimum density after transport
	assure_minimum_value(data[t_data::SIGMA],
			 parameters::sigma_floor * parameters::sigma0);

	if (parameters::Adiabatic) {
	// assure minimum temperature after transport. it
	// is crucial the check minimum density before!
	SetTemperatureFloorCeilValues(data, __FILE__, __LINE__);
	}
}

void OneWindRad(t_data &data, PolarGrid *Density, PolarGrid *VRadial,
		PolarGrid *Energy, double dt)
{
    compute_star_radial(Density, VRadial, DensityStar, dt);

    // boundary layer:
    if (CPU_Rank == CPU_Highest) {

	if (parameters::boundary_outer ==
	    parameters::boundary_condition_boundary_layer) {
	    boundary_layer_mass_influx(DensityStar, VRadial);
	}

	// prescribed time variable boundary
	if (parameters::boundary_outer ==
	    parameters::boundary_condition_precribed_time_variable) {
	    boundary_conditions::
		boundary_condition_precribed_time_variable_outer(data,
								 DensityStar);
	}

	if (parameters::massoverflow) {
	    boundary_conditions::mass_overflow_willy(data, DensityStar, true);
	}
    }

    copy_polargrid(data[t_data::DENSITY_INT], *Density);

    VanLeerRadial(data, VRadial, &radial_momentum_plus, dt);
    VanLeerRadial(data, VRadial, &radial_momentum_minus, dt);
    VanLeerRadial(data, VRadial, &angular_momentum_plus, dt);
    VanLeerRadial(data, VRadial, &angular_momentum_minus, dt);

    if (parameters::Adiabatic) {
	VanLeerRadial(data, VRadial, Energy, dt);
    }

    VanLeerRadial(data, VRadial, Density, dt); /* MUST be the last line */
}

/* Hereafter are the new specific procedures to the fast algorithm */

/**
	Compute average azimuthal velocities
*/
void compute_average_azimuthal_velocity(t_polargrid &v_azimuthal)
{
    for (unsigned int n_radial = 0; n_radial < v_azimuthal.get_size_radial();
	 ++n_radial) {
	double v_azimuthal_sum = 0.0;

	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < v_azimuthal.get_size_azimuthal(); ++n_azimuthal) {
	    v_azimuthal_sum += v_azimuthal(n_radial, n_azimuthal);
	}

	v_azimuthal_mean(n_radial) =
	    v_azimuthal_sum / (double)v_azimuthal.get_size_azimuthal();
    }
}

/**
	Compute (azimuthal) residual velocities
*/
void compute_residual_velocity(t_polargrid &v_azimuthal)
{
    for (unsigned int n_radial = 0; n_radial < v_azimuthal.get_size_radial();
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < v_azimuthal.get_size_azimuthal(); ++n_azimuthal) {
	    v_azimuthal_res(n_radial, n_azimuthal) =
		v_azimuthal(n_radial, n_azimuthal) - v_azimuthal_mean(n_radial);
	}
    }
}

void ComputeConstantResidual(PolarGrid *VAzimuthal, double dt)
{
    int i, j, l, nr, ns;
    long nitemp;
    double *vt, *vres, Ntilde, Nround, invdt;
    nr = VAzimuthal->Nrad;
    ns = VAzimuthal->Nsec;
    vt = VAzimuthal->Field;
    vres = v_azimuthal_res.Field;
    invdt = 1.0 / dt;

    for (i = 0; i < nr; i++) {
	Ntilde = v_azimuthal_mean(i) * InvRmed[i] * dt * invdphi;
	Nround = floor(Ntilde + 0.5);
	nitemp = (long)Nround;
	Nshift[i] = (long)nitemp;
	for (j = 0; j < ns; j++) {
	    l = j + i * ns;
	    vt[l] = (Ntilde - Nround) * Rmed[i] * invdt * dphi;
	}
	if (!parameters::fast_transport) {
	    NoSplitAdvection[i] = YES;
	    for (j = 0; j < ns; j++) {
		l = j + i * ns;
		vres[l] = vt[l] + vres[l];
		vt[l] = 0.0;
	    }
	} else {
	    NoSplitAdvection[i] = NO;
	}
    }
}

void AdvectSHIFT(t_polargrid &array)
{
    int i, j, ji, l, li, nr, ns;
    double *val;
    val = array.Field;
    nr = array.Nrad;
    ns = array.Nsec;

    for (i = 0; i < nr; i++) {
	for (j = 0; j < ns; j++) {
	    ji = j - Nshift[i];
	    while (ji < 0)
		ji += ns;
	    while (ji >= ns)
		ji -= ns;
	    l = j + i * ns;
	    li = ji + i * ns;
	    TempShift[l] = val[li];
	}
	for (j = 0; j < ns; j++) {
	    l = j + i * ns;
	    val[l] = TempShift[l];
	}
    }
}

void OneWindTheta(t_data &data, PolarGrid *Density, PolarGrid *VAzimuthal,
		  PolarGrid *Energy, double dt)
{
    compute_average_azimuthal_velocity(*VAzimuthal);
    compute_residual_velocity(*VAzimuthal);
    ComputeConstantResidual(
	VAzimuthal, dt); /* Constant residual is in VAzimuthal from now on */
    UniformTransport = NO;
    QuantitiesAdvection(data, Density, &v_azimuthal_res, Energy, dt);
    UniformTransport = YES;
    QuantitiesAdvection(data, Density, VAzimuthal, Energy, dt);
    AdvectSHIFT(radial_momentum_plus);
    AdvectSHIFT(radial_momentum_minus);
    AdvectSHIFT(angular_momentum_plus);
    AdvectSHIFT(angular_momentum_minus);
    if (parameters::Adiabatic)
	AdvectSHIFT(*Energy);
    AdvectSHIFT(*Density);
}

/* End of new specific procedures to the fast algorithm */

void QuantitiesAdvection(t_data &data, PolarGrid *Density,
			 PolarGrid *VAzimuthal, PolarGrid *Energy, double dt)
{
    ComputeStarTheta(Density, VAzimuthal, DensityStar, dt);
    copy_polargrid(data[t_data::DENSITY_INT], *Density);
    VanLeerTheta(data, VAzimuthal, &radial_momentum_plus, dt);
    VanLeerTheta(data, VAzimuthal, &radial_momentum_minus, dt);
    VanLeerTheta(data, VAzimuthal, &angular_momentum_plus, dt);
    VanLeerTheta(data, VAzimuthal, &angular_momentum_minus, dt);
    if (parameters::Adiabatic)
	VanLeerTheta(data, VAzimuthal, Energy, dt);
    VanLeerTheta(data, VAzimuthal, Density, dt); /* MUST be the last line */
}

static double van_leer_lim(const double a, const double b){
	if(a*b > 0.0){
		return 2.0*a*b / (a+b);
	} else {
		return 0;
	}
}

static double minmod(const double a, const double b){
	if(a*b > 0.0)
		return std::fabs(a) < std::fabs(b) ? a : b;
	else
		return 0.0;
}

static double MC_lim(const double a, const double b){
	return minmod(0.5*(a+b), 2.0*minmod(a, b));
}

static double flux_limiter(const double a, const double b){
	switch(flux_limiter_type){
		case 0:
			return van_leer_lim(a, b);
			break;
		case 1:
			return MC_lim(a, b);
			break;
		default:
			return van_leer_lim(a, b);
			break;
	}
}

static unsigned int cell_number(const unsigned int nrad, 
								const unsigned int naz, 
								const unsigned int Naz_tot) {
    return nrad * Naz_tot + naz;
}

/**
*/
// void compute_star_radial(t_polargrid* base, t_polargrid* V_Radial,
// t_polargrid* star, double dt)
void compute_star_radial(t_polargrid *Qbase, t_polargrid *VRadial,
			 t_polargrid *QStar, double dt)
{
    unsigned int nRadial, nAzimuthal;
    double dqp, dqm;

    for (nRadial = 0; nRadial < Qbase->Nrad; ++nRadial) {
	for (nAzimuthal = 0; nAzimuthal < Qbase->Nsec; ++nAzimuthal) {
	    unsigned int cell = cell_number(nRadial, nAzimuthal, Qbase->Nsec);
	    unsigned int cellPrevRadial = cell - Qbase->Nsec;
	    unsigned int cellNextRadial = cell + Qbase->Nsec;
	    if ((nRadial == 0) || (nRadial == Qbase->Nrad - 1)) {
		dq[cell] = 0.0;
	    } else {
		dqm = (Qbase->Field[cell] - Qbase->Field[cellPrevRadial]) *
		      InvDiffRmed[nRadial];
		dqp = (Qbase->Field[cellNextRadial] - Qbase->Field[cell]) *
		      InvDiffRmed[nRadial + 1];
		dq[cell] = flux_limiter(dqp, dqm);
		}
	}
    }
    // TODO: changed to nRadial =1 because of nRadial-1
    // TODO: potential problem: Using Rmed[nRadial] - Rmed[nRadial-1] for
    // a-mesh Qties (v_rad,..) as well as for b-mesh Qties (Density,...)
    for (nRadial = 1; nRadial < Qbase->Nrad; ++nRadial) {
	for (nAzimuthal = 0; nAzimuthal < Qbase->Nsec; ++nAzimuthal) {
	    unsigned int cell = cell_number(nRadial, nAzimuthal, Qbase->Nsec);
	    unsigned int cellPrevRadial = cell - Qbase->Nsec;
	    if (VRadial->Field[cell] > 0.0)
		QStar->Field[cell] = Qbase->Field[cellPrevRadial] +
				     (Rmed[nRadial] - Rmed[nRadial - 1] -
				      VRadial->Field[cell] * dt) *
					 0.5 * dq[cellPrevRadial];
	    else
		QStar->Field[cell] =
		    Qbase->Field[cell] - (Rmed[nRadial + 1] - Rmed[nRadial] +
					  VRadial->Field[cell] * dt) *
					     0.5 * dq[cell];
	}
    }
    // TODO: check here
    for (nAzimuthal = 0; nAzimuthal < Qbase->Nsec; ++nAzimuthal) {
	(*QStar)(0, nAzimuthal) = 0.0;
	if ((parameters::boundary_outer !=
	     parameters::boundary_condition_evanescent) &&
	    (parameters::boundary_outer !=
	     parameters::boundary_condition_boundary_layer) &&
	    (parameters::boundary_outer !=
	     parameters::boundary_condition_precribed_time_variable)) {
	    (*QStar)(QStar->get_max_radial(), nAzimuthal) = 0.0;
	}
    }
}

/**
 * @brief ComputeStarTheta update the quantity Qbase with advection in Azimuthal
 * direction, results are written to Qstar
 * @param Qbase
 * @param VAzimuthal
 * @param QStar
 * @param dt
 */
void ComputeStarTheta(PolarGrid *Qbase, PolarGrid *VAzimuthal, PolarGrid *QStar,
		      double dt)
{
    unsigned int nRadial, nAzimuthal, cell;

    int ljp, ljm, jm;
    double dqp, dqm, dxtheta, ksi, invdxtheta;

    for (nRadial = 0; nRadial < Qbase->Nrad; ++nRadial) {
	dxtheta = dphi * Rmed[nRadial];
	invdxtheta = 1.0 / dxtheta;
	for (nAzimuthal = 0; nAzimuthal < Qbase->Nsec; ++nAzimuthal) {
	    cell = nAzimuthal + nRadial * Qbase->Nsec;

	    ljp = cell + 1;
	    ljm = cell - 1;
	    if (nAzimuthal == 0) {
		ljm = nRadial * Qbase->Nsec + Qbase->Nsec - 1;
	    }
	    if (nAzimuthal == Qbase->Nsec - 1) {
		ljp = nRadial * Qbase->Nsec;
	    }
	    dqm = (Qbase->Field[cell] - Qbase->Field[ljm]);
		dqp = (Qbase->Field[ljp] - Qbase->Field[cell]);
		dq[cell] = 0.5*flux_limiter(dqp, dqm) * invdxtheta;
	}
	for (nAzimuthal = 0; nAzimuthal < Qbase->Nsec; ++nAzimuthal) {
	    cell = cell_number(nRadial, nAzimuthal, Qbase->Nsec);

	    jm = nAzimuthal - 1;
	    if (nAzimuthal == 0)
		jm = Qbase->Nsec - 1;
	    ljm = jm + nRadial * Qbase->Nsec;
	    ksi = VAzimuthal->Field[cell] * dt;
	    if (ksi > 0.0) {
		QStar->Field[cell] =
		    Qbase->Field[ljm] + (dxtheta - ksi) * dq[ljm];
	    } else {
		QStar->Field[cell] =
		    Qbase->Field[cell] - (dxtheta + ksi) * dq[cell];
	    }
	}
    }
}

/**
	Calculates radial/angular momenta from radial/azimuthal velocities
*/
void compute_momenta_from_velocities(t_polargrid &density,
				     t_polargrid &v_radial,
				     t_polargrid &v_azimuthal)
{
    for (unsigned int n_radial = 0; n_radial < density.get_size_radial();
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < density.get_size_azimuthal(); ++n_azimuthal) {
	    radial_momentum_plus(n_radial, n_azimuthal) =
		density(n_radial, n_azimuthal) *
		v_radial(n_radial + 1, n_azimuthal);
	    radial_momentum_minus(n_radial, n_azimuthal) =
		density(n_radial, n_azimuthal) *
		v_radial(n_radial, n_azimuthal);
	    angular_momentum_plus(n_radial, n_azimuthal) =
		density(n_radial, n_azimuthal) *
		(v_azimuthal(n_radial,
			     n_azimuthal == v_azimuthal.get_max_azimuthal()
				 ? 0
				 : n_azimuthal + 1) +
		 Rb[n_radial] * refframe::OmegaFrame) *
		Rb[n_radial];
	    angular_momentum_minus(n_radial, n_azimuthal) =
		density(n_radial, n_azimuthal) *
		(v_azimuthal(n_radial, n_azimuthal) +
		 Rb[n_radial] * refframe::OmegaFrame) *
		Rb[n_radial];
	}
    }
}

/**
	Calculates radial/azimuthal velocities from radial/angular momenta
*/
void compute_velocities_from_momenta(t_polargrid &density,
				     t_polargrid &v_radial,
				     t_polargrid &v_azimuthal)
{
    for (unsigned int n_radial = 0; n_radial < density.get_size_radial();
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < density.get_size_azimuthal(); ++n_azimuthal) {
	    unsigned int n_azimuthal_minus =
		(n_azimuthal == 0 ? density.get_max_azimuthal()
				  : n_azimuthal - 1);

	    // radial velocity from radial momentum
	    if (n_radial == 0) {
		v_radial(n_radial, n_azimuthal) = 0.0;
	    } else {
		v_radial(n_radial, n_azimuthal) =
		    (radial_momentum_plus(n_radial - 1, n_azimuthal) +
		     radial_momentum_minus(n_radial, n_azimuthal)) /
		    (density(n_radial - 1, n_azimuthal) +
		     density(n_radial, n_azimuthal));
	    }

	    // azimuthal velocity from angular momentum
	    v_azimuthal(n_radial, n_azimuthal) =
		(angular_momentum_plus(n_radial, n_azimuthal_minus) +
		 angular_momentum_minus(n_radial, n_azimuthal)) /
		    (density(n_radial, n_azimuthal_minus) +
		     density(n_radial, n_azimuthal)) *
		    InvRb[n_radial] -
		Rb[n_radial] * refframe::OmegaFrame;
	}
    }
}

/**
 * @brief VanLeerRadial perform advection of quantity Qbase via transforming it
 * into specific quantity
 * @param data
 * @param VRadial
 * @param Qbase
 * @param dt
 */
void VanLeerRadial(t_data &data, PolarGrid *VRadial, PolarGrid *Qbase,
		   double dt)
{
    unsigned int nRadial, nAzimuthal;
    bool is_density = false;
    units::t_unit *unit = Qbase->get_unit();
    if (unit != nullptr) {
	is_density = strcmp(unit->get_cgs_symbol(),
			    units::surface_density.get_cgs_symbol()) == 0;
    }

    divise_polargrid(*Qbase, data[t_data::DENSITY_INT],
		     *Work); // work = qbase/densityint
    compute_star_radial(Work, VRadial, QRStar, dt);

    for (nRadial = 0; nRadial < Qbase->get_size_radial(); ++nRadial) {
	for (nAzimuthal = 0; nAzimuthal < Qbase->get_size_azimuthal();
	     ++nAzimuthal) {
	    const unsigned int cell = nAzimuthal + nRadial * Qbase->Nsec;
	    const int lip = cell + Qbase->Nsec;
	    // mass crossing inf interface
	    const double varq_inf =
		dt * dphi * Rinf[nRadial] * QRStar->Field[cell] *
		DensityStar->Field[cell] * VRadial->Field[cell];
	    // mass crossing sup interface
	    const double varq_sup =
		dt * dphi * Rsup[nRadial] * QRStar->Field[lip] *
		DensityStar->Field[lip] * VRadial->Field[lip];
	    // update density
	    Qbase->Field[cell] += (varq_inf - varq_sup) * InvSurf[nRadial];

	    if (is_density) {
		if (parameters::write_disk_quantities) {
		    // TODO: boundary
		    // if ((nRadial == 0) && (parameters::boundary_inner ==
		    // parameters::boundary_condition_open)) if ((nRadial == 0)
		    // && (OpenInner))
		    if (CPU_Rank == 0 && nRadial == 1) {
			if (varq_inf > 0) {
			    sum_without_ghost_cells(MassDelta.InnerPositive,
						    varq_inf, nRadial);
			} else {
			    sum_without_ghost_cells(MassDelta.InnerNegative,
						    varq_inf, nRadial);
			}
		    } else if (CPU_Rank == CPU_Highest &&
			       nRadial == Qbase->get_max_radial() - 1) {
			if (varq_sup > 0) {
			    sum_without_ghost_cells(MassDelta.OuterNegative,
							-varq_sup, nRadial);
			} else {
			    sum_without_ghost_cells(MassDelta.OuterPositive,
							-varq_sup, nRadial);
			}
		    }
		}
		if (parameters::write_massflow) {
		    data[t_data::MASSFLOW](nRadial, nAzimuthal) += varq_inf;
		    if (CPU_Rank == CPU_Highest &&
			nRadial == Qbase->get_max_radial()) {
			data[t_data::MASSFLOW](nRadial, nAzimuthal) += varq_sup;
		    }
		}
	    }
	}
    }
}

/**
 * @brief VanLeerTheta perform advection on quantity Qbase by transforming it
 * into specific quantity
 * @param data
 * @param VAzimuthal
 * @param Qbase
 * @param dt
 */
void VanLeerTheta(t_data &data, PolarGrid *VAzimuthal, PolarGrid *Qbase,
		  double dt)
{
    unsigned int nRadial, nAzimuthal, cell;

    int ljp;

    divise_polargrid(*Qbase, data[t_data::DENSITY_INT],
		     *Work); // work = qbase/densityint

    ComputeStarTheta(Work, VAzimuthal, QRStar, dt);

    for (nRadial = 0; nRadial < Qbase->Nrad; ++nRadial) {
	double dxrad = (Rsup[nRadial] - Rinf[nRadial]) * dt;
	double invsurf = 1.0 / Surf[nRadial];

	if ((UniformTransport == NO) || (NoSplitAdvection[nRadial] == NO)) {
	    for (nAzimuthal = 0; nAzimuthal < Qbase->Nsec; ++nAzimuthal) {
		double varq;

		cell = nAzimuthal + nRadial * Qbase->Nsec;
		ljp = cell + 1;

		if (nAzimuthal == Qbase->Nsec - 1)
		    ljp = nRadial * Qbase->Nsec;

		varq = dxrad * QRStar->Field[cell] * DensityStar->Field[cell] *
		       VAzimuthal->Field[cell];
		varq -= dxrad * QRStar->Field[ljp] * DensityStar->Field[ljp] *
			VAzimuthal->Field[ljp];

		Qbase->Field[cell] += varq * invsurf;
	    }
	}
    }
}

/**
	boundary_layer_mass_influx modifies outermost value of DensityStar to
   account for a constant mass accretion rate M_dot = -2\pi\Sigma v_rad r
*/

void boundary_layer_mass_influx(PolarGrid *QStar, PolarGrid *VRadial)
{

    // printf("QStar->get_max_radial() = %d, VRadial->get_max_radial() = %d\n",
    // QStar->get_max_radial(), VRadial->get_max_radial());

    for (unsigned int n_azimuthal = 0;
	 n_azimuthal < QStar->get_size_azimuthal(); ++n_azimuthal) {
	(*QStar)(QStar->get_max_radial() - 1, n_azimuthal) =
	    -parameters::mass_accretion_rate /
	    ((*VRadial)(VRadial->get_max_radial() - 1, n_azimuthal) *
	     NAzimuthal);
    }
}
