	/**
	\file TransportEuler.c

	Functions that handle the transport substep of a hydrodynamical time step. The FARGO algorithm is implemented here. The transport is performed in a manner similar to what is done for the ZEUS code (Stone & Norman, 1992), except for the momenta transport (we define a left and right momentum for each zone, which we declare zone centered; we then transport then normally, and deduce the new velocity in each zone by a proper averaging).
*/

#include <math.h>
#include <string.h>

#include "polargrid.h"
#include "radialgrid.h"

#include "TransportEuler.h"
#include "SideEuler.h"
#include "constants.h"
#include "macros.h"
#include "LowTasks.h"
#include "global.h"
#include "SourceEuler.h"
#include "parameters.h"

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

extern double OmegaFrame;
//static double VMed[MAX1D];
static int Nshift[MAX1D];
static boolean NoSplitAdvection[MAX1D];
static boolean UniformTransport;
extern int TimeStep;
extern boolean FastTransport;
extern bool Adiabatic;
extern boolean OpenInner;

/**
	Initializes (allocates) all variables needed.
*/
void InitTransport()
{
	radial_momentum_plus.set_scalar(true);
	radial_momentum_plus.set_size(NRadial, NAzimuthal);
	radial_momentum_minus.set_scalar(true);
	radial_momentum_minus.set_size(NRadial, NAzimuthal);
	angular_momentum_plus.set_scalar(true);
	angular_momentum_plus.set_size(NRadial, NAzimuthal);
	angular_momentum_minus.set_scalar(true);
	angular_momentum_minus.set_size(NRadial, NAzimuthal);

	v_azimuthal_res.set_vector(false);
	v_azimuthal_res.set_size(NRadial, NAzimuthal);

	v_azimuthal_mean.set_vector(false);
	v_azimuthal_mean.set_size(NRadial);

	Work = CreatePolarGrid(NRadial, NAzimuthal, "WorkGrid");

	QRStar = CreatePolarGrid(NRadial, NAzimuthal, "QRStar");
	// TODO: This has to be checked!!!
	QRStar->set_vector(true);
	QRStar->set_size(NRadial, NAzimuthal);

	DensityStar = CreatePolarGrid(NRadial, NAzimuthal, "DensityStar");
	// TODO: This has to be checked! set_vector(true) is wrong, since there cannot be flow through outermost ghost cell
	DensityStar->set_vector(true);
	DensityStar->set_size(NRadial, NAzimuthal);

	TempShift = (double*)malloc(NRadial*NAzimuthal*sizeof(double));

	dq = (double*)malloc(NRadial*NAzimuthal*sizeof(double));
}

/**
	Frees all variables allocate before by InitTransport.
*/
void FreeTransport()
{
	delete Work;
	delete QRStar;

	free(TempShift);
	free(dq);
}

void Transport(t_data &data, PolarGrid* Density, PolarGrid* VRadial, PolarGrid* VAzimuthal, PolarGrid* Energy, double dt)
{
	compute_momenta_from_velocities(*Density, *VRadial, *VAzimuthal);

	/* No-Alternate Directionnal Splitting */
	OneWindRad(data, Density, VRadial, Energy, dt);
	OneWindTheta(data, Density, VAzimuthal, Energy, dt);

	compute_velocities_from_momenta(*Density, *VRadial, *VAzimuthal);
}

void OneWindRad(t_data &data, PolarGrid* Density, PolarGrid* VRadial, PolarGrid* Energy, double dt)
{
	compute_star_radial(Density, VRadial, DensityStar, dt);

	// boundary layer:
	if( parameters::boundary_outer == parameters::boundary_condition_boundary_layer && CPU_Rank == CPU_Highest )
		boundary_layer_mass_influx(DensityStar, VRadial);

 	copy_polargrid(data[t_data::DENSITY_INT], *Density);

	VanLeerRadial(data, VRadial, &radial_momentum_plus, dt);
	VanLeerRadial(data, VRadial, &radial_momentum_minus, dt);
	VanLeerRadial(data, VRadial, &angular_momentum_plus, dt);
	VanLeerRadial(data, VRadial, &angular_momentum_minus, dt);

	if (Adiabatic)
		VanLeerRadial(data, VRadial, Energy, dt);

	LostMass += VanLeerRadial(data, VRadial, Density, dt); /* MUST be the last line */
}

/* Hereafter are the new specific procedures to the fast algorithm */

/**
	Compute average azimuthal velocities
*/
void compute_average_azimuthal_velocity(t_polargrid &v_azimuthal, double dt)
{
	for (unsigned int n_radial = 0; n_radial <= v_azimuthal.get_max_radial(); ++n_radial) {
		double v_azimuthal_sum = 0.0;

		for (unsigned int n_azimuthal = 0; n_azimuthal <= v_azimuthal.get_max_azimuthal(); ++n_azimuthal) {
			v_azimuthal_sum += v_azimuthal(n_radial, n_azimuthal);
		}

		v_azimuthal_mean(n_radial) = v_azimuthal_sum/(double)v_azimuthal.get_size_azimuthal();
	}
}

/**
	Compute (azimuthal) residual velocities
*/
void compute_residual_velocity(t_polargrid &v_azimuthal, double dt)
{
	for (unsigned int n_radial = 0; n_radial <= v_azimuthal.get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= v_azimuthal.get_max_azimuthal(); ++n_azimuthal) {
			v_azimuthal_res(n_radial, n_azimuthal) = v_azimuthal(n_radial, n_azimuthal)-v_azimuthal_mean(n_radial);
		}
	}
}

void ComputeConstantResidual(PolarGrid* VAzimuthal, double dt)
{
	int i,j,l,nr,ns;
	long nitemp;
	double *vt, *vres, Ntilde, Nround, maxfrac, invdt, dpinvns;
	nr = VAzimuthal->Nrad;
	ns = VAzimuthal->Nsec;
	vt = VAzimuthal->Field;
	vres = v_azimuthal_res.Field;
	invdt = 1.0/dt;
	dpinvns=2.0*PI/(double)ns;

	if (FastTransport)
		maxfrac = 1.0;		/* Fast algorithm */
	else
		maxfrac = 0.0;

	for (i = 0; i < nr; i++) {
		Ntilde = v_azimuthal_mean(i)*InvRmed[i]*dt*(double)ns/2.0/PI;
		Nround = floor(Ntilde+0.5);
		nitemp = (long)Nround;
		Nshift[i] = (long)nitemp;
		for (j = 0; j < ns; j++) {
			l=j+i*ns;
			vt[l] = (Ntilde-Nround)*Rmed[i]*invdt*dpinvns;
		}
		if (maxfrac < 0.5) {
			NoSplitAdvection[i] = YES;
			for (j = 0; j < ns; j++) {
				l=j+i*ns;
				vres[l] = vt[l]+vres[l];
				vt[l] = 0.0;
			}
		} else {
			NoSplitAdvection[i] = NO;
		}
	}
}

void AdvectSHIFT(t_polargrid &array)
{
	int i,j,ji,l,li,nr,ns;
	double *val;
	val = array.Field;
	nr  = array.Nrad;
	ns  = array.Nsec;

	for (i = 0; i < nr; i++) {
		for (j = 0; j < ns; j++) {
			ji = j-Nshift[i];
			while (ji < 0)
				ji += ns;
			while (ji >= ns)
				ji -= ns;
			l = j+i*ns;
			li= ji+i*ns;
			TempShift[l]=val[li];
		}
		for (j = 0; j < ns; j++) {
			l = j+i*ns;
			val[l] = TempShift[l];
		}
	}
}

void OneWindTheta(t_data &data, PolarGrid* Density, PolarGrid* VAzimuthal, PolarGrid* Energy, double dt)
{
	compute_average_azimuthal_velocity(*VAzimuthal, dt);
	compute_residual_velocity(*VAzimuthal, dt);
	ComputeConstantResidual(VAzimuthal, dt);	/* Constant residual is in VAzimuthal from now on */
	UniformTransport = NO;
	QuantitiesAdvection(data, Density, &v_azimuthal_res, Energy, dt);
	UniformTransport = YES;
	QuantitiesAdvection(data, Density, VAzimuthal, Energy, dt);
	AdvectSHIFT(radial_momentum_plus);
	AdvectSHIFT(radial_momentum_minus);
	AdvectSHIFT(angular_momentum_plus);
	AdvectSHIFT(angular_momentum_minus);
	if (Adiabatic)
		AdvectSHIFT(*Energy);
	AdvectSHIFT(*Density);
}

/* End of new specific procedures to the fast algorithm */

void QuantitiesAdvection(t_data &data, PolarGrid* Density, PolarGrid* VAzimuthal, PolarGrid* Energy, double dt)
{
	ComputeStarTheta(Density, VAzimuthal, DensityStar, dt);
	copy_polargrid(data[t_data::DENSITY_INT], *Density);
	VanLeerTheta(data, VAzimuthal, &radial_momentum_plus, dt);
	VanLeerTheta(data, VAzimuthal, &radial_momentum_minus, dt);
	VanLeerTheta(data, VAzimuthal, &angular_momentum_plus, dt);
	VanLeerTheta(data, VAzimuthal, &angular_momentum_minus, dt);
	if (Adiabatic)
		VanLeerTheta(data, VAzimuthal, Energy, dt);
	VanLeerTheta(data, VAzimuthal, Density, dt); /* MUST be the last line */
}

/**

*/
//void compute_star_radial(t_polargrid* base, t_polargrid* V_Radial, t_polargrid* star, double dt)
void compute_star_radial(t_polargrid* Qbase, t_polargrid* VRadial, t_polargrid* QStar, double dt)
{
/*	// calculate all dq
	for (unsigned int n_radial = 0; n_radial <= base.get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= base.get_max_azimuthal(); ++n_azimuthal) {
			// inner and outer most cells need dq = 0
			if ((n_radial == 0) || (n_radial == base.get_max_radial())) {
				dq[n_radial] = 0.0;
			} else {
				double dq_minus = base(n_radial, n_azimuthal)-base(n_azimuthal-1,n_azimuthal)*InvDiffRmed[nRadial];
				double dq_plus = base(n_radial+1, n_azimuthal)-base(n_radial,n_azimuthal)*InvDiffRmed[nRadial+1];

				if (dq_minus * dq_plus > 0.0) {
					dq[n_radial] = 2.0*dq_minus*dq_plus/(dq_minus+dq_plus);
				} else {
					dq[n_radial] = 0.0;
				}
			}
		}
	}

	// calculate star
	for (unsigned int n_radial = 0; n_radial <= base.get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= base.get_max_azimuthal(); ++n_azimuthal) {
			if (v_radial(n_radial, n_azimuthal) > 0.0) {
				star(n_radial, n_azimuthal) = base(n_radial-1, n_azimuthal)+(Rmed[nRadial]-Rmed[nRadial-1]-v_radial(n_radial,n_azimuthal)*dt)*0.5*dq[nRadial-1];
			} else {
				star(n_radial, n_azimuthal) = base(n_radial, n_azimuthal)-(Rmed[nRadial+1]-Rmed[nRadial]+v_radial(n_radial,n_azimuthal)*dt)*0.5*dq[nRadial];
			}
		}
	}*/

	unsigned int nRadial, nAzimuthal;
	double dqp, dqm;

	for (nAzimuthal = 0; nAzimuthal < Qbase->Nsec; ++nAzimuthal) {

		for (nRadial = 0; nRadial < Qbase->Nrad; ++nRadial) {
			unsigned int cell = CELL(nRadial,nAzimuthal,Qbase->Nsec);
			unsigned int cellPrevRadial = cell-Qbase->Nsec;
			unsigned int cellNextRadial = cell+Qbase->Nsec;

			if ((nRadial == 0) || (nRadial == Qbase->Nrad-1)) {
				dq[nRadial] = 0.0;
			} else {
				dqm = (Qbase->Field[cell]-Qbase->Field[cellPrevRadial])*InvDiffRmed[nRadial];
				dqp = (Qbase->Field[cellNextRadial]-Qbase->Field[cell])*InvDiffRmed[nRadial+1];
				if (dqp * dqm > 0.0)
					dq[nRadial] = 2.0*dqp*dqm/(dqp+dqm);
				else
					dq[nRadial] = 0.0;
			}
		}

		// TODO: changed to nRadial =1 because of nRadial-1
		// TODO: potential problem: Using Rmed[nRadial] - Rmed[nRadial-1] for a-mesh Qties (v_rad,..) as well as for b-mesh Qties (Density,...)
		for (nRadial = 1; nRadial < Qbase->Nrad; ++nRadial) {
			unsigned int cell = CELL(nRadial, nAzimuthal, Qbase->Nsec);
			unsigned int cellPrevRadial = cell-Qbase->Nsec;

			if (VRadial->Field[cell] > 0.0)
				QStar->Field[cell] = Qbase->Field[cellPrevRadial]+(Rmed[nRadial]-Rmed[nRadial-1]-VRadial->Field[cell]*dt)*0.5*dq[nRadial-1];
			else
				QStar->Field[cell] = Qbase->Field[cell]-(Rmed[nRadial+1]-Rmed[nRadial]+VRadial->Field[cell]*dt)*0.5*dq[nRadial];
		}
		// TODO: check here
		(*QStar)(0,nAzimuthal) = 0.0;
		(*QStar)(QStar->get_max_radial(),nAzimuthal) = 0.0;
	}
}

void ComputeStarTheta(PolarGrid* Qbase, PolarGrid* VAzimuthal, PolarGrid* QStar, double dt)
{
	unsigned int nRadial, nAzimuthal,cell;

	int ljp,ljm,jm;
	double dqp, dqm,dxtheta,ksi,invdxtheta;

	for (nRadial = 0; nRadial < Qbase->Nrad; ++nRadial) {
		dxtheta = 2.0*PI/(double)Qbase->Nsec*Rmed[nRadial];
		invdxtheta = 1.0/dxtheta;
		for (nAzimuthal = 0; nAzimuthal < Qbase->Nsec; ++nAzimuthal) {
			cell = nAzimuthal+nRadial*Qbase->Nsec;

			ljp = cell+1;
			ljm = cell-1;
			if (nAzimuthal == 0)
				ljm = nRadial*Qbase->Nsec+Qbase->Nsec-1;
			if (nAzimuthal == Qbase->Nsec-1)
				ljp = nRadial*Qbase->Nsec;
			dqm = (Qbase->Field[cell]-Qbase->Field[ljm]);
			dqp = (Qbase->Field[ljp]-Qbase->Field[cell]);
			if (dqp * dqm > 0.0)
				dq[cell] = dqp*dqm/(dqp+dqm)*invdxtheta;
			else
				dq[cell] = 0.0;
		}
		for (nAzimuthal = 0; nAzimuthal < Qbase->Nsec; ++nAzimuthal) {
			//cell = nAzimuthal+nRadial*Qbase->Nsec;
			cell = CELL(nRadial, nAzimuthal, Qbase->Nsec);

			jm = nAzimuthal-1;
			if (nAzimuthal == 0)
				jm = Qbase->Nsec-1;
			ljm = jm+nRadial*Qbase->Nsec;
			ksi=VAzimuthal->Field[cell]*dt;
			if (ksi > 0.0)
				QStar->Field[cell] = Qbase->Field[ljm]+(dxtheta-ksi)*dq[ljm];
			else
				QStar->Field[cell] = Qbase->Field[cell]-(dxtheta+ksi)*dq[cell];
		}
	}
}

/**
	Calculates radial/angular momenta from radial/azimuthal velocities
*/
void compute_momenta_from_velocities(t_polargrid &density, t_polargrid &v_radial, t_polargrid &v_azimuthal)
{
	for (unsigned int n_radial = 0; n_radial <= density.get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= density.get_max_azimuthal(); ++n_azimuthal) {
			radial_momentum_plus(n_radial, n_azimuthal) = density(n_radial, n_azimuthal)*v_radial(n_radial+1, n_azimuthal);
			radial_momentum_minus(n_radial, n_azimuthal) = density(n_radial, n_azimuthal)*v_radial(n_radial, n_azimuthal);
			angular_momentum_plus(n_radial, n_azimuthal) = density(n_radial, n_azimuthal)*(v_azimuthal(n_radial, n_azimuthal == v_azimuthal.get_max_azimuthal() ? 0 : n_azimuthal+1)+Rb[n_radial]*OmegaFrame)*Rb[n_radial];
			angular_momentum_minus(n_radial, n_azimuthal) = density(n_radial, n_azimuthal)*(v_azimuthal(n_radial, n_azimuthal)+Rb[n_radial]*OmegaFrame)*Rb[n_radial];
		}
	}
}

/**
	Calculates radial/azimuthal velocities from radial/angular momenta
*/
void compute_velocities_from_momenta(t_polargrid &density, t_polargrid &v_radial, t_polargrid &v_azimuthal)
{
	for (unsigned int n_radial = 0; n_radial <= density.get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= density.get_max_azimuthal(); ++n_azimuthal) {
			unsigned int n_azimuthal_minus = (n_azimuthal == 0 ? density.get_max_azimuthal() : n_azimuthal-1);

			// radial velocity from radial momentum
			if (n_radial == 0) {
				v_radial(n_radial, n_azimuthal) = 0.0;
			} else {
				v_radial(n_radial, n_azimuthal) = (radial_momentum_plus(n_radial-1, n_azimuthal)+radial_momentum_minus(n_radial, n_azimuthal))/(density(n_radial-1, n_azimuthal)+density(n_radial, n_azimuthal));
			}

			// azimuthal velocity from angular momentum
			v_azimuthal(n_radial, n_azimuthal) = (angular_momentum_plus(n_radial, n_azimuthal_minus)+angular_momentum_minus(n_radial, n_azimuthal))/(density(n_radial,n_azimuthal_minus)+density(n_radial,n_azimuthal))*InvRb[n_radial]-Rb[n_radial]*OmegaFrame;
		}
	}
}

double VanLeerRadial(t_data &data, PolarGrid* VRadial, PolarGrid* Qbase, double dt)
{
	unsigned int nRadial, nAzimuthal, cell;
	int lip;
	double dtheta;
	double LostByDisk=0.0;

	divise_polargrid(*Qbase, data[t_data::DENSITY_INT], *Work); // work = qbase/densityint
	compute_star_radial(Work, VRadial, QRStar, dt);

	dtheta = 2.0*PI/(double)Qbase->Nsec;

	for (nRadial = 0; nRadial <= Qbase->get_max_radial(); ++nRadial) {
		for (nAzimuthal = 0; nAzimuthal <= Qbase->get_max_azimuthal(); ++nAzimuthal) {
			double varq;

			cell=nAzimuthal+nRadial*Qbase->Nsec;
			lip=cell+Qbase->Nsec;
			varq =dt*dtheta*Rinf[nRadial]*QRStar->Field[cell]*DensityStar->Field[cell]*VRadial->Field[cell];
			varq-=dt*dtheta*Rsup[nRadial]*QRStar->Field[lip]*DensityStar->Field[lip]*VRadial->Field[lip];
			Qbase->Field[cell] += varq*InvSurf[nRadial];

			// TODO: boundary
			if ((nRadial == 0) && (parameters::boundary_inner == parameters::boundary_condition_open))
			//if ((nRadial == 0) && (OpenInner))
				LostByDisk += varq;
		}
	}

	return LostByDisk;
}

void VanLeerTheta(t_data &data, PolarGrid* VAzimuthal, PolarGrid* Qbase, double dt)
{
	unsigned int nRadial, nAzimuthal, cell;

	int ljp;

	divise_polargrid(*Qbase, data[t_data::DENSITY_INT], *Work);

	ComputeStarTheta(Work, VAzimuthal, QRStar, dt);

	for (nRadial = 0; nRadial < Qbase->Nrad; ++nRadial) {
		double dxrad = (Rsup[nRadial]-Rinf[nRadial])*dt;
		double invsurf = 1.0/Surf[nRadial];

		if ((UniformTransport == NO) || (NoSplitAdvection[nRadial] == NO)) {
			for (nAzimuthal = 0; nAzimuthal < Qbase->Nsec; ++nAzimuthal) {
				double varq;

				cell=nAzimuthal+nRadial*Qbase->Nsec;
				ljp=cell+1;

				if (nAzimuthal == Qbase->Nsec-1)
					ljp=nRadial*Qbase->Nsec;

				varq  = dxrad*QRStar->Field[cell]*DensityStar->Field[cell]*VAzimuthal->Field[cell];
				varq -= dxrad*QRStar->Field[ljp]*DensityStar->Field[ljp]*VAzimuthal->Field[ljp];

				Qbase->Field[cell] += varq*invsurf;
			}
		}
	}

}

/**
	boundary_layer_mass_influx modifies outermost value of DensityStar to account for a constant mass accretion rate M_dot = -2\pi\Sigma v_rad r
*/

void boundary_layer_mass_influx(PolarGrid* QStar, PolarGrid* VRadial) {

	//printf("QStar->get_max_radial() = %d, VRadial->get_max_radial() = %d\n", QStar->get_max_radial(), VRadial->get_max_radial());

	for (unsigned int n_azimuthal = 0; n_azimuthal <= QStar->get_max_azimuthal(); ++n_azimuthal) {
		(*QStar)(QStar->get_max_radial()-1, n_azimuthal) = - parameters::mass_accretion_rate / ( (*VRadial)(VRadial->get_max_radial()-1, n_azimuthal) * NAzimuthal);
	}
}