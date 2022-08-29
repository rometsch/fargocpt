/*
	\file SideEuler.c

	Total mass and angular momentum monitoring, and boundary conditions. In
   addition, this file contains a few low-level functions that manipulate
   PolarGrid 's or initialize the forces evaluation.

*/

#include <float.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>

#include "Force.h"
#include "Interpret.h"
#include "LowTasks.h"
#include "SideEuler.h"
#include "SourceEuler.h"
#include "Theo.h"
#include "axilib.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "macros.h"
#include "nongnu.h"
#include "parameters.h"
#include "quantities.h"
#include "selfgravity.h"
#include "util.h"
#include "viscosity.h"
#include <cassert>

extern boolean Damping;
extern boolean OuterSourceMass;


/**
	Divide each cell of a polar grid with a cell of another polargrid

	\param Num numerator
	\param Denom denominator
	\param Result result
*/
void divise_polargrid(t_polargrid &num, t_polargrid &denom, t_polargrid &result)
{
    const unsigned int Nmax =
	result.get_size_radial() * result.get_size_azimuthal();
    for (unsigned int n = 0; n < Nmax; n++) {
	assert(denom.Field[n] > 0.0);
	/// denom + DBL_EPSILON can cause problems because DBL_EPSILON can be
	/// bigger than denom, depending on units.
	result.Field[n] =
	    num.Field[n] /
	    denom.Field[n]; /// in case of crash, use something like
			    /// (denom.Field[n] + 1.0e-200) instead.
    }
}

/**

*/
void InitCellCenterCoordinates()
{
    unsigned int nRadial, nAzimuthal, cell;

    delete CellCenterY;
    delete CellCenterX;

    CellCenterX = CreatePolarGrid(NRadial + 1, NAzimuthal, "cell_center_x");
    CellCenterY = CreatePolarGrid(NRadial + 1, NAzimuthal, "cell_center_y");

    for (nRadial = 0; nRadial < CellCenterX->Nrad; ++nRadial) {
	for (nAzimuthal = 0; nAzimuthal < CellCenterX->Nsec; ++nAzimuthal) {
	    cell = nAzimuthal + nRadial * CellCenterX->Nsec;
	    CellCenterX->Field[cell] =
		Rmed[nRadial] * std::cos(dphi * (double)nAzimuthal);
	    CellCenterY->Field[cell] =
		Rmed[nRadial] * std::sin(dphi * (double)nAzimuthal);
	}
    }
}

void FreeCellCenterCoordinates()
{
    delete CellCenterY;
    delete CellCenterX;
}

/**
	\param VRadial radial velocity polar grid
	\param Density density polar grid
	\param Energy energy polar grid
*/
void NonReflectingBoundary_inner(t_data &data, t_polargrid *VRadial,
				 t_polargrid *Density, t_polargrid *Energy)
{
    unsigned int nRadial, nAzimuthal, cell;

    unsigned int jp;
    int lp, i_angle;
    double *rho, *vr;
    double dangle, mean, vr_med;
    double cs0, cs1;

    rho = Density->Field;
    vr = VRadial->Field;

    if (CPU_Rank == 0) {
	cs0 = 0.0;
	cs1 = 0.0;
	for (nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
	    cs0 += data[t_data::SOUNDSPEED].Field[nAzimuthal];
	    cs1 += data[t_data::SOUNDSPEED].Field[Density->Nsec + nAzimuthal];
	}
	cs0 /= (double)(Density->Nsec);
	cs1 /= (double)(Density->Nsec);
	nRadial = 1;
	// The expression below should be refined
	// We need to know the orbital frequency of the nearest planet
	dangle = (pow(Rinf[1], -1.5) - 1.0) / (.5 * (cs0 + cs1));
	dangle *= (Rmed[1] - Rmed[0]);
	i_angle = (int)(dangle / 2.0 / M_PI * (double)NAzimuthal + .5);

	for (nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
	    cell = nAzimuthal + nRadial * Density->Nsec;
	    jp = nAzimuthal + i_angle;
	    if (jp >= Density->Nsec)
		jp -= Density->Nsec;
	    /*if (jp < 0)
		    jp += Density->Nsec;*/
	    lp = jp;
	    rho[lp] = rho[cell]; // copy first ring into ghost ring
	    Energy->Field[lp] =
		Energy->Field[cell]; // copy first ring into ghost ring
	    vr_med = -data[t_data::SOUNDSPEED].Field[cell] *
		     (rho[cell] - SigmaMed[1]) / SigmaMed[1];
	    vr[cell] = 2. * vr_med - vr[cell + Density->Nsec];
	}
	mean = 0.0;
	for (nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
	    mean += rho[nAzimuthal];
	}
	mean /= (double)(Density->Nsec);
	for (nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
	    rho[nAzimuthal] += SigmaMed[0] - mean;
	}
	mean = 0.0;
	for (nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
	    mean += Energy->Field[nAzimuthal];
	}
	mean /= (double)Density->Nsec;
	for (nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
	    Energy->Field[nAzimuthal] += EnergyMed[0] - mean;
	}
    }
}

/**
	\param VRadial radial velocity polar grid
	\param Density density polar grid
	\param Energy energy polar grid
 */
void NonReflectingBoundary_outer(t_data &data, t_polargrid *VRadial,
				 t_polargrid *Density, t_polargrid *Energy)
{
    unsigned int nRadial, nAzimuthal, cell;

    unsigned int jp;
    int lp, i_angle;
    double *rho, *vr;
    double dangle, mean, vr_med;
    double csnrm1, csnrm2;

    rho = Density->Field;
    vr = VRadial->Field;

    if (CPU_Rank == CPU_Highest) {
	csnrm2 = 0.0;
	csnrm1 = 0.0;
	for (nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
	    csnrm2 +=
		data[t_data::SOUNDSPEED]
		    .Field[(Density->Nrad - 2) * Density->Nsec + nAzimuthal];
	    csnrm1 +=
		data[t_data::SOUNDSPEED]
		    .Field[(Density->Nrad - 1) * Density->Nsec + nAzimuthal];
	}
	csnrm1 /= (double)(Density->Nsec);
	csnrm2 /= (double)(Density->Nsec);
	nRadial = Density->Nrad - 1; // The expression below should be refined
	// We need to know the orbital frequency of the nearest planet
	dangle = (pow(Rinf[Density->Nrad - 2], -1.5) - 1.0) /
		 (.5 * (csnrm1 + csnrm2));
	dangle *= (Rmed[Density->Nrad - 1] - Rmed[Density->Nrad - 2]);
	i_angle = (int)(dangle / 2.0 / M_PI * (double)NAzimuthal + .5);

	for (nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
	    cell = nAzimuthal + nRadial * Density->Nsec;
	    jp = nAzimuthal - i_angle;
	    if (jp >= Density->Nsec)
		jp -= Density->Nsec;
	    /*if (jp < 0)
		    jp += Density->Nsec;*/
	    lp = jp + (nRadial - 1) * Density->Nsec;
	    rho[cell] = rho[lp]; // copy first ring into ghost ring
	    Energy->Field[cell] =
		Energy->Field[lp]; // copy first ring into ghost ring
	    vr_med = data[t_data::SOUNDSPEED].Field[cell] *
		     (rho[cell - Density->Nsec] - SigmaMed[Density->Nrad - 2]) /
		     SigmaMed[Density->Nrad - 2];
	    vr[cell] = 2. * vr_med - vr[cell - Density->Nsec];
	}
	mean = 0.0;
	for (nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
	    mean += rho[nAzimuthal + Density->Nsec * (Density->Nrad - 1)];
	}
	mean /= (double)(Density->Nsec);
	for (nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
	    rho[nAzimuthal + (Density->Nrad - 1) * Density->Nsec] +=
		SigmaMed[Density->Nrad - 1] - mean;
	}
	mean = 0.0;
	for (nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
	    mean +=
		Energy->Field[nAzimuthal + Density->Nsec * (Density->Nrad - 1)];
	}
	mean /= (double)(Density->Nsec);
	for (nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
	    Energy->Field[nAzimuthal + (Density->Nrad - 1) * Density->Nsec] +=
		EnergyMed[Density->Nrad - 1] - mean;
	}
    }
}

/**

\param VRadial radial velocity polar grid
\param VAzimuthal azimuthal velocity polar grid
\param Density density polar grid
\param Energy energy polar grid
\param step
\param sys planetary system
*/
void EvanescentBoundary(t_data &data, double step)
{
    unsigned int nRadial, nAzimuthal, cell;

    double xp, yp, rp;
    double *vrad, *vtheta, *dens, *energ;
    double vrad0, vtheta0, /*viscosity,*/ dens0, energ0;
    double DRMIN, DRMAX, damping, Tin, Tout, lambda;

    vrad = data[t_data::V_RADIAL].Field;
    vtheta = data[t_data::V_AZIMUTHAL].Field;
    dens = data[t_data::SIGMA].Field;
    energ = data[t_data::ENERGY].Field;

    /* Orbital period at inner and outer boundary */
    Tin = 2.0 * M_PI * pow(GlobalRmed[0], 3. / 2);
    Tout = 2.0 * M_PI * pow(GlobalRmed[GlobalNRadial - 1], 3. / 2);

    /* DRMIN AND DRMAX are global Radii boundaries of killing wave zones */
    xp = data.get_planetary_system().get_planet(0).get_x();
    yp = data.get_planetary_system().get_planet(0).get_y();
    rp = sqrt(xp * xp + yp * yp);
    DRMIN = GlobalRmed[0] + 0.1667 * (rp - GlobalRmed[0]);
    DRMAX = GlobalRmed[GlobalNRadial - 1] -
	    0.2667 * (GlobalRmed[GlobalNRadial - 1] - rp);

    lambda = 0.0;
    for (nRadial = Zero_or_active; nRadial < Max_or_active; ++nRadial) {
	if ((Rmed[nRadial] < DRMIN) || (Rmed[nRadial] > DRMAX)) {
	    /* Damping operates only inside the wave killing zones */
	    if (Rmed[nRadial] < DRMIN) {
		damping = (Rmed[nRadial] - DRMIN) / (GlobalRmed[0] - DRMIN);
		lambda = damping * damping * 10.0 * step / Tin;
	    }
	    if (Rmed[nRadial] > DRMAX) {
		damping = (Rmed[nRadial] - DRMAX) /
			  (GlobalRmed[GlobalNRadial - 1] - DRMAX);
		lambda = damping * damping * 10.0 * step / Tout;
	    }
	    /*
	    if (ViscosityAlpha || (VISCOSITY != 0.0) )
		    viscosity = FViscosity (Rmed[i]);
	    if (!ViscosityAlpha && (VISCOSITY == 0.0) )
		    viscosity = 0.0;
	    if (!SelfGravity) {
		    vtheta0 = sqrt ( G*1.0/Rmed[i] * ( 1.0 -
	    (1.0+SIGMASLOPE-2.0*FLARINGINDEX)*
	    pow(AspectRatio(Rmed[i]),2.0)*pow(Rmed[i],2.0*FLARINGINDEX) ) );
	    }
	    if (SelfGravity) {
		    vtheta0 = sqrt ( G*1.0/Rmed[i] * ( 1.0 -
	    (1.0+SIGMASLOPE-2.0*FLARINGINDEX)*
	    pow(AspectRatio(Rmed[i]),2.0)*pow(Rmed[i],2.0*FLARINGINDEX) ) -
	    Rmed[i]*GLOBAL_AxiSGAccr[i+IMIN] );
	    }
	    // this could be refined if CentrifugalBalance is used...
	    vtheta0 -= Rmed[i]*OmegaFrame;
	    vrad0 = -3.0*viscosity/Rmed[i]*(-SIGMASLOPE+.5);
	    dens0 = SigmaMed[i];
	    energ0 = EnergyMed[i];
	    */
	    vrad0 = 0.0;
	    vtheta0 = 0.0;
	    dens0 = 0.0;
	    energ0 = 0.0;
	    for (nAzimuthal = 0; nAzimuthal < data[t_data::SIGMA].Nsec;
		 ++nAzimuthal) {
		cell = nRadial * data[t_data::SIGMA].Nsec + nAzimuthal;
		vrad0 += vrad[cell];
		vtheta0 += vtheta[cell];
		dens0 += dens[cell];
		energ0 += energ[cell];
	    }
	    vrad0 /= (double)(data[t_data::SIGMA].Nsec);
	    vtheta0 /= (double)(data[t_data::SIGMA].Nsec);
	    dens0 /= (double)(data[t_data::SIGMA].Nsec);
	    energ0 /= (double)(data[t_data::SIGMA].Nsec);

	    for (nAzimuthal = 0; nAzimuthal < data[t_data::SIGMA].Nsec;
		 ++nAzimuthal) {
		cell = nRadial * data[t_data::SIGMA].Nsec + nAzimuthal;
		vrad[cell] = (vrad[cell] + lambda * vrad0) / (1.0 + lambda);
		vtheta[cell] =
		    (vtheta[cell] + lambda * vtheta0) / (1.0 + lambda);
		dens[cell] = (dens[cell] + lambda * dens0) / (1.0 + lambda);
		if (parameters::Adiabatic)
		    energ[cell] =
			(energ[cell] + lambda * energ0) / (1.0 + lambda);
	    }
	}
    }
}

/**
	\param Density
	\param VRadial
*/
void ApplyOuterSourceMass(t_polargrid *Density, t_polargrid *VRadial)
{
    unsigned int nRadial, nAzimuthal, cell;

    double averageRho = 0.0;
    double penul_vr;

    if (CPU_Rank != CPU_Highest)
	return;

    nRadial = Density->Nrad - 1;
    for (nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
	cell = nAzimuthal + nRadial * Density->Nsec;
	averageRho += Density->Field[cell];
    }

    averageRho /= (double)(Density->Nsec);
    averageRho = SigmaMed[Density->Nrad - 1] - averageRho;

    for (nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
	cell = nAzimuthal + nRadial * Density->Nsec;
	Density->Field[cell] += averageRho;
    }

    penul_vr =
	IMPOSEDDISKDRIFT * pow((Rinf[Density->Nrad - 1] / 1.0), -SIGMASLOPE);
    for (nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
	cell = nAzimuthal + nRadial * Density->Nsec;
	VRadial->Field[cell] = penul_vr;
    }
}

/**
	\param VAzimuthal azimuthal velocity polar grid
*/
void ApplySubKeplerianBoundaryInner(t_polargrid &v_azimuthal)
{
    double VKepIn = 0.0;
    if (!parameters::self_gravity) {
	/* (3.4) on page 44 */
	VKepIn = sqrt(constants::G * hydro_center_mass / Rb[0] *
		      (1.0 - (1.0 + SIGMASLOPE - 2.0 * FLARINGINDEX) *
				 pow(ASPECTRATIO_REF, 2.0) *
				 pow(Rb[0], 2.0 * FLARINGINDEX)));
    } else {
	mpi_make1Dprofile(selfgravity::g_radial, GLOBAL_AxiSGAccr);

	/* (3.42) on page 55 */
	/* VKepIn is only needed on innermost CPU */
	if (CPU_Rank == 0) {
	    // viscosity::aspect_ratio(Rmed[0])
	    VKepIn = sqrt(constants::G * hydro_center_mass / Rb[0] *
			      (1.0 - (1.0 + SIGMASLOPE - 2.0 * FLARINGINDEX) *
					 pow(ASPECTRATIO_REF, 2.0) *
					 pow(Rb[0], 2.0 * FLARINGINDEX)) -
			  Rb[0] * GLOBAL_AxiSGAccr[0]);
	}
    }

    if (CPU_Rank == 0) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < v_azimuthal.get_size_azimuthal(); n_azimuthal++) {
	    v_azimuthal(0, n_azimuthal) = VKepIn - Rb[0] * OmegaFrame;
	}
    }
}

/**
	\param VAzimuthal azimuthal velocity polar grid
*/
void ApplySubKeplerianBoundaryOuter(t_polargrid &v_azimuthal, const bool did_sg)
{
    double VKepOut = 0.0;
    if (!parameters::self_gravity) {
	/* (3.4) on page 44 */
	VKepOut = sqrt(constants::G * hydro_center_mass /
		       Rb[v_azimuthal.get_max_radial()] *
		       (1.0 - (1.0 + SIGMASLOPE - 2.0 * FLARINGINDEX) *
				  pow(ASPECTRATIO_REF, 2.0) *
				  pow(Rb[v_azimuthal.get_max_radial()],
				      2.0 * FLARINGINDEX)));
    } else {

	if (!did_sg) {
	    mpi_make1Dprofile(selfgravity::g_radial, GLOBAL_AxiSGAccr);
	}

	/* (3.42) on page 55 */
	/* VKepOut is only needed on outermost CPU */
	if (CPU_Rank == CPU_Highest) {
	    // viscosity::aspect_ratio(Rmed[VAzimuthal->Nrad-1])
	    VKepOut =
		sqrt(constants::G * hydro_center_mass /
			 Rb[v_azimuthal.get_max_radial()] *
			 (1.0 - (1.0 + SIGMASLOPE - 2.0 * FLARINGINDEX) *
				    pow(ASPECTRATIO_REF, 2.0) *
				    pow(Rb[v_azimuthal.get_max_radial()],
					2.0 * FLARINGINDEX)) -
		     Rb[v_azimuthal.get_max_radial()] *
			 GLOBAL_AxiSGAccr[v_azimuthal.get_max_radial() + IMIN]);
	}
    }

    if (CPU_Rank == CPU_Highest) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < v_azimuthal.get_size_azimuthal(); n_azimuthal++) {
	    v_azimuthal(v_azimuthal.get_max_radial(), n_azimuthal) =
		VKepOut - Rb[v_azimuthal.get_max_radial()] * OmegaFrame;
	}
    }
}

void correct_v_azimuthal(t_polargrid &v_azimuthal, double dOmega)
{
    // TODO: maybe think about max-1 here as this might be alread set in
    // boundary condition or change in ApplySubKeplerianBoundary
    for (unsigned int n_radial = 0; n_radial < v_azimuthal.get_size_radial();
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < v_azimuthal.get_size_azimuthal(); ++n_azimuthal) {
	    v_azimuthal(n_radial, n_azimuthal) -= dOmega * Rb[n_radial];
	}
    }
}

/*
void ApplyNoTorqueBoundaryInner(t_polargrid &v_azimuthal)
{
	if (CPU_Rank == 0) {
		for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= v_azimuthal.get_max_azimuthal();
		 ++n_azimuthal) {
		// this is a work around as long as V_AZIMUTHAL is defined as a
		// vector
			const double R_outer = Rmed[1];
			const double R_inner = Rmed[0];
			v_azimuthal(0, n_azimuthal) =
			std::sqrt(R_inner / R_outer) *
			v_azimuthal(1, n_azimuthal);
		}
	}
}

void ApplyNoTorqueBoundaryOuter(t_polargrid &v_azimuthal)
{
	if (CPU_Rank == CPU_Highest) {
		for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= v_azimuthal.get_max_azimuthal();
		 ++n_azimuthal) {
		// this is a work around as long as V_AZIMUTHAL is defined as a
		// vector
			const double R_outer =
Rmed[v_azimuthal.get_max_radial()]; const double R_inner =
Rmed[v_azimuthal.get_max_radial() - 1]; v_azimuthal(
v_azimuthal.get_max_radial(), n_azimuthal) = std::sqrt(R_outer / R_inner) *
			v_azimuthal(v_azimuthal.get_max_radial() - 1,
n_azimuthal);
		}
	}
}
*/
