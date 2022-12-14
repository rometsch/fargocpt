/*
	\file SideEuler.c

	Total mass and angular momentum monitoring, and boundary conditions. In
   addition, this file contains a few low-level functions that manipulate
   PolarGrid 's or initialize the forces evaluation.

*/
#ifdef _OPENMP
#include <omp.h>
#endif
#include <mpi.h>

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
#include "nongnu.h"
#include "parameters.h"
#include "quantities.h"
#include "selfgravity.h"
#include "util.h"
#include "viscosity/viscosity.h"
#include "frame_of_reference.h"
#include <cassert>

extern boolean Damping;
extern boolean OuterSourceMass;


/**
	Divide each cell of a polar grid with a cell of another polargrid

	\param Num numerator
	\param Denom denominator
	\param Result result
*/
void divise_polargrid(const t_polargrid &num, const t_polargrid &denom, t_polargrid &result)
{
    const unsigned int Nmax =
	result.get_size_radial() * result.get_size_azimuthal();

	#pragma omp parallel for
    for (unsigned int n = 0; n < Nmax; n++) {
	assert(denom.Field[n] > 0.0);
	/// denom + DBL_EPSILON can cause problems because DBL_EPSILON can be
	/// bigger than denom, depending on units.
	result.Field[n] = num.Field[n] / denom.Field[n]; /// in case of crash, use something like
			    /// (denom.Field[n] + 1.0e-200) instead.
    }
}

/**

*/
void InitCellCenterCoordinates()
{
    delete CellCenterY;
    delete CellCenterX;

    CellCenterX = CreatePolarGrid(NRadial + 1, NAzimuthal, "cell_center_x");
    CellCenterY = CreatePolarGrid(NRadial + 1, NAzimuthal, "cell_center_y");

	const unsigned int Nr = CellCenterX->Nrad;
	const unsigned int Nphi = CellCenterX->Nsec;

	#pragma omp parallel for collapse(2)
	for (unsigned int nRadial = 0; nRadial < Nr; ++nRadial) {
	for (unsigned int nAzimuthal = 0; nAzimuthal < Nphi; ++nAzimuthal) {
		unsigned int cell = nAzimuthal + nRadial * Nphi;
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

	#pragma omp parallel for
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
	#pragma omp parallel for
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

    auto & vrad = data[t_data::V_RADIAL];
    auto & vtheta = data[t_data::V_AZIMUTHAL];
    auto & dens = data[t_data::SIGMA];
    auto & energ = data[t_data::ENERGY];

    /* Orbital period at inner and outer boundary */
    const double Tin = 2.0 * M_PI * pow(GlobalRmed[0], 3. / 2);
    const double Tout = 2.0 * M_PI * pow(GlobalRmed[GlobalNRadial - 1], 3. / 2);

    /* DRMIN AND DRMAX are global Radii boundaries of killing wave zones */
    const double xp = data.get_planetary_system().get_planet(1).get_x();
    const double yp = data.get_planetary_system().get_planet(1).get_y();
    const double rp = sqrt(xp * xp + yp * yp);
    const double R_damp_in = GlobalRmed[0] + 0.1667 * (rp - GlobalRmed[0]);
    const double R_damp_out = GlobalRmed[GlobalNRadial - 1] - 0.2667 * (GlobalRmed[GlobalNRadial - 1] - rp);

	#pragma omp parallel for
    for (unsigned int nrad = Zero_or_active; nrad < Max_or_active; ++nrad) {
	if ((Rmed[nrad] < R_damp_in) || (Rmed[nrad] > R_damp_out)) {
		double lambda = 0.0;

	    /* Damping operates only inside the wave killing zones */
	    if (Rmed[nrad] < R_damp_in) {
		const double Rin = GlobalRmed[0];
		const double damping = (Rmed[nrad] - R_damp_in) / (Rin - R_damp_in);
		lambda = damping * damping * 10.0 * step / Tin;
	    }
	    if (Rmed[nrad] > R_damp_out) {
		const double Rout = GlobalRmed[GlobalNRadial - 1];
		const double damping = (Rmed[nrad] - R_damp_out) / (Rout - R_damp_out);
		lambda = damping * damping * 10.0 * step / Tout;
	    }

	    double vrad_mean = 0.0;
	    double vaz_mean = 0.0;
	    double dens_mean = 0.0;
	    double energ_mean = 0.0;
		// Calculate azimuthal averages
	    for (unsigned int naz = 0; naz < data[t_data::SIGMA].Nsec;
		 ++naz) {
		vrad_mean += vrad(nrad, naz);
		vaz_mean += vtheta(nrad, naz);
		dens_mean += dens(nrad, naz);
		energ_mean += energ(nrad, naz);
	    }
	    vrad_mean /= (double)(data[t_data::SIGMA].Nsec);
	    vaz_mean /= (double)(data[t_data::SIGMA].Nsec);
	    dens_mean /= (double)(data[t_data::SIGMA].Nsec);
	    energ_mean /= (double)(data[t_data::SIGMA].Nsec);

		// Apply damping
	    for (unsigned int naz = 0; naz < data[t_data::SIGMA].Nsec;
		 ++naz) {
			vrad(nrad, naz) = (vrad(nrad, naz) + lambda * vrad_mean) / (1.0 + lambda);
			vtheta(nrad, naz) = (vtheta(nrad, naz) + lambda * vaz_mean) / (1.0 + lambda);
			dens(nrad, naz) = (dens(nrad, naz) + lambda * dens_mean) / (1.0 + lambda);
			if (parameters::Adiabatic){
				energ(nrad, naz) = (energ(nrad, naz) + lambda * energ_mean) / (1.0 + lambda);
			}
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
    double averageRho = 0.0;

	if (CPU_Rank != CPU_Highest){
	return;
	}

	const unsigned int Nphi = Density->Nsec;
	const unsigned int nRadial = Density->Nrad - 1;

	#pragma omp parallel for reduction(+ : averageRho)
	for (unsigned int nAzimuthal = 0; nAzimuthal < Nphi; ++nAzimuthal) {
	const unsigned int cell = nAzimuthal + nRadial * Density->Nsec;
	averageRho += Density->Field[cell];
    }

	averageRho /= (double)(Nphi);
	averageRho = SigmaMed[nRadial] - averageRho;

	#pragma omp parallel for
	for (unsigned int nAzimuthal = 0; nAzimuthal < Nphi; ++nAzimuthal) {
	const unsigned int cell = nAzimuthal + nRadial * Density->Nsec;
	Density->Field[cell] += averageRho;
    }

	const double penul_vr =
	parameters::IMPOSEDDISKDRIFT * std::pow((Rinf[nRadial] / 1.0), -parameters::SIGMASLOPE);

	#pragma omp parallel for
	for (unsigned int nAzimuthal = 0; nAzimuthal < Nphi; ++nAzimuthal) {
	const unsigned int cell = nAzimuthal + nRadial * Density->Nsec;
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
	// if we have a cutoff, we assume the pressure is low enough that kepler velocity is appropriate
	if (parameters::profile_cutoff_inner) {
	VKepIn = compute_v_kepler(Rb[0], hydro_center_mass);
	} else {
	/* (3.4) on page 44 */
	VKepIn = initial_locally_isothermal_v_az(Rb[0], hydro_center_mass);
	}

    } else {
	mpi_make1Dprofile(selfgravity::g_radial, GLOBAL_AxiSGAccr);

	/* (3.42) on page 55 */
	/* VKepIn is only needed on innermost CPU */
	if (CPU_Rank == 0) {
		// if we have a cutoff, we assume the pressure is low enough that kepler velocity is appropriate
		const double R = Rb[0];
		if (parameters::profile_cutoff_inner) {
			const double vk_2 = constants::G * hydro_center_mass / R;
			VKepIn = std::sqrt(vk_2 - R * GLOBAL_AxiSGAccr[0]);
		} else {
		const double h0 = parameters::ASPECTRATIO_REF;
		const double F = parameters::FLARINGINDEX;
		const double S = parameters::SIGMASLOPE;
		const double h = h0 * std::pow(R, F);
		//const double eps = parameters::thickness_smoothing;
		const double vk_2 = constants::G * hydro_center_mass / R;
		const double pressure_support_2 = (2.0 * F - 1.0 - S) * std::pow(h, 2);

		//const double smoothing_derivative_2 = (1.0 + (F+1.0) * std::pow(h * eps, 2))
		//		/ std::sqrt(1 + std::pow(h * eps, 2));
		const double smoothing_derivative_2 = 1.0;

		VKepIn = std::sqrt(vk_2 * (smoothing_derivative_2 + pressure_support_2) - R * GLOBAL_AxiSGAccr[0]);
		}
	}
    }

    if (CPU_Rank == 0) {
	const unsigned int Nphi = v_azimuthal.get_size_azimuthal();
	#pragma omp parallel for
	for (unsigned int naz = 0; naz < Nphi; naz++) {
		v_azimuthal(0, naz) = VKepIn - Rb[0] * refframe::OmegaFrame;
	}
    }
}

/**
	\param VAzimuthal azimuthal velocity polar grid
*/
void ApplySubKeplerianBoundaryOuter(t_polargrid &v_azimuthal, const bool did_sg)
{
    double VKepOut = 0.0;
	const unsigned int nr = v_azimuthal.get_max_radial();

    if (!parameters::self_gravity) {
	// if we have a cutoff, we assume the pressure is low enough that kepler velocity is appropriate
	if(parameters::profile_cutoff_outer){
	VKepOut = compute_v_kepler(Rb[nr], hydro_center_mass);
	} else {
	/* (3.4) on page 44 */
	VKepOut = initial_locally_isothermal_v_az(Rb[nr], hydro_center_mass);
	}
    } else {

	if (!did_sg) {
	    mpi_make1Dprofile(selfgravity::g_radial, GLOBAL_AxiSGAccr);
	}

	/* (3.42) on page 55 */
	/* VKepOut is only needed on outermost CPU */
	if (CPU_Rank == CPU_Highest) {
		// if we have a cutoff, we assume the pressure is low enough that kepler velocity is appropriate
		const double R = Rb[nr];
		if(parameters::profile_cutoff_outer){
			const double vk_2 = constants::G * hydro_center_mass / R;
			VKepOut = std::sqrt(vk_2 - R * GLOBAL_AxiSGAccr[nr + IMIN]);
		} else {
		const double h0 = parameters::ASPECTRATIO_REF;
		const double F = parameters::FLARINGINDEX;
		const double S = parameters::SIGMASLOPE;
		const double h = h0 * std::pow(R, F);
		//const double eps = parameters::thickness_smoothing;
		const double vk_2 = constants::G * hydro_center_mass / R;
		const double pressure_support_2 = (2.0 * F - 1.0 - S) * std::pow(h, 2);
		//const double smoothing_derivative_2 = (1.0 + (F+1.0) * std::pow(h * eps, 2))
		//		/ std::pow(std::sqrt(1 + std::pow(h * eps, 2)), 3);
		const double smoothing_derivative_2 = 1.0;

		VKepOut = std::sqrt(vk_2 * (smoothing_derivative_2 + pressure_support_2) - R * GLOBAL_AxiSGAccr[nr + IMIN]);
		}
	}
    }

    if (CPU_Rank == CPU_Highest) {
	const unsigned int Nphi = v_azimuthal.get_size_azimuthal();

	#pragma omp parallel for
	for (unsigned int naz = 0; naz < Nphi; naz++) {
		v_azimuthal(nr, naz) =	VKepOut - Rb[nr] * refframe::OmegaFrame;
	}
    }
}

void correct_v_azimuthal(t_polargrid &v_azimuthal, double dOmega)
{
    // TODO: maybe think about max-1 here as this might be alread set in
    // boundary condition or change in ApplySubKeplerianBoundary

	// We update velocities from old OmegaFrame to new OmegaFrame;
	// As the ghost cells should have the old OmegaFrame, we need to update them here aswell.
	const unsigned int Nr = v_azimuthal.get_size_radial();
	const unsigned int Nphi = v_azimuthal.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		v_azimuthal(nr, naz) -= dOmega * Rb[nr];
	}
    }
}

/*
void ApplyKeplerExtraplationBoundaryInner(t_polargrid &v_azimuthal)
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
			std::sqrt(R_outer / R_inner) *
			v_azimuthal(1, n_azimuthal);
		}
	}
}

void ApplyKeplerExtrapolationBoundaryOuter(t_polargrid &v_azimuthal)
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
v_azimuthal.get_max_radial(), n_azimuthal) = std::sqrt(R_inner / R_outer) *
			v_azimuthal(v_azimuthal.get_max_radial() - 1,
n_azimuthal);
		}
	}
}
*/
