/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"
#include "../global.h"
#include "../parameters.h"

namespace boundary_conditions
{


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

    /// TODO OMP cannot share variables between threads
    //#pragma omp parallel for
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


} // namespace boundary_conditions
