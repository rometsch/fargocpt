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
#include "../simulation.h"

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
    if (CPU_Rank != 0) {
		return;
	}

    unsigned int nRadial, cell;

    unsigned int jp;
    int lp;
    double *rho, *vr;
    double mean, vr_med;
    double cs0, cs1;

	t_polargrid &cs = data[t_data::SOUNDSPEED];
	t_polargrid &Sigma = data[t_data::SIGMA];
	unsigned int Naz = Sigma.get_max_azimuthal();


    rho = Density->Field;
    vr = VRadial->Field;

	// Compute the average sound speed in the first ring
	cs0 = 0.0;
	cs1 = 0.0;
	for (unsigned int naz = 0; naz <= Naz; ++naz) {
	    cs0 += cs(0, naz);
	    cs1 += cs(1, naz);
	}
	cs0 /= (Naz + 1);
	cs1 /= (Naz + 1);
	nRadial = 1;

	// The expression below should be refined
	// We need to know the orbital frequency of the nearest planet
	double dangle = (pow(Rinf[1], -1.5) - 1.0) / (.5 * (cs0 + cs1));
	dangle *= (Rmed[1] - Rmed[0]);
	double i_angle = (int)(dangle / 2.0 / M_PI * (double)NAzimuthal + .5);

    /// TODO OMP cannot share variables between threads
    //#pragma omp parallel for
	for (unsigned int naz = 0; naz <= Naz; ++naz) {
	    cell = naz + nRadial * Density->Nsec;
	    jp = naz + i_angle;
	    if (jp >= Density->Nsec)
		jp -= Density->Nsec;
	    /*if (jp < 0)
		    jp += Density->Nsec;*/
	    lp = jp;
	    rho[lp] = rho[cell]; // copy first ring into ghost ring
	    Energy->Field[lp] =
		Energy->Field[cell]; // copy first ring into ghost ring
	    vr_med = -cs.Field[cell] *
		     (rho[cell] - SigmaMed[1]) / SigmaMed[1];
	    vr[cell] = 2. * vr_med - vr[cell + Density->Nsec];
	}
	mean = 0.0;
	for (unsigned int naz = 0; naz <= Naz; ++naz) {
	    mean += rho[naz];
	}
	mean /= (double)(Density->Nsec);
	for (unsigned int naz = 0; naz <= Naz; ++naz) {
	    rho[naz] += SigmaMed[0] - mean;
	}
	mean = 0.0;
	for (unsigned int naz = 0; naz <= Naz; ++naz) {
	    mean += Energy->Field[naz];
	}
	mean /= (double)Density->Nsec;
	for (unsigned int naz = 0; naz <= Naz; ++naz) {
	    Energy->Field[naz] += EnergyMed[0] - mean;
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

	t_polargrid &Sigma = data[t_data::SIGMA];
	t_polargrid &cs = data[t_data::SOUNDSPEED];
	unsigned int Naz = Sigma.get_max_azimuthal();


    unsigned int nRadial, cell;

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
	for (unsigned int naz = 0; naz <= Naz; ++naz) {
	    csnrm2 +=
		cs
		    .Field[(Density->Nrad - 2) * Density->Nsec + naz];
	    csnrm1 +=
		cs
		    .Field[(Density->Nrad - 1) * Density->Nsec + naz];
	}
	csnrm1 /= (double)(Density->Nsec);
	csnrm2 /= (double)(Density->Nsec);
	nRadial = Density->Nrad - 1; // The expression below should be refined
	// We need to know the orbital frequency of the nearest planet
	dangle = (pow(Rinf[Density->Nrad - 2], -1.5) - 1.0) /
		 (.5 * (csnrm1 + csnrm2));
	dangle *= (Rmed[Density->Nrad - 1] - Rmed[Density->Nrad - 2]);
	i_angle = (int)(dangle / 2.0 / M_PI * (double)NAzimuthal + .5);

	for (unsigned int naz = 0; naz <= Naz; ++naz) {
	    cell = naz + nRadial * Density->Nsec;
	    jp = naz - i_angle;
	    if (jp >= Density->Nsec)
		jp -= Density->Nsec;
	    /*if (jp < 0)
		    jp += Density->Nsec;*/
	    lp = jp + (nRadial - 1) * Density->Nsec;
	    rho[cell] = rho[lp]; // copy first ring into ghost ring
	    Energy->Field[cell] =
		Energy->Field[lp]; // copy first ring into ghost ring
	    vr_med = cs.Field[cell] *
		     (rho[cell - Density->Nsec] - SigmaMed[Density->Nrad - 2]) /
		     SigmaMed[Density->Nrad - 2];
	    vr[cell] = 2. * vr_med - vr[cell - Density->Nsec];
	}
	mean = 0.0;
	for (unsigned int naz = 0; naz <= Naz; ++naz) {
	    mean += rho[naz + Density->Nsec * (Density->Nrad - 1)];
	}
	mean /= (double)(Density->Nsec);
	for (unsigned int naz = 0; naz <= Naz; ++naz) {
	    rho[naz + (Density->Nrad - 1) * Density->Nsec] +=
		SigmaMed[Density->Nrad - 1] - mean;
	}
	mean = 0.0;
	for (unsigned int naz = 0; naz <= Naz; ++naz) {
	    mean +=
		Energy->Field[naz + Density->Nsec * (Density->Nrad - 1)];
	}
	mean /= (double)(Density->Nsec);
	for (unsigned int naz = 0; naz <= Naz; ++naz) {
	    Energy->Field[naz + (Density->Nrad - 1) * Density->Nsec] +=
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
    auto & Sigma = data[t_data::SIGMA];
    auto & Energy = data[t_data::ENERGY];

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
		dens_mean += Sigma(nrad, naz);
		energ_mean += Energy(nrad, naz);
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
			Sigma(nrad, naz) = (Sigma(nrad, naz) + lambda * dens_mean) / (1.0 + lambda);
			if (parameters::Adiabatic){
				Energy(nrad, naz) = (Energy(nrad, naz) + lambda * energ_mean) / (1.0 + lambda);
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
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
	const unsigned int cell = naz + nRadial * Density->Nsec;
	averageRho += Density->Field[cell];
    }

	averageRho /= (double)(Nphi);
	averageRho = SigmaMed[nRadial] - averageRho;

	#pragma omp parallel for
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
	const unsigned int cell = naz + nRadial * Density->Nsec;
	Density->Field[cell] += averageRho;
    }

	const double penul_vr =
	parameters::IMPOSEDDISKDRIFT * std::pow((Rinf[nRadial] / 1.0), -parameters::SIGMASLOPE);

	#pragma omp parallel for
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
	const unsigned int cell = naz + nRadial * Density->Nsec;
	VRadial->Field[cell] = penul_vr;
    }
}


void zero_gradient_boundary_inner(t_data &data)
{
    if (CPU_Rank != 0){
		return;
	}

	t_polargrid &sig = data[t_data::SIGMA];
	t_polargrid &eng = data[t_data::ENERGY];
	t_polargrid &vrad = data[t_data::V_RADIAL];

	const unsigned int Naz = sig.get_max_azimuthal();

	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {
		// copy first ring into ghost ring
		sig(0, naz) = sig(1, naz);
		eng(0, naz) = eng(1, naz);

		vrad(1, naz) = vrad(2, naz);
		vrad(0, naz) = vrad(2, naz);
    }
}

void zero_gradient_boundary_outer(t_data &data)
{
    if (CPU_Rank != CPU_Highest) {
		return;
	}

	t_polargrid &sig = data[t_data::SIGMA];
	t_polargrid &eng = data[t_data::ENERGY];
	t_polargrid &vrad = data[t_data::V_RADIAL];

	const unsigned int Nrad = sig.get_max_radial();
	const unsigned int Nrad_v = vrad.get_max_radial();
	const unsigned int Naz = sig.get_max_azimuthal();


	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {
	// copy last ring into ghost ring
	sig(Nrad, naz) = sig(Nrad - 1, naz);
	eng(Nrad, naz) = eng(Nrad - 1, naz);

	vrad(Nrad_v - 1, naz) = vrad(Nrad_v - 2, naz);
	vrad(Nrad_v, naz) = vrad(Nrad_v - 2, naz);
    }
}


void viscous_outflow_boundary_inner(t_data &data)
{
    if (CPU_Rank != 0) {
		return;
	}

	t_polargrid &sig = data[t_data::SIGMA];
	t_polargrid &eng = data[t_data::ENERGY];
	t_polargrid &vrad = data[t_data::V_RADIAL];
	t_polargrid &visc = data[t_data::VISCOSITY];

	const unsigned int Naz = sig.get_max_azimuthal();

	for (unsigned int naz = 0; naz <= Naz; ++naz) {
	    sig(0, naz) = sig(1, naz);
	    eng(0, naz) = eng(1, naz);
	}

	const double s = parameters::viscous_outflow_speed;

	#pragma omp parallel for
	for (unsigned int naz = 0; naz <= Naz; ++naz) {

	    const double Nu0 = visc(0, naz);
	    const double Nu1 = visc(1, naz);
	    const double Nu = 0.5 * (Nu0 + Nu1);

	    // V_rad =  - 1.5 / r * Nu (Kley, Papaloizou and Ogilvie, 2008)
	    vrad(1, naz) = -1.5 * s / Rinf[1] * Nu;
	    vrad(0, naz) = -1.5 * s / Rinf[0] * Nu;
	}
}




/**
	inner reflecting boundary condition
*/
void reflecting_boundary_inner(t_data &data)
{
    if (CPU_Rank != 0) {
		return;
	}

	t_polargrid &sig = data[t_data::SIGMA];
	t_polargrid &eng = data[t_data::ENERGY];
	t_polargrid &vr = data[t_data::V_RADIAL];

	const unsigned int Naz = sig.get_max_azimuthal();

	#pragma omp parallel for
	for (unsigned int naz = 0; naz <= Naz; ++naz) {
	    // copy first ring into ghost ring
	    sig(0, naz) = sig(1, naz);
	    eng(0, naz) = eng(1, naz);

	    vr(1, naz) = 0.0;
	    vr(0, naz) = -vr(2, naz);
	}
}

/**
	outer reflecting boundary condition
*/
void reflecting_boundary_outer(t_data &data)
{
    if (CPU_Rank != CPU_Highest) {
		return;
	}

	t_polargrid &sig = data[t_data::SIGMA];
	t_polargrid &eng = data[t_data::ENERGY];
	t_polargrid &vr = data[t_data::V_RADIAL];

	const unsigned int Nrad = sig.get_max_radial();
	const unsigned int Nrad_v = vr.get_max_radial();
	const unsigned int Naz = sig.get_max_azimuthal();


	#pragma omp parallel for
	for (unsigned int naz = 0; naz <= Naz; ++naz) {
	    // copy last ring into ghost ring
	    sig(Nrad, naz) = sig(Nrad - 1, naz);
	    eng(Nrad, naz) = eng(Nrad - 1, naz);

	    vr(Nrad_v - 1, naz) = 0.0;
	    vr(Nrad_v, naz) = -vr(Nrad_v - 2, naz);
	}
}




/**
	set inner boundary values to initial values
 */
void reference_value_boundary_inner(t_data &data)
{
    if (CPU_Rank != 0) {
	return;
    }
    if (sim::N_hydro_iter == 0) {
	return;
    }

    auto &sig = data[t_data::SIGMA];
    auto &eng = data[t_data::ENERGY];
    auto &vr = data[t_data::V_RADIAL];

    auto &sig0 = data[t_data::SIGMA0];
    auto &eng0 = data[t_data::ENERGY0];
    auto &vr0 = data[t_data::V_RADIAL0];

    const unsigned int Naz = sig.get_max_azimuthal();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {
	sig(0, naz) = sig0(0, naz);
	eng(0, naz) = eng0(0, naz);
	vr(0, naz) = vr0(0, naz);
	vr(1, naz) = vr0(1, naz);
    }
}

/**
	set outer boundary values to initial values
 */
void reference_value_boundary_outer(t_data &data)
{
    if (CPU_Rank != CPU_Highest) {
	return;
    }
    if (sim::N_hydro_iter == 0) {
	return;
    }

    auto &sig = data[t_data::SIGMA];
    auto &eng = data[t_data::ENERGY];
    auto &vr = data[t_data::V_RADIAL];

    auto &sig0 = data[t_data::SIGMA0];
    auto &eng0 = data[t_data::ENERGY0];
    auto &vr0 = data[t_data::V_RADIAL0];

    const unsigned int Naz = sig.get_max_azimuthal();

#pragma omp parallel for
    for (unsigned int naz = 0; naz <= Naz; ++naz) {
	const unsigned int nr = sig.get_max_radial();
	sig(nr, naz) = sig0(nr, naz);
	eng(nr, naz) = eng0(nr, naz);
	const unsigned int nrv = vr.get_max_radial();
	vr(nrv, naz) = vr0(nrv, naz);
	vr(nrv - 1, naz) = vr0(nrv - 1, naz);
    }
}


} // namespace boundary_conditions
