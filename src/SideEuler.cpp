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

#include "LowTasks.h"
#include "SideEuler.h"
#include "global.h"
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
	\param VAzimuthal azimuthal velocity polar grid
*/

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
