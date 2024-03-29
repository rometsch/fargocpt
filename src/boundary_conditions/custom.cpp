/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"
#include "../global.h"
#include "../Theo.h"
#include "../frame_of_reference.h"


namespace boundary_conditions
{

/* 
// Place to declare custom variables.
*/

// Declare them here.

/* In this function, you can initialize custom variables for your boundary condition.
// E.g., you might precalculate values or load files.
*/
void init_custom([[maybe_unused]] t_data& data) {

}


/* In this function, you can clean up whatevery you initialized in the custom_init function.
// E.g., you might free memory or close files.
*/
void cleanup_custom() {

}




/* Template to modify for implementing custom boundary conditions.
//
// How to use this template:
// -------------------------
// 1) Comment out all assinments of variables that are set by by other conditions you want to keep.
// 2) Modify the values asigned to the variables you want to change.
//
// For the azimuthal velocity, keep the -r*OmegaF term. This corrects for the rotating frame.
//
*/
void custom_inner([[maybe_unused]] t_data& data, [[maybe_unused]] const double t, [[maybe_unused]] const double dt) {
	if (CPU_Rank != 0) {
		return;
	}

	t_polargrid &Sigma = data[t_data::SIGMA];
	t_polargrid &Energy = data[t_data::ENERGY];
	t_polargrid &vrad = data[t_data::V_RADIAL];
	t_polargrid &vaz = data[t_data::V_AZIMUTHAL];

	// Set inner boundary

	// get highest index in the azimuthal direction
	const unsigned int Iaz = Sigma.get_max_azimuthal();
	const double vKep = compute_v_kepler(Rmed[0], hydro_center_mass);
	const double OmegaF = refframe::OmegaFrame;
	const double r = Rmed[0];

	#pragma omp parallel for
	for (unsigned int naz = 0; naz <= Iaz; ++naz) {
		// set cell interfaces of the ghost cell
		Sigma(0, naz) = 0;
		Energy(0, naz) = 0;
		vaz(0, naz) = vKep - r*OmegaF;
		// set both interfaces of the ghost cell
		vrad(0, naz) = 0;
		vrad(1, naz) = 0; // domain boundary
	}

}

void custom_outer([[maybe_unused]] t_data& data, [[maybe_unused]] const double t, [[maybe_unused]] const double dt) {

	if (CPU_Rank != CPU_Highest) {
		return;
	}

	t_polargrid &Sigma = data[t_data::SIGMA];
	t_polargrid &Energy = data[t_data::ENERGY];
	t_polargrid &vrad = data[t_data::V_RADIAL];
	t_polargrid &vaz = data[t_data::V_AZIMUTHAL];

	// Set outer boundary

	// get highest index in the both directions
	const unsigned int Iaz = Sigma.get_max_azimuthal();
	const unsigned int Irad = Sigma.get_max_radial();

	const double vKep = compute_v_kepler(Rmed[Irad], hydro_center_mass);
	const double OmegaF = refframe::OmegaFrame;
	const double r = Rmed[Irad];

	#pragma omp parallel for
	for (unsigned int naz = 0; naz <= Iaz; ++naz) {
		// set cell interfaces of the ghost cell
		Sigma(Irad, naz) = 0;
		Energy(Irad, naz) = 0;
		vaz(Irad, naz) = vKep - r*OmegaF;
		// set both interfaces of the ghost cell
		vrad(Irad, naz) = 0;
		vrad(Irad-1, naz) = 0; // domain boundary
	}

}

} // namespace boundary_conditions
