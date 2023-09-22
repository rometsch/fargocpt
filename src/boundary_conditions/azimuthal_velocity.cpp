/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"

#include "../Theo.h"
#include "../global.h"
#include "../parameters.h"
#include "../frame_of_reference.h"
#include "../constants.h"
#include "../axilib.h"
#include "../selfgravity.h"


namespace boundary_conditions
{

/****************************************
/ Enforce zero shear at the boundaries
****************************************/
void zero_shear_boundary(t_polargrid &vaz) {
	/// Vphi_i / r_i - Vphi_(i-1) / r_(i-1) = 0
	const unsigned int Naz = vaz.get_max_azimuthal();
	const unsigned int Nrad = vaz.get_max_radial();

	if (CPU_Rank == CPU_Highest) {
		#pragma omp parallel for
	    for (unsigned int naz = 0; naz <= Naz; ++naz) {
			const double v_active = vaz(Nrad - 1, naz);
			const double r = Rmed[Nrad];
			const double ri = Rmed[Nrad - 1];

			vaz(Nrad, naz) = r / ri * v_active;
	    }
	}

	if (CPU_Rank == 0) {
		#pragma omp parallel for
	    for (unsigned int naz = 0; naz <= Naz; ++naz) {
			const double v_active = vaz(1, naz);
			const double ro = Rmed[1];
			const double r = Rmed[0];

			vaz(0, naz) = r / ro * v_active;
	    }
	}

}


void ApplyKeplerianBoundaryInner(t_polargrid &vaz)
{
    double vaz_in = 0.0;
    if (!parameters::self_gravity) {
    // we assume the pressure is low enough that kepler velocity is appropriate
	vaz_in = compute_v_kepler(Rb[0], hydro_center_mass);

    } else {
	mpi_make1Dprofile(selfgravity::g_radial, GLOBAL_AxiSGAccr);

	/* (3.42) on page 55 */
	/* vaz_in is only needed on innermost CPU */
	if (CPU_Rank == 0) {
        // we assume the pressure is low enough that kepler velocity is appropriate
		const double R = Rb[0];
		const double vk_2 = constants::G * hydro_center_mass / R;
        vaz_in = std::sqrt(vk_2 - R * GLOBAL_AxiSGAccr[0]);
	}
    }

    if (CPU_Rank == 0) {
	const unsigned int Nphi = vaz.get_size_azimuthal();
	#pragma omp parallel for
	for (unsigned int naz = 0; naz < Nphi; naz++) {
		vaz(0, naz) = vaz_in - Rb[0] * refframe::OmegaFrame;
	}
    }
}

/**
	\param VAzimuthal azimuthal velocity polar grid
*/
void ApplySubKeplerianBoundaryOuter(t_polargrid &vaz, const bool did_sg)
{
    if (CPU_Rank != CPU_Highest) {
		return;
	}

	const unsigned int nr = vaz.get_max_radial();
	const double R = Rb[nr];
	const double vk_2 = constants::G * hydro_center_mass / R;

	double support = 0.0;
	
	if (!parameters::profile_cutoff_outer) {
		// if we have a cutoff, we assume the pressure is low enough that kepler velocity is appropriate
		// TODO: why do we do this?
		support += support_azi_pressure(R);
		support += support_azi_smoothing_derivative(R);
	}
	
	if(parameters::v_azimuthal_with_quadropole_support){
		/* (3.4) on page 44 */
		support += support_azi_quadrupole(R);
	}

	// compute square of equilibrium azimuthal velocity
	double v_sq = vk_2 * support;

	// correct for self-gravity if enabled
	if (parameters::self_gravity) {
		if (!did_sg) {
			// TODO: why do we need to communicate something when we only need the outer boundary?
			mpi_make1Dprofile(selfgravity::g_radial, GLOBAL_AxiSGAccr);
		}
		/* (3.42) on page 55 */
		const double sg_accel_radial = GLOBAL_AxiSGAccr[nr + IMIN];
		const double correction_sg = -R * sg_accel_radial;
		v_sq += correction_sg;
	}

	double vaz_out = std::sqrt(v_sq);

	// Correct for rotating frame;
	vaz_out -= Rb[nr] * refframe::OmegaFrame;

	// Apply boundary condition
	const unsigned int Nphi = vaz.get_size_azimuthal();
	#pragma omp parallel for
	for (unsigned int naz = 0; naz < Nphi; naz++) {
		vaz(nr, naz) =	vaz_out;
	}

}

} // namespace boundary_conditions
