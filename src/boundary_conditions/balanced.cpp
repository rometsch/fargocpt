/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"

#include "../Theo.h"
#include "../frame_of_reference.h"
#include "../global.h"
#include "../parameters.h"

namespace boundary_conditions
{

/**
	\param VAzimuthal azimuthal velocity polar grid
*/
static void balanced_boundary(t_polargrid &vaz, const unsigned int nr)
{
    // Reference in this function refer to Clement Baruteau's PhD thesis 
    // https://ui.adsabs.harvard.edu/abs/2008PhDT.......292B

    const double R = Rb[nr];

    const double vk_2 = std::pow(compute_v_kepler(R, hydro_center_mass), 2);

    double support = 0.0;

    // if (!parameters::profile_cutoff_outer) {
    // if we have a cutoff, we assume the pressure is low enough that kepler
    // velocity is appropriate
    // TODO: why do we do this?
    support += support_azi_pressure(R);
    support += support_azi_smoothing_derivative(R);
    // }

    if (parameters::v_azimuthal_with_quadropole_support) {
	/* (3.4) on page 44 */
	support += support_azi_quadrupole(R);
    }

    // compute square of equilibrium azimuthal velocity
    double v_sq = vk_2 * support;

    // correct for self-gravity if enabled
    if (parameters::self_gravity) {
	// if (!did_sg) {
	// 	// TODO: why do we need to communicate something when we only
	// need the outer boundary? 	mpi_make1Dprofile(selfgravity::g_radial,
	// GLOBAL_AxiSGAccr);
	// }
	/* (3.42) on page 55 */
	const double sg_accel_radial = GLOBAL_AxiSGAccr[nr + IMIN];
	const double correction_sg = -R * sg_accel_radial;
	v_sq += correction_sg;
    }

    double vaz_balanced = std::sqrt(v_sq);

    // Correct for rotating frame;
    vaz_balanced -= Rb[nr] * refframe::OmegaFrame;

    // Apply boundary condition
    const unsigned int Nphi = vaz.get_size_azimuthal();
#pragma omp parallel for
    for (unsigned int naz = 0; naz < Nphi; naz++) {
	vaz(nr, naz) = vaz_balanced;
    }
}

void balanced_inner(t_polargrid &vaz, [[maybe_unused]] t_polargrid &dummy,
		    [[maybe_unused]] t_data &ddummy)
{
    if (CPU_Rank != 0) {
	return;
    }

    balanced_boundary(vaz, 0);
}

void balanced_outer(t_polargrid &vaz, [[maybe_unused]] t_polargrid &dummy,
		    [[maybe_unused]] t_data &ddummy)
{
    if (CPU_Rank != CPU_Highest) {
	return;
    }

    const unsigned int Irad = vaz.get_max_radial();
    balanced_boundary(vaz, Irad);
}

} // namespace boundary_conditions
