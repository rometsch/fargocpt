/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"

#include "../Theo.h"
#include "../find_cell_id.h"
#include "../global.h"
#include "../logging.h"
#include "../parameters.h"
#include "../util.h"
#include "../frame_of_reference.h"
#include "../simulation.h"
#include "../constants.h"
#include "../quantities.h"
#include "../axilib.h"
#include "../selfgravity.h"

#include <algorithm>
#include <cstring>
#include <cmath>
#include <vector>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include "viscosity/viscous_radial_speed.h"


// temporary
#include "SourceEuler.h"
#include "LowTasks.h"
#include "SideEuler.h"

namespace boundary_conditions
{

/**
	Boundary conditions for calculation of the boundary layer (BL) starting
   here:
	TODO: Stellar radiative flux into disk via implicit routine!
*/

/**
	Inner boundary: zero gradient & fixed velocities
*/

void boundary_layer_inner_boundary(t_data &data)
{

    if (CPU_Rank != 0)
	return;

	#pragma omp parallel for
    for (unsigned int n_azimuthal = 0;
	 n_azimuthal <= data[t_data::SIGMA].get_max_azimuthal();
	 ++n_azimuthal) {
	// zero gradient
	data[t_data::SIGMA](0, n_azimuthal) =
	    data[t_data::SIGMA](1, n_azimuthal);
	data[t_data::ENERGY](0, n_azimuthal) =
	    data[t_data::ENERGY](1, n_azimuthal);

	// set vrad to fraction of Keplerian velocity
	data[t_data::V_RADIAL](1, n_azimuthal) =
	    -1. * parameters::vrad_fraction_of_kepler * std::sqrt(1. / Ra[1]);
	data[t_data::V_RADIAL](0, n_azimuthal) =
	    data[t_data::V_RADIAL](1, n_azimuthal);

	// set vphi to stellar rotation rate 
	data[t_data::V_AZIMUTHAL](0, n_azimuthal) =
	    parameters::stellar_rotation_rate * std::sqrt(1. / Rb[0]);
    }
}

/**
	Outer boundary: floating boundary conditions & pressure correction for
   Omega
*/

void boundary_layer_outer_boundary(t_data &data)
{

    if (CPU_Rank != CPU_Highest)
	return;

	#pragma omp parallel for
    for (unsigned int n_azimuthal = 0;
	 n_azimuthal <= data[t_data::SIGMA].get_max_azimuthal();
	 ++n_azimuthal) {
	// floating BCs
	data[t_data::SIGMA](data[t_data::SIGMA].get_max_radial(), n_azimuthal) =
	    data[t_data::SIGMA](data[t_data::SIGMA].get_max_radial() - 1,
				n_azimuthal) *
	    std::sqrt(Ra[data[t_data::SIGMA].get_max_radial() - 1] /
		      Ra[data[t_data::SIGMA].get_max_radial()]);
	data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial(),
			     n_azimuthal) =
	    data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial() - 1,
				 n_azimuthal) *
	    std::pow(Ra[data[t_data::ENERGY].get_max_radial() - 1] /
			 Ra[data[t_data::ENERGY].get_max_radial()],
		     1.25);

	data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial() - 1,
			       n_azimuthal) =
	    -1. *
	    fabs(data[t_data::V_RADIAL](
		data[t_data::V_RADIAL].get_max_radial() - 2, n_azimuthal)) *
	    std::sqrt(Ra[data[t_data::V_RADIAL].get_max_radial() - 2] /
		      Ra[data[t_data::V_RADIAL].get_max_radial() - 1]);

	data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial(),
			       n_azimuthal) =
	    -1. *
	    fabs(data[t_data::V_RADIAL](
		data[t_data::V_RADIAL].get_max_radial() - 2, n_azimuthal)) *
	    std::sqrt(Ra[data[t_data::V_RADIAL].get_max_radial() - 2] /
		      Ra[data[t_data::V_RADIAL].get_max_radial()]);

	// Omega at outer boundary equals calculate_omega_kepler (plus leading
	// order pressure correction)
	data[t_data::V_AZIMUTHAL](data[t_data::V_AZIMUTHAL].get_max_radial(),
				  n_azimuthal) =
	    1. / std::sqrt(Rb[data[t_data::SIGMA].get_max_radial()]);
	// TODO: Include pressure correction, like in uphi[*jN]
	// = 1./sqrt(Rb[*jN]) +
	// 0.5/Sigma[*jN]*sqrt(pow3(Rb[*jN])*pow2(Rb[*jN]))*.5/Rb[*jN]*(P[*jN+1]-P[*jN-1])/DeltaRa[*jN+1];
    }
}

} // namespace boundary_conditions
