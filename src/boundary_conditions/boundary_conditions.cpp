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
extern boolean OuterSourceMass;

namespace boundary_conditions
{

void apply_boundary_condition(t_data &data, const double current_time, const double dt, const bool final)
{
    // if this is the final boundary condition and damping is enable, do it
    if (final && parameters::damping) {
	damping(data, dt);

	if(ECC_GROWTH_MONITOR){
		quantities::calculate_disk_delta_ecc_peri(data, delta_ecc_damp, delta_peri_damp);
	}
    }

    // inner boundary
    switch (parameters::boundary_inner) {
	case parameters::boundary_condition_initial:
	initial_boundary_inner(data);
	break;
    case parameters::boundary_condition_open:
	open_boundary_inner(data);
	break;
    case parameters::boundary_condition_reflecting:
	reflecting_boundary_inner(data);
	break;
    case parameters::boundary_condition_zero_gradient:
	zero_gradient_boundary_inner(data);
	break;
    case parameters::boundary_condition_boundary_layer:
	boundary_layer_inner_boundary(data);
	break;
    case parameters::boundary_condition_nonreflecting:
	NonReflectingBoundary_inner(data, &data[t_data::V_RADIAL],
				    &data[t_data::SIGMA],
				    &data[t_data::ENERGY]);
	break;
    case parameters::boundary_condition_viscous_outflow:
	viscous_outflow_boundary_inner(data);
	break;
	case parameters::boundary_condition_center_of_mass_initial:
	if (final) {
		damping_initial_center_of_mass_inner(data, dt);
	}
	initial_center_of_mass_boundary_inner(data);
	break;
    case parameters::boundary_condition_jibin_spreading_ring:
	spreading_ring_inner(data);
	break;
    case parameters::boundary_condition_evanescent: // evanescent works only for
						    // inner and outer together
						    // until now
	if (parameters::boundary_outer ==
	    parameters::boundary_condition_evanescent) {
	    EvanescentBoundary(data, dt);
	} else {
	    logging::print_master(
		LOG_ERROR
		"Different EvanescentBoundary Parameters. Old parameter file?\n");
	    die("inner/outer evanescent boundary not implemented yet");
	}
	break;
    case parameters::boundary_condition_keplerian:
	keplerian2d_boundary_inner(data);
	break;
    case parameters::boundary_condition_precribed_time_variable:
	die("Inner precribed time variable boundary condition is not implemented yet!\n");
	break;
    }

    // outer boundary
    switch (parameters::boundary_outer) {
	case parameters::boundary_condition_initial:
	initial_boundary_outer(data);
	break;
    case parameters::boundary_condition_open:
	open_boundary_outer(data);
	break;
    case parameters::boundary_condition_reflecting:
	reflecting_boundary_outer(data);
	break;
    case parameters::boundary_condition_center_of_mass_initial: {
	if (final) {
	    damping_initial_center_of_mass_outer(data, dt);
	}
	initial_center_of_mass_boundary_outer(data);
	break;
    }
    case parameters::boundary_condition_zero_gradient:
	zero_gradient_boundary_outer(data);
	break;
    case parameters::boundary_condition_boundary_layer:
	boundary_layer_outer_boundary(data);
	break;
    case parameters::boundary_condition_nonreflecting:
	NonReflectingBoundary_outer(data, &data[t_data::V_RADIAL],
				    &data[t_data::SIGMA],
				    &data[t_data::ENERGY]);
	break;
    case parameters::boundary_condition_jibin_spreading_ring:
	spreading_ring_outer(data);
	break;
    case parameters::boundary_condition_precribed_time_variable: {
	boundary_condition_precribed_time_variable_outer(data,
							 &data[t_data::SIGMA], current_time);
    } break;
    case parameters::boundary_condition_viscous_outflow:
	die("outer viscous outflow boundary not implemented");
	break;
    case parameters::boundary_condition_evanescent:
	// EvanescentBoundary_outer(VRadial, VAzimuthal, Density, Energy, dt,
	// sys);
	break;
    case parameters::boundary_condition_keplerian:
	keplerian2d_boundary_outer(data);
	break;
    }

    /// d Omega / dr = 0 has really bad effects on the massflow test
    /// not recommended
    if (parameters::domegadr_zero) {
		zero_shear_boundary(data[t_data::V_RADIAL]);
    }

    if (OuterSourceMass) {
	ApplyOuterSourceMass(&data[t_data::SIGMA], &data[t_data::V_RADIAL]);
    }

    if (parameters::massoverflow) {
	boundary_conditions::mass_overflow_willy(data, nullptr, false);
    }
}

} // namespace boundary_conditions
