/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"
#include "../parameters.h"
#include "../find_cell_id.h"
#include "../global.h"
#include "../Theo.h"
#include "../util.h"
#include "../SourceEuler.h"
#include "../logging.h"
#include "../LowTasks.h"

#include <vector>
#include <cassert>

namespace boundary_conditions
{

std::vector<t_DampingType> damping_vector;

bool damping_enabled;
bool is_damping_reference = false;
double damping_inner_limit;
double damping_outer_limit;
double damping_time_factor;
double damping_time_radius_outer;

int damping_energy_id;

static t_DampingType write_damping_type(t_damping_type type_inner,
					t_damping_type type_outer,
					t_data::t_polargrid_type quantity,
					t_data::t_polargrid_type quantity0,
					std::string description)
{
    t_DampingType damping_type;
    damping_type.array_to_damp = quantity;
    damping_type.array_with_damping_values = quantity0;
    damping_type.type_inner = type_inner;
    damping_type.type_outer = type_outer;
    std::string description_inner;
    std::string description_outer;

    switch (type_inner) {
    case damping_none:
	damping_type.inner_damping_function = nullptr;
	description_inner =
	    "Damping " + description + " is disabled at inner boundary.";
	break;
    case damping_mean:
	damping_type.inner_damping_function =
	    &boundary_conditions::damping_single_inner_mean;
	description_inner =
	    "Damping " + description + " to mean value at inner boundary.";
	break;
    case damping_reference:
	damping_type.inner_damping_function =
	    &boundary_conditions::damping_single_inner;
	description_inner =
	    "Damping " + description + " to reference value at inner boundary.";
	break;
    case damping_zero:
	damping_type.inner_damping_function =
	    &boundary_conditions::damping_single_inner_zero;
	description_inner =
	    "Damping " + description + " to zero at inner boundary.";
	break;
    case damping_visc:
	damping_type.inner_damping_function =
	    &boundary_conditions::damping_vradial_inner_visc;
	description_inner = "Damping " + description +
			    " to viscous radial speed at inner boundary.";
	break;
    }

    switch (type_outer) {
    case damping_none:
	damping_type.outer_damping_function = nullptr;
	description_outer =
	    "Damping " + description + " is disabled at outer boundary.";
	break;
    case damping_mean:
	damping_type.outer_damping_function =
	    &boundary_conditions::damping_single_outer_mean;
	description_outer =
	    "Damping " + description + " to mean value at outer boundary.";
	break;
    case damping_reference:
	damping_type.outer_damping_function =
	    &boundary_conditions::damping_single_outer;
	description_outer =
	    "Damping " + description + " to reference value at outer boundary.";
	break;
    case damping_zero:
	damping_type.outer_damping_function =
	    &boundary_conditions::damping_single_outer_zero;
	description_outer =
	    "Damping " + description + " to zero at outer boundary.";
	break;
    case damping_visc:
	die(("Damping " + description +
	     " to viscous radial speed at outer boundary not implemented!\n")
		.c_str());
	break;
    }

    logging::print_master(LOG_INFO "%s\n", description_inner.c_str());
    logging::print_master(LOG_INFO "%s\n", description_outer.c_str());

    return damping_type;
}

/**
	Get a value as t_damping_type to a corresponding key  if
   available, else set to default

	\param key key
	\param defvalue default value
	\returns t_damping_type
*/
static t_damping_type value_as_boudary_damping_default(const char *key,
						const char *defvalue)
{
    const std::string ret = config::cfg.get<std::string>(key, defvalue);

    t_damping_type boundary_condition;
    switch (tolower(ret[0])) {
    case 'n':
	boundary_condition = t_damping_type::damping_none;
	break;
	case 'r':
	boundary_condition = t_damping_type::damping_reference;
	break;
    case 'i':
	boundary_condition = t_damping_type::damping_reference;
	break;
    case 'y': // for legacy compatibility
	boundary_condition = t_damping_type::damping_reference;
	break;
    case 'm':
	boundary_condition = t_damping_type::damping_mean;
	break;
    case 'z':
	boundary_condition = t_damping_type::damping_zero;
	break;
    case 'v':
	boundary_condition = t_damping_type::damping_visc;
	break;
    default:
	boundary_condition = t_damping_type::damping_none;
    }
    return boundary_condition;

}

void damping_config() {

    damping_enabled = config::cfg.get_flag("Damping", false);

    damping_inner_limit = config::cfg.get<double>("DampingInnerLimit", 1.05);
    if (damping_inner_limit < 1) {
	die("DampingInnerLimit must not be <1\n");
    }
    damping_outer_limit = config::cfg.get<double>("DampingOuterLimit", 0.95);
    if (damping_outer_limit > 1) {
	die("DampingOuterLimit must not be >1\n");
    }
    damping_time_factor = config::cfg.get<double>("DampingTimeFactor", 1.0);

    damping_time_radius_outer =
	config::cfg.get<double>("DampingTimeRadiusOuter", RMAX);

    logging::print_master(
	"DampingTimeFactor: %.5e Outer damping time is computed at radius of %.5e\n",
	damping_time_factor, damping_time_radius_outer);

    t_damping_type tmp_damping_inner;
    t_damping_type tmp_damping_outer;

    if (config::cfg.contains("DampingVRadial"))
	die("DampingVRadial flag is decrepated used DampingVRadialInner and DampingVRadialOuter instead!");

    tmp_damping_inner =
	value_as_boudary_damping_default("DampingVRadialInner", "None");
    tmp_damping_outer =
	value_as_boudary_damping_default("DampingVRadialOuter", "None");

    if (tmp_damping_inner == t_damping_type::damping_visc) {
	damping_vector.push_back(
	    write_damping_type(tmp_damping_inner, tmp_damping_outer,
			       t_data::V_RADIAL, t_data::VISCOSITY, "VRadial"));
    } else {
	damping_vector.push_back(
	    write_damping_type(tmp_damping_inner, tmp_damping_outer,
			       t_data::V_RADIAL, t_data::V_RADIAL0, "VRadial"));
    }

    if (config::cfg.contains("DampingVAzimuthal"))
	die("DampingVRadial flag is decrepated used DampingVAzimuthalInner and DampingVAzimuthalOuter instead!");

    tmp_damping_inner =
	value_as_boudary_damping_default("DampingVAzimuthalInner", "None");
    tmp_damping_outer =
	value_as_boudary_damping_default("DampingVAzimuthalOuter", "None");

    damping_vector.push_back(write_damping_type(
	tmp_damping_inner, tmp_damping_outer, t_data::V_AZIMUTHAL,
	t_data::V_AZIMUTHAL0, "VAzimuthal"));

    if (config::cfg.contains("DampingSurfaceDensity"))
	die("DampingSurfaceDensity flag is decrepated used DampingSurfaceDensityInner and DampingSurfaceDensityOuter instead!");

    tmp_damping_inner =
	value_as_boudary_damping_default("DampingSurfaceDensityInner", "None");
    tmp_damping_outer =
	value_as_boudary_damping_default("DampingSurfaceDensityOuter", "None");

    damping_vector.push_back(
	write_damping_type(tmp_damping_inner, tmp_damping_outer, t_data::SIGMA,
			   t_data::SIGMA0, "SurfaceDensity"));

    if (config::cfg.contains("DampingEnergy"))
	die("DampingEnergy flag is decrepated used DampingEnergyInner and DampingEnergyOuter instead!");

    tmp_damping_inner =
	value_as_boudary_damping_default("DampingEnergyInner", "None");
    tmp_damping_outer =
	value_as_boudary_damping_default("DampingEnergyOuter", "None");

    damping_vector.push_back(
	write_damping_type(tmp_damping_inner, tmp_damping_outer, t_data::ENERGY,
			   t_data::ENERGY0, "Energy"));
    damping_energy_id = (int)damping_vector.size() - 1;
}

void describe_damping() {
	if (!damping_enabled) {
		logging::print_master(LOG_INFO "Damping at boundaries is disabled.\n");
		is_damping_reference = false;
    } else {
		is_damping_reference = false;
		for (unsigned int i = 0; i < damping_vector.size(); ++i) {
			is_damping_reference =
			is_damping_reference ||
			(damping_vector[i].type_inner == damping_reference);
			is_damping_reference =
			is_damping_reference ||
			(damping_vector[i].type_outer == damping_reference);
		}
    }
}



// Determine whether initial values must be stored.
bool initial_values_needed() {
	const bool damping = is_damping_reference;
	const bool betacooling = parameters::cooling_beta_reference;
	const bool boundary = boundary_conditions::reference_values_needed();
	return damping || betacooling || boundary;
}


void copy_initial_values(t_data &data) {
    if (initial_values_needed()) {
	// save starting values (needed for damping)
		copy_polargrid(data[t_data::V_RADIAL0], data[t_data::V_RADIAL]);
		copy_polargrid(data[t_data::V_AZIMUTHAL0], data[t_data::V_AZIMUTHAL]);
		copy_polargrid(data[t_data::SIGMA0], data[t_data::SIGMA]);
		copy_polargrid(data[t_data::ENERGY0], data[t_data::ENERGY]);
    }
}

void damping_single_inner(t_polargrid &quantity, t_polargrid &quantity0,
			  double dt)
{
    // use the correct radius array corresponding to quantity
    t_radialarray &radius = quantity.is_scalar() ? Rb : Ra;
    const bool is_density = strcmp(quantity.get_name(), "Sigma") == 0;

    // is this CPU in the inner damping domain?
    if ((damping_inner_limit > 1.0) &&
	(radius[0] < RMIN * damping_inner_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMIN * damping_inner_limit), false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMIN * damping_inner_limit), true);
	}

    const double tau = damping_time_factor * 2.0 * M_PI /
		     calculate_omega_kepler(RMIN);

	// Needed for OpenMP to work
	double &InnerWaveDampingMassCreation = MassDelta.InnerWaveDampingMassCreation;
	double &InnerWaveDampingMassRemoval = MassDelta.InnerWaveDampingMassRemoval;

	#pragma omp parallel for reduction (+ : InnerWaveDampingMassCreation, InnerWaveDampingMassRemoval)
	for (unsigned int n_radial = 0; n_radial <= limit; ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMIN * damping_inner_limit) /
		    (RMIN - RMIN * damping_inner_limit),
		2);
	    double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
        const double X = quantity(n_radial, n_azimuthal);
        const double X0 = quantity0(n_radial, n_azimuthal);
        const double Xnew = (X - X0) * exp_factor + X0;
		quantity(n_radial, n_azimuthal) = Xnew;
        const double delta = Xnew - X;
		if (is_density) {
		    if (delta > 0) {
			sum_without_ghost_cells(InnerWaveDampingMassCreation,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(InnerWaveDampingMassRemoval,
						-delta * Surf[n_radial],
						n_radial);
		    }
		}
	    }
	}
    }
}

void damping_single_outer(t_polargrid &quantity, t_polargrid &quantity0,
			  double dt)
{
    // use the correct radius array corresponding to quantity
    t_radialarray &radius = quantity.is_scalar() ? Rb : Ra;
    const bool is_density = strcmp(quantity.get_name(), "Sigma") == 0;

    // is this CPU in the outer damping domain?
    if ((damping_outer_limit < 1.0) &&
	(radius[quantity.get_max_radial()] >
	 RMAX * damping_outer_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMAX * damping_outer_limit) + 1, false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMAX * damping_outer_limit) + 1, true);
	}

	double tau = damping_time_factor * 2.0 * M_PI /
			 calculate_omega_kepler(damping_time_radius_outer);

	// Needed for OpenMP to work
	double &OuterWaveDampingMassCreation = MassDelta.OuterWaveDampingMassCreation;
	double &OuterWaveDampingMassRemoval = MassDelta.OuterWaveDampingMassRemoval;

	#pragma omp parallel for reduction (+ : OuterWaveDampingMassCreation, OuterWaveDampingMassRemoval)
	for (unsigned int n_radial = limit;
	     n_radial < quantity.get_size_radial(); ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMAX * damping_outer_limit) /
		    (RMAX - RMAX * damping_outer_limit),
		2);
	    double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
        const double X = quantity(n_radial, n_azimuthal);
        const double X0 = quantity0(n_radial, n_azimuthal);
        const double Xnew = (X - X0) * exp_factor + X0;
		quantity(n_radial, n_azimuthal) = Xnew;
        const double delta = Xnew - X;
		if (is_density) {
		    if (delta > 0) {
			sum_without_ghost_cells(OuterWaveDampingMassCreation,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(OuterWaveDampingMassRemoval,
						-delta * Surf[n_radial],
						n_radial);
		    }
		}
	    }
	}
    }
}

void damping_single_inner_zero(t_polargrid &quantity, t_polargrid &quantity0,
			       double dt)
{
    (void)quantity0;
    // use the correct radius array corresponding to quantity
    t_radialarray &radius = quantity.is_scalar() ? Rb : Ra;
    const bool is_density = strcmp(quantity.get_name(), "Sigma") == 0;

    // is this CPU in the inner damping domain?
    if ((damping_inner_limit > 1.0) &&
	(radius[0] < RMIN * damping_inner_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMIN * damping_inner_limit), false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMIN * damping_inner_limit), true);
	}

    const double tau = damping_time_factor * 2.0 * M_PI /
		     calculate_omega_kepler(RMIN);

	// Needed for OpenMP to work
	double &InnerWaveDampingMassCreation = MassDelta.InnerWaveDampingMassCreation;
	double &InnerWaveDampingMassRemoval = MassDelta.InnerWaveDampingMassRemoval;

	#pragma omp parallel for reduction (+ : InnerWaveDampingMassCreation, InnerWaveDampingMassRemoval)
	for (unsigned int n_radial = 0; n_radial <= limit; ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMIN * damping_inner_limit) /
		    (RMIN - RMIN * damping_inner_limit),
		2);
	    double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
		const double X = quantity(n_radial, n_azimuthal);
		double X0;
		if (is_density) {
		    X0 = parameters::sigma_floor * parameters::sigma0;
		} else {
		    X0 = 0.0;
		}
		const double Xnew = (X - X0) * exp_factor + X0;
		quantity(n_radial, n_azimuthal) = Xnew;
		const double delta = Xnew - X;
		if (is_density) {
		    if (delta > 0) {
			sum_without_ghost_cells(InnerWaveDampingMassCreation,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(InnerWaveDampingMassRemoval,
						-delta * Surf[n_radial],
						n_radial);
		    }
		}
	    }
	}
    }
}

void damping_single_outer_zero(t_polargrid &quantity, t_polargrid &quantity0,
			       double dt)
{
    (void)quantity0;
    // use the correct radius array corresponding to quantity
    t_radialarray &radius = quantity.is_scalar() ? Rb : Ra;

    const bool is_density = strcmp(quantity.get_name(), "Sigma") == 0;

    // is this CPU in the outer damping domain?
    if ((damping_outer_limit < 1.0) &&
	(radius[quantity.get_max_radial()] >
	 RMAX * damping_outer_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMAX * damping_outer_limit) + 1, false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMAX * damping_outer_limit) + 1, true);
	}

    const double tau = damping_time_factor * 2.0 * M_PI /
			 calculate_omega_kepler(damping_time_radius_outer);

	// Needed for OpenMP to work
	double &OuterWaveDampingMassCreation = MassDelta.OuterWaveDampingMassCreation;
	double &OuterWaveDampingMassRemoval = MassDelta.OuterWaveDampingMassRemoval;

	#pragma omp parallel for reduction (+ : OuterWaveDampingMassCreation, OuterWaveDampingMassRemoval)
	for (unsigned int n_radial = limit;
	     n_radial < quantity.get_size_radial(); ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMAX * damping_outer_limit) /
		    (RMAX - RMAX * damping_outer_limit),
		2);
	    double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
		const double X = quantity(n_radial, n_azimuthal);
		const double X0 = 0.0;
		const double Xnew = (X - X0) * exp_factor + X0;
		quantity(n_radial, n_azimuthal) = Xnew;
		const double delta = Xnew - X;
		if (is_density) {
		    if (delta > 0) {
			sum_without_ghost_cells(OuterWaveDampingMassCreation,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(OuterWaveDampingMassRemoval,
						-delta * Surf[n_radial],
						n_radial);
		    }
		}
	    }
	}
    }
}

void damping_single_inner_mean(t_polargrid &quantity, t_polargrid &quantity0,
			       double dt)
{
    // use the correct radius array corresponding to quantity
    t_radialarray &radius = quantity.is_scalar() ? Rb : Ra;
    const bool is_density = strcmp(quantity.get_name(), "Sigma") == 0;

    // is this CPU in the inner damping domain?
    if ((damping_inner_limit > 1.0) &&
	(radius[0] < RMIN * damping_inner_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMIN * damping_inner_limit), false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMIN * damping_inner_limit), true);
	}

    const double tau = damping_time_factor * 2.0 * M_PI /
		     calculate_omega_kepler(RMIN);

	// get mean quantity
	for (unsigned int n_radial = 0; n_radial <= limit; ++n_radial) {
	    quantity0(n_radial, 0) = 0.0;
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
		quantity0(n_radial, 0) += quantity(n_radial, n_azimuthal);
	    }
	    quantity0(n_radial, 0) /= quantity.get_size_azimuthal();
	}

	// Needed for OpenMP to work, OpenMP want a local variable for reduction
	double &InnerWaveDampingMassCreation = MassDelta.InnerWaveDampingMassCreation;
	double &InnerWaveDampingMassRemoval = MassDelta.InnerWaveDampingMassRemoval;

	#pragma omp parallel for reduction (+ : InnerWaveDampingMassCreation, InnerWaveDampingMassRemoval)
	for (unsigned int n_radial = 0; n_radial <= limit; ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMIN * damping_inner_limit) /
		    (RMIN - RMIN * damping_inner_limit),
		2);
        const double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
        const double X = quantity(n_radial, n_azimuthal);
        const double X0 = quantity0(n_radial, 0);
        const double Xnew = (X - X0) * exp_factor + X0;
		quantity(n_radial, n_azimuthal) = Xnew;
        const double delta = Xnew - X;
		if (is_density) {
		    if (delta > 0) {
			sum_without_ghost_cells(InnerWaveDampingMassCreation,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(InnerWaveDampingMassRemoval,
						-delta * Surf[n_radial],
						n_radial);
		    }
		}
	    }
	}
    }
}

void damping_vradial_inner_visc(t_polargrid &vrad, t_polargrid &viscosity,
				double dt)
{
    bool is_vrad = strcmp(vrad.get_name(), "vrad") == 0;
    assert(is_vrad);
    (void)is_vrad; /// removes IDE 'unused' warning

    t_radialarray &rinf = Ra;

    // is this CPU in the inner damping domain?
    if ((damping_inner_limit > 1.0) &&
	(rinf[0] < RMIN * damping_inner_limit)) {
	// find range

	const static unsigned int limit =
	    get_rinf_id(RMIN * damping_inner_limit);

	const double tau = damping_time_factor * 2.0 * M_PI /
			   calculate_omega_kepler(RMIN);

	const double s = viscous_outflow_speed;

	#pragma omp parallel for
	for (unsigned int n_radial = Zero_no_ghost; n_radial <= limit;
	     ++n_radial) {
	    double factor = std::pow(
		(rinf[n_radial] - RMIN * damping_inner_limit) /
		    (RMIN - RMIN * damping_inner_limit),
		2);
        const double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < vrad.get_size_azimuthal(); ++n_azimuthal) {
		const double viscosity_above = viscosity(n_radial, n_azimuthal);

		double Nu;

		if (n_radial == 0) {
		    Nu = viscosity_above;
		} else {
		    const double viscosity_below =
			viscosity(n_radial - 1, n_azimuthal);
		    Nu = 0.5 * (viscosity_above + viscosity_below);
		}

		// V_rad =  - 1.5 / r * Nu (Kley, Papaloizou and Ogilvie, 2008)
		const double V_rad = -1.5 * s / Rinf[n_radial] * Nu;

		const double X = vrad(n_radial, n_azimuthal);
		const double X0 = V_rad;
		const double Xnew = (X - X0) * exp_factor + X0;
		vrad(n_radial, n_azimuthal) = Xnew;
	    }
	}
    }
}

void damping_single_outer_mean(t_polargrid &quantity, t_polargrid &quantity0,
			       double dt)
{
    // use the correct radius array corresponding to quantity
    t_radialarray &radius = quantity.is_scalar() ? Rb : Ra;
    bool is_density = strcmp(quantity.get_name(), "Sigma") == 0;

    // is this CPU in the outer damping domain?
    if ((damping_outer_limit < 1.0) &&
	(radius[quantity.get_max_radial()] >
	 RMAX * damping_outer_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMAX * damping_outer_limit) + 1, false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMAX * damping_outer_limit) + 1, true);
	}

    const double tau = damping_time_factor * 2.0 * M_PI /
			 calculate_omega_kepler(damping_time_radius_outer);

	for (unsigned int n_radial = limit;
	     n_radial < quantity.get_size_radial(); ++n_radial) {
	    quantity0(n_radial, 0) = 0.0;
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
		quantity0(n_radial, 0) += quantity(n_radial, n_azimuthal);
	    }
	    quantity0(n_radial, 0) /= quantity.get_size_azimuthal();
	}

	// Needed for OpenMP to work
	double &OuterWaveDampingMassCreation = MassDelta.OuterWaveDampingMassCreation;
	double &OuterWaveDampingMassRemoval = MassDelta.OuterWaveDampingMassRemoval;

	#pragma omp parallel for reduction (+ : OuterWaveDampingMassCreation, OuterWaveDampingMassRemoval)
	for (unsigned int n_radial = limit;
	     n_radial < quantity.get_size_radial(); ++n_radial) {
        const double factor = std::pow(
		(radius[n_radial] - RMAX * damping_outer_limit) /
		    (RMAX - RMAX * damping_outer_limit),
		2);
        const double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
        const double X = quantity(n_radial, n_azimuthal);
        const double X0 = quantity0(n_radial, 0);
        const double Xnew = (X - X0) * exp_factor + X0;
		quantity(n_radial, n_azimuthal) = Xnew;
        const double delta = Xnew - X;
		if (is_density) {
		    if (delta > 0) {
			sum_without_ghost_cells(OuterWaveDampingMassCreation,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(OuterWaveDampingMassRemoval,
						-delta * Surf[n_radial],
						n_radial);
		    }
		}
	    }
	}
    }
}


/**
	damping of all selected quantities
*/
void damping(t_data &data, const double dt)
{
	// Iterate over all function pointers previously stored
    if (damping_enabled) {
	const unsigned int number_of_quantities_to_damp =
	    damping_vector.size();

	for (unsigned int i = 0; i < number_of_quantities_to_damp; ++i) {
	    const t_DampingType *damper =
		&damping_vector[i];
	    if (damper->inner_damping_function != nullptr)
		(damper->inner_damping_function)(
		    data[damper->array_to_damp],
		    data[damper->array_with_damping_values], dt);
	    if (damper->outer_damping_function != nullptr)
		(damper->outer_damping_function)(
		    data[damper->array_to_damp],
		    data[damper->array_with_damping_values], dt);
	}
    }
}

} // namespace boundary_conditions
