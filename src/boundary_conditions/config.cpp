/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"

#include "../config.h"
#include "../logging.h"
#include "../units.h"
#include "../LowTasks.h"

#include <string>

extern boolean OuterSourceMass;

namespace boundary_conditions
{

void (*sigma_inner_func)(t_polargrid &, t_polargrid &, t_data &);
void (*sigma_outer_func)(t_polargrid &, t_polargrid &, t_data &);
void (*energy_inner_func)(t_polargrid &, t_polargrid &, t_data &);
void (*energy_outer_func)(t_polargrid &, t_polargrid &, t_data &);
void (*vrad_inner_func)(t_polargrid &, t_polargrid &, t_data &);
void (*vrad_outer_func)(t_polargrid &, t_polargrid &, t_data &);
void (*vaz_inner_func)(t_polargrid &, t_polargrid &, t_data &);
void (*vaz_outer_func)(t_polargrid &, t_polargrid &, t_data &);
void (*special_inner_func)(t_polargrid &, t_polargrid &, t_data &);
void (*special_outer_func)(t_polargrid &, t_polargrid &, t_data &);

std::string sigma_inner_name = "";
std::string sigma_outer_name = "";
std::string energy_inner_name = "";
std::string energy_outer_name = "";
std::string vrad_inner_name = "";
std::string vrad_outer_name = "";
std::string vaz_inner_name = "";
std::string vaz_outer_name = "";
std::string composite_inner_name = "";
std::string composite_outer_name = "";

double keplerian_azimuthal_outer_factor = 1.0;
double keplerian_azimuthal_inner_factor = 1.0;
double keplerian_radial_outer_factor = 1.0;
double keplerian_radial_inner_factor = 1.0;


bool domegadr_zero;

double viscous_outflow_speed;
bool massoverflow;
bool mof_variableTransfer;
unsigned int mof_planet;
double mof_temperature;
double mof_value;
double mof_rampingtime;
double mof_averaging_time;
double mof_gamma;


/* 
 * Dummy to call when all variables are handled by one special function.
*/



/*
 * Get the type of an individual variable's boundary condition for a given key.
 * This function takes into account types set by a composite boundary condition.
 * See function for composites below.
 */
static std::string get_type(const std::string key, std::string &name)
{
    std::string str = config::cfg.get_lowercase(key, "infer");
    if (str == "infer") {
	if (name == "") { // nothing to infer here
	    throw std::runtime_error(
		"Can not infer '" + key +
		"' when 'InnerBoundary/OuterBoundary' is set to 'individual'");
	} else {
	    return name;
	}
    } else { // set individual boundary directly
	if (name != "") {
	    logging::print_master(LOG_INFO "BC: Overwriting %s from %s set via composite to %s.\n",
			  key.c_str(), name.c_str(), str.c_str());
	}
	name = str;
	return str;
    }
}

static void sigma_inner()
{
    std::string str = get_type("InnerBoundarySigma", sigma_inner_name);

    if (str == "zerogradient") {
	sigma_inner_func = zero_gradient_inner;
    } else if (str == "diskmodel") {
	sigma_inner_func = diskmodel_inner_sigma;
    } else if (str == "reference") {
	sigma_inner_func = reference_inner;
    } else if (str == "none") {
	sigma_inner_func = no_operation;
    if (composite_inner_name == "") {
        throw std::runtime_error(
            "Can not set InnerBoundarySigma to 'none' without selecting special boundary (e.g. 'custom').");
    }
    } else {
	throw std::runtime_error(
	    "Unknown boundary condition for sigma inner: " + str);
    }

    logging::print_master(LOG_INFO "BC: Sigma inner = %s\n",
			  sigma_inner_name.c_str());
}

static void sigma_outer()
{
    const std::string str = get_type("OuterBoundarySigma", sigma_outer_name);

    if (str == "zerogradient") {
	sigma_outer_func = zero_gradient_outer;
    } else if (str == "diskmodel") {
	sigma_outer_func = diskmodel_outer_sigma;
    } else if (str == "reference") {
	sigma_outer_func = reference_outer;
    } else if (str == "none") {
    sigma_outer_func = no_operation;
    if (composite_outer_name == "") {
        throw std::runtime_error(
            "Can not set OuterBoundarySigma to 'none' without selecting special boundary (e.g. 'custom').");
    }
    } else {
	throw std::runtime_error(
	    "Unknown boundary condition for sigma outer: " + str);
    }
    logging::print_master(LOG_INFO "BC: Sigma outer = %s\n",
			  sigma_outer_name.c_str());
}

static void energy_inner()
{
    const std::string str = get_type("InnerBoundaryEnergy", energy_outer_name);

    if (str == "zerogradient") {
	energy_inner_func = zero_gradient_inner;
    } else if (str == "diskmodel") {
	energy_inner_func = diskmodel_inner_energy;
    } else if (str == "reference") {
	energy_inner_func = reference_inner;
    } else if (str == "none") {
    energy_inner_func = no_operation;
    if (composite_inner_name == "") {
        throw std::runtime_error(
            "Can not set InnerBoundaryEnergy to 'none' without selecting special boundary (e.g. 'custom').");
    }
    } else {
	throw std::runtime_error(
	    "Unknown boundary condition for energy inner: " + str);
    }
    logging::print_master(LOG_INFO "BC: Energy inner = %s\n",
			  energy_inner_name.c_str());
}

static void energy_outer()
{
    const std::string str = get_type("OuterBoundaryEnergy", energy_outer_name);

    if (str == "zerogradient") {
	energy_outer_func = zero_gradient_outer;
    } else if (str == "diskmodel") {
	energy_outer_func = diskmodel_outer_energy;
    } else if (str == "reference") {
	energy_outer_func = reference_outer;
    } else if (str == "none") {
    energy_outer_func = no_operation;
    if (composite_outer_name == "") {
        throw std::runtime_error(
            "Can not set OuterBoundaryEnergy to 'none' without selecting special boundary (e.g. 'custom').");
    }
    } else {
	throw std::runtime_error(
	    "Unknown boundary condition for energy outer: " + str);
    }
    logging::print_master(LOG_INFO "BC: Energy outer = %s\n",
			  energy_outer_name.c_str());
}

static void vrad_inner()
{
    const std::string str = get_type("InnerBoundaryVrad", vrad_inner_name);

    if (str == "zerogradient") {
	vrad_inner_func = zero_gradient_inner;
    } else if (str == "reference") {
	vrad_inner_func = reference_inner;
    } else if (str == "reflecting") {
	vrad_inner_func = reflecting_inner;
    } else if (str == "outflow") {
	vrad_inner_func = outflow_inner;
    } else if (str == "viscous") {
	vrad_inner_func = viscous_outflow_inner;
    } else if (str == "keplerian") {
	vrad_inner_func = keplerian_radial_inner;
    } else if (str == "none") {
    vrad_inner_func = no_operation;
    if (composite_inner_name == "") {
        throw std::runtime_error(
            "Can not set InnerBoundaryVrad to 'none' without selecting special boundary (e.g. 'custom').");
    }
    } else {
	throw std::runtime_error("Unknown boundary condition for vrad inner: " +
				 str);
    }
    // parse parameters
    keplerian_radial_inner_factor =
	config::cfg.get<double>("InnerBoundaryVradKeplerianFactor", 0.1);
    // log
    logging::print_master(LOG_INFO "BC: Vrad inner = %s\n",
			  vrad_inner_name.c_str());
}

static void vrad_outer()
{
    const std::string str = get_type("OuterBoundaryVrad", vrad_outer_name);

    if (str == "zerogradient") {
	vrad_outer_func = zero_gradient_outer;
    } else if (str == "reference") {
	vrad_outer_func = reference_outer;
    } else if (str == "reflecting") {
	vrad_outer_func = reflecting_outer;
    } else if (str == "outflow") {
	vrad_outer_func = outflow_outer;
    } else if (str == "viscous") {
	vrad_outer_func = viscous_inflow_outer;
    } else if (str == "keplerian") {
	vrad_outer_func = keplerian_radial_outer;
    } else if (str == "none") {
    vrad_outer_func = no_operation;
    if (composite_outer_name == "") {
        throw std::runtime_error(
            "Can not set OuterBoundaryVrad to 'none' without selecting special boundary (e.g. 'custom').");
    }
    } else {
	throw std::runtime_error("Unknown boundary condition for vrad outer: " +
				 str);
    }
    // parse parameters
    keplerian_radial_outer_factor =
	config::cfg.get<double>("OuterBoundaryVradKeplerianFactor", 0.1);
    // log
    logging::print_master(LOG_INFO "BC: Vrad outer = %s\n",
			  vrad_outer_name.c_str());
}

static void vaz_inner()
{
    const std::string str =
	config::cfg.get_lowercase("InnerBoundaryVazi", "keplerian");
    vaz_inner_name = str;
    if (str == "zerogradient") {
	vaz_inner_func = zero_gradient_inner;
    } else if (str == "reference") {
	vaz_inner_func = reference_inner;
    } else if (str == "zeroshear") {
	vaz_inner_func = zero_shear_inner;
    } else if (str == "balanced") {
	vaz_inner_func = balanced_inner;
    } else if (str == "keplerian") {
	vaz_inner_func = keplerian_azimuthal_inner;
    } else if (str == "none") {
    vaz_inner_func = no_operation;
    if (composite_inner_name == "") {
        throw std::runtime_error(
            "Can not set InnerBoundaryVazi to 'none' without selecting special boundary (e.g. 'custom').");
    }
    } else {
	throw std::runtime_error("Unknown boundary condition for vaz inner: " +
				 str);
    }
    // parse parameters
    keplerian_azimuthal_inner_factor =
	config::cfg.get<double>("InnerBoundaryVaziKeplerianFactor", 1.0);
    // log
    logging::print_master(LOG_INFO "BC: Vaz inner = %s\n",
			  vaz_inner_name.c_str());
}

static void vaz_outer()
{
    const std::string str =
	config::cfg.get_lowercase("OuterBoundaryVazi", "keplerian");
    vaz_outer_name = str;
    if (str == "zerogradient") {
	vaz_outer_func = zero_gradient_outer;
    } else if (str == "reference") {
	vaz_outer_func = reference_outer;
    } else if (str == "zeroshear") {
	vaz_outer_func = zero_shear_outer;
    } else if (str == "balanced") {
	vaz_outer_func = balanced_outer;
    } else if (str == "keplerian") {
	vaz_outer_func = keplerian_azimuthal_outer;
    } else if (str == "none") {
    vaz_outer_func = no_operation;
    if (composite_outer_name == "") {
        throw std::runtime_error(
            "Can not set OuterBoundaryVazi to 'none' without selecting special boundary (e.g. 'custom').");
    }
    } else {
	throw std::runtime_error("Unknown boundary condition for vaz outer: " +
				 str);
    }
    // parse parameters
    keplerian_azimuthal_outer_factor =
	config::cfg.get<double>("OuterBoundaryVaziKeplerianFactor", 1.0);
    // log
    logging::print_master(LOG_INFO "BC: Vaz outer = %s\n",
			  vaz_outer_name.c_str());
}


static void composite_inner()
{
    const std::string str =
	config::cfg.get_lowercase("InnerBoundary", "individual");
    composite_inner_name = str;
    if (str == "individual") {
	return;
    } else if (str == "zerogradient") {
	sigma_inner_name = "zerogradient";
	energy_inner_name = "zerogradient";
	vrad_inner_name = "zerogradient";
    } else if (str == "outflow") {
	sigma_inner_name = "zerogradient";
	energy_inner_name = "zerogradient";
	vrad_inner_name = "outflow";
    } else if (str == "reflecting") {
	sigma_inner_name = "zerogradient";
	energy_inner_name = "zerogradient";
	vrad_inner_name = "reflecting";
    } else if (str == "reference") {
	sigma_inner_name = "reference";
	energy_inner_name = "reference";
	vrad_inner_name = "reference";
	} else if (str == "centerofmass") {
	sigma_inner_name = "none";
    energy_inner_name = "none";
    vrad_inner_name = "none";
    vaz_inner_name = "none";
    } else if (str == "custom") {
	sigma_inner_name = "none";
    energy_inner_name = "none";
    vrad_inner_name = "none";
    vaz_inner_name = "none";
	} else {
	throw std::runtime_error(
	    "Unknown boundary condition for inner boundary: " + str);
    }
    if (str != "individual") {
	logging::print_master(LOG_INFO "BC: Inner composite = %s\n",
			      composite_inner_name.c_str());
    }
}

static void composite_outer()
{
    const std::string str =
	config::cfg.get_lowercase("OuterBoundary", "individual");
    composite_outer_name = str;
    if (str == "individual") {
	return;
    } else if (str == "zerogradient") {
	sigma_outer_name = "zerogradient";
	energy_outer_name = "zerogradient";
	vrad_outer_name = "zerogradient";
    } else if (str == "outflow") {
	sigma_outer_name = "zerogradient";
	energy_outer_name = "zerogradient";
	vrad_outer_name = "outflow";
    } else if (str == "reflecting") {
	sigma_outer_name = "zerogradient";
	energy_outer_name = "zerogradient";
	vrad_outer_name = "reflecting";
    } else if (str == "reference") {
	sigma_outer_name = "reference";
	energy_outer_name = "reference";
	vrad_outer_name = "reference";
	} else if (str == "massoverflow") {
		throw std::logic_error("TODO: Selection of massoverflow not yet implemented.");
	} else if (str == "centerofmass") {
	sigma_outer_name = "none";
    energy_outer_name = "none";
    vrad_outer_name = "none";
    vaz_outer_name = "none";
    } else if (str == "custom") {
	sigma_outer_name = "none";
    energy_outer_name = "none";
    vrad_outer_name = "none";
    vaz_outer_name = "none";
	} else {
	throw std::runtime_error(
	    "Unknown boundary condition for outer boundary: " + str);
    }
    if (str != "individual") {
	logging::print_master(LOG_INFO "BC: Outer composite = %s\n",
			      composite_outer_name.c_str());
    }
}

void parse_config()
{

    composite_inner();
    composite_outer();

    sigma_inner();
    sigma_outer();

    energy_inner();
    energy_outer();

    vrad_inner();
    vrad_outer();

    vaz_inner();
    vaz_outer();

    domegadr_zero = config::cfg.get_flag("DomegaDrZero", false);

    if (domegadr_zero)
	die("DomegaDrZero is deprecated!");

    viscous_outflow_speed = config::cfg.get<double>("ViscousOutflowSpeed", 1.0);


    // mass overflow
    massoverflow = config::cfg.get_flag("massoverflow", "no");
	mof_variableTransfer = config::cfg.get_flag("variableTransfer", "no");
    mof_planet = config::cfg.get<int>("mofplanet", 1);
    mof_temperature = config::cfg.get<double>("moftemperature", "1000.0 K", units::Temp0);
    mof_value = config::cfg.get<double>("mofvalue", 10E-9, units::M0/units::T0);
    mof_rampingtime = config::cfg.get<double>("moframpingtime", 30.0);
	mof_averaging_time = config::cfg.get<double>("mofaveragingtime", 10.0);
	mof_gamma = config::cfg.get<double>("mofgamma", 0.5);

	damping_config();

}


} // namespace boundary_conditions
