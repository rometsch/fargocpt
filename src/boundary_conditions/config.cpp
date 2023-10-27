/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"

#include "../config.h"
#include "../global.h"
#include "../logging.h"
#include "../units.h"

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
std::string special_name = "";
std::string composite_inner_name = "";
std::string composite_outer_name = "";

double keplerian_azimuthal_outer_factor = 1.0;
double keplerian_azimuthal_inner_factor = 1.0;
double keplerian_radial_outer_factor = 1.0;
double keplerian_radial_inner_factor = 1.0;


std::string PRESCRIBED_BOUNDARY_OUTER_FILE;

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
 * Get the type of an individual variabl;es boundary condition for a given key.
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
	    throw std::runtime_error(
		"Can not set '" + key +
		"' via composite when type is different from 'infer'");
	} else {
	    name = str;
	    return str;
	}
    }
}

static void sigma_inner()
{
    std::string str = get_type("BoundarySigmaInner", sigma_inner_name);

    if (str == "zerogradient") {
	sigma_inner_func = zero_gradient_inner;
    } else if (str == "diskmodel") {
	sigma_inner_func = diskmodel_inner_sigma;
    } else if (str == "reference") {
	sigma_inner_func = reference_inner;
    } else {
	throw std::runtime_error(
	    "Unknown boundary condition for sigma inner: " + str);
    }

    logging::print_master(LOG_INFO "BC: Sigma inner = %s\n",
			  sigma_inner_name.c_str());
}

static void sigma_outer()
{
    const std::string str = get_type("BoundarySigmaOuter", sigma_outer_name);

    if (str == "zerogradient") {
	sigma_outer_func = zero_gradient_outer;
    } else if (str == "diskmodel") {
	sigma_outer_func = diskmodel_outer_sigma;
    } else if (str == "reference") {
	sigma_outer_func = reference_outer;
    } else {
	throw std::runtime_error(
	    "Unknown boundary condition for sigma outer: " + str);
    }
    logging::print_master(LOG_INFO "BC: Sigma outer = %s\n",
			  sigma_outer_name.c_str());
}

static void energy_inner()
{
    const std::string str = get_type("BoundaryEnergyInner", energy_outer_name);

    if (str == "zerogradient") {
	energy_inner_func = zero_gradient_inner;
    } else if (str == "diskmodel") {
	energy_inner_func = diskmodel_inner_energy;
    } else if (str == "reference") {
	energy_inner_func = reference_inner;
    } else {
	throw std::runtime_error(
	    "Unknown boundary condition for energy inner: " + str);
    }
    logging::print_master(LOG_INFO "BC: Energy inner = %s\n",
			  energy_inner_name.c_str());
}

static void energy_outer()
{
    const std::string str = get_type("BoundaryEnergyOuter", energy_outer_name);

    if (str == "zerogradient") {
	energy_outer_func = zero_gradient_outer;
    } else if (str == "diskmodel") {
	energy_outer_func = diskmodel_outer_energy;
    } else if (str == "reference") {
	energy_outer_func = reference_outer;
    } else {
	throw std::runtime_error(
	    "Unknown boundary condition for energy outer: " + str);
    }
    logging::print_master(LOG_INFO "BC: Energy outer = %s\n",
			  energy_outer_name.c_str());
}

static void vrad_inner()
{
    const std::string str = get_type("BoundaryVradInner", vrad_inner_name);

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
    } else {
	throw std::runtime_error("Unknown boundary condition for vrad inner: " +
				 str);
    }
    // parse parameters
    keplerian_radial_inner_factor =
	config::cfg.get<double>("BoundaryVradInnerKeplerianFactor", 0.1);
    // log
    logging::print_master(LOG_INFO "BC: Vrad inner = %s\n",
			  vrad_inner_name.c_str());
}

static void vrad_outer()
{
    const std::string str = get_type("BoundaryVradOuter", vrad_outer_name);

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
    } else {
	throw std::runtime_error("Unknown boundary condition for vrad outer: " +
				 str);
    }
    // parse parameters
    keplerian_radial_outer_factor =
	config::cfg.get<double>("BoundaryVradOuterKeplerianFactor", 0.1);
    // log
    logging::print_master(LOG_INFO "BC: Vrad outer = %s\n",
			  vrad_outer_name.c_str());
}

static void vaz_inner()
{
    const std::string str =
	config::cfg.get_lowercase("BoundaryVazInner", "keplerian");
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
    } else {
	throw std::runtime_error("Unknown boundary condition for vaz inner: " +
				 str);
    }
    // parse parameters
    keplerian_azimuthal_inner_factor =
	config::cfg.get<double>("BoundaryVazInnerKeplerianFactor", 1.0);
    // log
    logging::print_master(LOG_INFO "BC: Vaz inner = %s\n",
			  vaz_inner_name.c_str());
}

static void vaz_outer()
{
    const std::string str =
	config::cfg.get_lowercase("BoundaryVazOuter", "keplerian");
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
    } else {
	throw std::runtime_error("Unknown boundary condition for vaz outer: " +
				 str);
    }
    // parse parameters
    keplerian_azimuthal_outer_factor =
	config::cfg.get<double>("BoundaryVazOuterKeplerianFactor", 1.0);
    // log
    logging::print_master(LOG_INFO "BC: Vaz outer = %s\n",
			  vaz_outer_name.c_str());
}

static void special()
{
    const std::string str =
	config::cfg.get_lowercase("BoundarySpecial", "none");
    special_name = str;
    if (str == "com") {
    } else if (str == "massoverflow") {
    } else if (str == "custom") {
    } else if (str == "none") {
    } else {
	throw std::runtime_error("Unknown boundary condition for special: " +
				 str);
    }
    logging::print_master(LOG_INFO "BC: Special = %s\n", str.c_str());
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
		throw std::logic_error("TODO: Selection of centerofmass not yet implemented.");
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
		throw std::logic_error("TODO: Selection of centerofmass not yet implemented.");
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

    special();

    // check if file for prescribed time variable boundary exists
    if (config::cfg.contains("PRESCRIBEDBOUNDARYFILEOUTER")) {
	if (config::cfg.get<std::string>("PRESCRIBEDBOUNDARYFILEOUTER")
		.length() > 0) {
	    PRESCRIBED_BOUNDARY_OUTER_FILE =
		config::cfg.get<std::string>("PRESCRIBEDBOUNDARYFILEOUTER");
	} else {
	    die("Error looking for data for the prescribed time variable boundary condition. Path could not be read!\n");
	}
    } else {
	PRESCRIBED_BOUNDARY_OUTER_FILE = "";
    }

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
