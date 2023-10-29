/**
	\file units.cpp
	\author Tobias Mueller <Tobias_Mueller@twam.info>

	This file handles all kinds of unit information to allow output in CGS
   units
*/

#include "units.h"
#include "LowTasks.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "options.h"
#include "parameters.h"
#include "units/units.hpp"
#include "config.h"
#include "output.h"
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unordered_map>

namespace units
{

t_unit not_yet_defined;
t_unit unitless;
t_unit length;
t_unit mass;
t_unit time;
t_unit temperature;
t_unit energy;
t_unit energy_density;
t_unit density;
t_unit surface_density;
t_unit opacity;
t_unit energy_flux;
t_unit velocity;
t_unit angular_momentum;
t_unit kinematic_viscosity;
t_unit dynamic_viscosity;
t_unit acceleration;
t_unit stress;
t_unit pressure;
t_unit power;
t_unit potential;
t_unit torque;
t_unit force;
t_unit mass_accretion_rate;

llnlunits::precise_unit L0 = llnlunits::precise::cm;
llnlunits::precise_unit M0 = llnlunits::precise::g;
llnlunits::precise_unit T0 = llnlunits::precise::second;
llnlunits::precise_unit Temp0 = llnlunits::precise::Kelvin;

bool has_unit(const std::string &val){
    auto q = llnlunits::measurement_from_string(val);
    const bool ret =  (q.units().unit_type_count() > 0);
    return ret;
}

static void check_unit(llnlunits::precise_measurement q, const std::string& val, const precise_unit &unit = llnlunits::precise::one) {
	if (!q.as_unit().is_convertible(unit)) {
		die("Invalid unit transformation from '%s' to unit '%s'!\n", val.c_str(), llnlunits::to_string(unit).c_str());
	}
}

template <typename T> T parse_units(const std::string &val) {
    const T rv = val;
    return rv;
}

template <> double parse_units(const std::string &val) {
    auto q = llnlunits::measurement_from_string(val);
	check_unit(q, val);
    return (double) q.convert_to_base().value();
}

template <> unsigned int parse_units(const std::string &val) {
    auto q = llnlunits::measurement_from_string(val);
	check_unit(q, val);
    return (unsigned int) q.convert_to_base().value();
}

template <> int parse_units(const std::string &val) {
    auto q = llnlunits::measurement_from_string(val);
	check_unit(q, val);
    return (int) q.convert_to_base().value();
}

template <> double parse_units(const std::string &val, const precise_unit& unit) {
    auto q = llnlunits::measurement_from_string(val);
	check_unit(q, val, unit);
    return (double) q.value_as(unit);
}

template <> unsigned int parse_units(const std::string &val, const precise_unit& unit) {
    auto q = llnlunits::measurement_from_string(val);
	check_unit(q, val, unit);
    return (unsigned int) q.value_as(unit);
}

template <> int parse_units(const std::string &val, const precise_unit& unit) {
    auto q = llnlunits::measurement_from_string(val);
	check_unit(q, val, unit);
    return (int) q.value_as(unit);
}

template std::string parse_units<std::string>(const std::string &val);

// TODO: pick an official system
// TODO: make these visible in output
llnlunits::precise_unit solMass(1.98847e30, llnlunits::precise::kilogram);
llnlunits::precise_unit solRadius(6.95700e8, llnlunits::precise::meter);
llnlunits::precise_unit au(1.495978707e11, llnlunits::precise::meter);
llnlunits::precise_unit jupiterMass(1.8982e27, llnlunits::precise::kilogram);
llnlunits::precise_unit jupiterRadius(69911000, llnlunits::precise::meter);
llnlunits::precise_unit earthMass(5.97217e24, llnlunits::precise::kilogram);
llnlunits::precise_unit earthRadius(6371000, llnlunits::precise::meter);


static void add_astro_units() {
	llnlunits::addUserDefinedUnit("au", au);
	llnlunits::addUserDefinedUnit("solMass",   solMass);
	llnlunits::addUserDefinedUnit("solRadius", solRadius);
	llnlunits::addUserDefinedUnit("jupiterMass",   jupiterMass);
	llnlunits::addUserDefinedUnit("jupiterRadius", jupiterRadius);
	llnlunits::addUserDefinedUnit("earthMass",   earthMass);
	llnlunits::addUserDefinedUnit("earthRadius", earthRadius);
}


void set_baseunits(
	const std::string &l0s, 
	const std::string &m0s,
	const std::string &t0s,
	const std::string &temp0s
) {

	add_astro_units();

	const bool l0hu = has_unit(l0s);
	const bool m0hu = has_unit(m0s);
	const bool t0hu = has_unit(t0s);
	const bool temp0hu = has_unit(temp0s);

	if (m0hu && l0hu) {
		L0 = llnlunits::measurement_from_string(l0s).as_unit();
		M0 = llnlunits::measurement_from_string(m0s).as_unit();

		if (!L0.is_convertible(llnlunits::meter)) {
			die("Baseunit of length '%s' not convertible to meter!", l0s.c_str());
		}
		if (!M0.is_convertible(llnlunits::kilogram)) {
			die("Baseunit of mass '%s' not convertible to kilogram!", m0s.c_str());
		}

	} else if ((!m0hu) && (!l0hu)) {
		L0 = (llnlunits::measurement_from_string(l0s)*au).as_unit();
		M0 = (llnlunits::measurement_from_string(m0s)*solMass).as_unit();
		logging::print_master(LOG_INFO "Physical units implicitly applied!\n");
		logging::print_master(LOG_INFO "L0 = %f au, M0 = %f solMass\n",
			(1*L0).value_as(au),
			(1*M0).value_as(solMass));
	} else {
		die("l0 and m0 need to either all have a unit or all have no unit!\n However, they are l0 = %s, m0 = %s!", l0s.c_str(), m0s.c_str());
	}

	if (t0hu) {
		T0 = llnlunits::measurement_from_string(t0s).as_unit();	
	} else {
		T0 = llnlunits::sqrt((1 * L0 * L0 * L0) / (1* M0 * llnlunits::constants::G)).as_unit();
	}
	if (!T0.is_convertible(llnlunits::second)) {
		die("Baseunit of time '%s' not convertible to second!", t0s.c_str());
	}

	if (temp0hu) {
		Temp0 = llnlunits::measurement_from_string(temp0s).as_unit();
	} else {
		const auto G = llnlunits::constants::G;
		const auto mu = llnlunits::constants::mu;
		const auto kB = llnlunits::constants::k;
		Temp0 = (1*G*mu/kB*M0/L0).as_unit();
	}
	if (!Temp0.is_convertible(llnlunits::Kelvin)) {
		die("Baseunit of temperature '%s' not convertible to Kelvin!", temp0s.c_str());
	}
}



t_unit::t_unit()
{
	m_name = NULL;
    m_cgs_symbol = NULL;
    m_cgs_factor = 1;
}

t_unit::~t_unit() { 
	delete[] m_cgs_symbol; 
	delete[] m_name;
}

/**
	set conversion factor to cgs system

	\param factor cgs conversion factor
*/
void t_unit::set_cgs_factor(double factor)
{
    m_cgs_factor = factor;
    m_inverse_cgs_factor = 1.0 / factor;
}

/**
	set unit symbol in cgs system
*/
void t_unit::set_cgs_symbol(const char *symbol)
{
    // delete old symbol
    delete[] m_cgs_symbol;

    // aquire memory for symbol
    const unsigned int length = strlen(symbol) + 1;
    m_cgs_symbol = new char[length];

    // copy symbol
    strncpy(m_cgs_symbol, symbol, length);
}

/**
 * set name of the unit
*/
void t_unit::set_name(const char *name)
{
	// delete old name
	delete[] m_name;

	// aquire memory for name
	const unsigned int length = strlen(name) + 1;
	m_name = new char[length];

	// copy name
	strncpy(m_name, name, length);
}

/**
	get unit symbol in cgs system
*/
const char *t_unit::get_cgs_symbol(void) const { return m_cgs_symbol; }

/**
	get name of the unit
*/
const char *t_unit::get_name(void) const { return m_name; }

std::string t_unit::get_cgs_factor_symbol()
{
    // a string containing the pair of value and unit as
    // a string such as '1.7823468234...e16 g'
    // i.e. the number with format #.16e
    std::stringstream us;
    us.precision(16);
    us << std::scientific << m_cgs_factor << " " << std::string(m_cgs_symbol);
    return us.str();
}


void calculate_unit_factors()
{

    not_yet_defined.set_cgs_factor(1);
    not_yet_defined.set_cgs_symbol("1");
	not_yet_defined.set_name("not_yet_defined");

    unitless.set_cgs_factor(1);
    unitless.set_cgs_symbol("1");
	unitless.set_name("unitless");

    length.set_cgs_factor((1*L0).value_as(llnlunits::precise::cm));
    length.set_cgs_symbol("cm");
	length.set_name("length");

    mass.set_cgs_factor((1*M0).value_as(llnlunits::precise::g));
    mass.set_cgs_symbol("g");
	mass.set_name("mass");

    time.set_cgs_factor((1*T0).value_as(llnlunits::precise::second));
    time.set_cgs_symbol("s");
	time.set_name("time");

    energy.set_cgs_factor(length * length * mass / (time * time));
    energy.set_cgs_symbol("erg");
	energy.set_name("energy");

    energy_density.set_cgs_factor(mass / (time * time));
    energy_density.set_cgs_symbol("erg cm^-2");
	energy_density.set_name("energy surface density");

	temperature.set_cgs_factor((1*Temp0).value_as(llnlunits::precise::Kelvin));
    temperature.set_cgs_symbol("K");
	temperature.set_name("temperature");

    density.set_cgs_factor(mass / (length * length * length));
    density.set_cgs_symbol("g cm^-3");
	density.set_name("density");

    surface_density.set_cgs_factor(mass / (length * length));
    surface_density.set_cgs_symbol("g cm^-2");
	surface_density.set_name("mass surface density");

    opacity.set_cgs_factor(length * length / mass);
    opacity.set_cgs_symbol("g^-1 cm^2");
	opacity.set_name("opacity");

    energy_flux.set_cgs_factor(energy / (length * length * time));
    energy_flux.set_cgs_symbol("erg cm^-2 s^-1");
	energy_flux.set_name("energy flux");

    velocity.set_cgs_factor(length / time);
    velocity.set_cgs_symbol("cm s^-1");
	velocity.set_name("velocity");

    acceleration.set_cgs_factor(length / (time * time));
    acceleration.set_cgs_symbol("cm s^-2");
	acceleration.set_name("acceleration");

    angular_momentum.set_cgs_factor(length * mass * velocity);
    angular_momentum.set_cgs_symbol("cm^2 g s^-1");
	angular_momentum.set_name("angular momentum");

    kinematic_viscosity.set_cgs_factor(length * length / time);
    kinematic_viscosity.set_cgs_symbol("cm^2 s^-1");
	kinematic_viscosity.set_name("kinematic viscosity");

    dynamic_viscosity.set_cgs_factor(mass / (length * time));
    dynamic_viscosity.set_cgs_symbol("P");
	dynamic_viscosity.set_name("dynamic viscosity");

    stress.set_cgs_factor(mass / (time * time));
    stress.set_cgs_symbol("g s^-2");
	stress.set_name("stress");

    pressure.set_cgs_factor(mass / (time * time));
    pressure.set_cgs_symbol("dyn cm^-1");
	pressure.set_name("pressure");

    power.set_cgs_factor(mass * length * length / (time * time * time));
    power.set_cgs_symbol("erg/s");
	power.set_name("power");

    potential.set_cgs_factor(length * length / (time * time));
    potential.set_cgs_symbol("erg/g");
	potential.set_name("potential");

    torque.set_cgs_factor(length * length * mass / (time * time));
    torque.set_cgs_symbol("erg");
	torque.set_name("torque");

    force.set_cgs_factor(mass * length / (time * time));
    force.set_cgs_symbol("dyn");
	force.set_name("force");

    mass_accretion_rate.set_cgs_factor(mass / time);
    mass_accretion_rate.set_cgs_symbol("g s^-1");
	mass_accretion_rate.set_name("mass accretion rate");

    // after all units have calculated, calculate constants in code units
    constants::calculate_constants_in_code_units();
}

void print_code_units()
{
    logging::print_master(LOG_VERBOSE "Code units:\n");
    logging::print_master(LOG_VERBOSE
			  "                     length:       l0 = %.16g %s\n",
			  length.get_code_to_cgs_factor(), length.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                       mass:       m0 = %.16g %s\n",
			  mass.get_code_to_cgs_factor(), mass.get_cgs_symbol());
    logging::print_master(
	LOG_VERBOSE
	"                       time:       t0 = %.16g %s\t\t= %15.10g a\n",
	time.get_code_to_cgs_factor(), time.get_cgs_symbol(),
	time.get_code_to_cgs_factor() / (24 * 60 * 60 * 365.2425));
    logging::print_master(
	LOG_VERBOSE "                temperature:       T0 = %.16g %s\n",
	temperature.get_code_to_cgs_factor(), temperature.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                     energy:       E0 = %.16g %s\n",
			  energy.get_code_to_cgs_factor(), energy.get_cgs_symbol());
    logging::print_master(
	LOG_VERBOSE "             energy density:       e0 = %.16g %s\n",
	energy_density.get_code_to_cgs_factor(), energy_density.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                    density:     rho0 = %.16g %s\n",
			  density.get_code_to_cgs_factor(), density.get_cgs_symbol());
    logging::print_master(
	LOG_VERBOSE "            surface density:   Sigma0 = %.16g %s\n",
	surface_density.get_code_to_cgs_factor(), surface_density.get_cgs_symbol());
    logging::print_master(
	LOG_VERBOSE "                energy flux:       S0 = %.16g %s\n",
	energy_flux.get_code_to_cgs_factor(), energy_flux.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                   velocity:       v0 = %.16g %s\n",
			  velocity.get_code_to_cgs_factor(), velocity.get_cgs_symbol());
    logging::print_master(
	LOG_VERBOSE "               acceleration:       a0 = %.16g %s\n",
	acceleration.get_code_to_cgs_factor(), acceleration.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "        kinematic viscosity:      nu0 = %.16g %s\n",
			  kinematic_viscosity.get_code_to_cgs_factor(),
			  kinematic_viscosity.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                     stress:   sigma0 = %.16g %s\n",
			  stress.get_code_to_cgs_factor(), stress.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                   pressure:       p0 = %.16g %s\n",
			  pressure.get_code_to_cgs_factor(), pressure.get_cgs_symbol());
    logging::print_master(
	LOG_VERBOSE "           angular momentum:      L0 = %.16g %s\n",
	angular_momentum.get_code_to_cgs_factor(), angular_momentum.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                      power:       P0 = %.16g %s\n",
			  power.get_code_to_cgs_factor(), power.get_cgs_symbol());

    logging::print_master(
	LOG_VERBOSE "                      potential:       V0 = %.16g %s\n",
	potential.get_code_to_cgs_factor(), potential.get_cgs_symbol());

    logging::print_master(LOG_VERBOSE
			  "                     torque:     tau0 = %.16g %s\n",
			  torque.get_code_to_cgs_factor(), torque.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                      force:       F0 = %.16g %s\n",
			  force.get_code_to_cgs_factor(), force.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "        mass accretion rate:    Mdot0 = %.16g %s\n",
			  mass_accretion_rate.get_code_to_cgs_factor(),
			  mass_accretion_rate.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                    opacity:   kappa0 = %.16g %s\n",
			  opacity.get_code_to_cgs_factor(), opacity.get_cgs_symbol());
}



void write_code_units_file()
{
    /* Write a file containing the base constants to the output folder. */

    if (CPU_Master) {

	const std::string filename = output::outdir + "units.yml";

	std::ofstream of(filename);

	of.precision(std::numeric_limits<double>::max_digits10);
	of << "# code units file" << std::endl;
	of << "# version 0.2" << std::endl << std::endl;

	const std::string indent = "  ";

	std::vector<t_unit*> units = {
		&length,
		&mass,
		&time,
		&temperature,
		&energy,
		&energy_density,
		&density,
		&surface_density,
		&opacity,
		&energy_flux,
		&velocity,
		&angular_momentum,
		&kinematic_viscosity,
		&dynamic_viscosity,
		&acceleration,
		&stress,
		&pressure,
		&power,
		&potential,
		&torque,
		&force,
		&mass_accretion_rate,
	};
	
	for (auto unit : units) {
		of << unit->get_name() << ":" << std::endl;
		of << indent << "cgs symbol: " << unit->get_cgs_symbol() << std::endl;
		of << indent << "cgs value: " << unit->get_code_to_cgs_factor() << std::endl;
		of << std::endl;
	}

	of.close();
	}
    MPI_Barrier(MPI_COMM_WORLD);
}


} // namespace units
