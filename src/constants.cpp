/**
	\file constants.cpp
	\author Tobias Mueller <Tobias_Mueller@twam.info>

	This file handles all kinds physical constants used in this code
*/

#include "constants.h"
#include "LowTasks.h"
#include "logging.h"
#include "units.h"


#include <cstring>

namespace constants
{
#ifdef DO_USE_PLUTO_UNITS
#undef M_PI
#define M_PI 3.14159265358979
/// Source pluto.h
const double cgs_G = 6.6726e-8;

/// molecular mass (dalton) in cgs
const double cgs_m_u = 1.66053886e-24;

/// Boltzmann constant in cgs
const double cgs_k_B = 1.3806505e-16;
/// Planck constant in cgs
const double cgs_h = 6.62606876e-27;
/// speed of light in cgs
const double cgs_c = 2.99792458e10;
#else
const auto g = llnlunits::precise::g;
const auto s = llnlunits::precise::second;
const auto cm = llnlunits::precise::cm;
const auto K = llnlunits::precise::K;

/// gravitational constant in cgs
// const double cgs_G = 6.6738480e-8;
const double cgs_G = llnlunits::constants::G.value_as(cm*cm*cm/(g*s*s));
double global_cgs_G = cgs_G;

/// Boltzmann constant in cgs
// const double cgs_k_B = 1.380650424e-16;
const double cgs_k_B = llnlunits::constants::k.value_as(g*cm*cm/(K*s*s));

/// molecular mass in cgs
// const double cgs_m_u = 1.66053878283e-24;
// const double cgs_m_u = 1.6737236e-24;
const double cgs_m_u = llnlunits::constants::mu.value_as(g);

/// Planck constant in cgs
// const double cgs_h = 6.6260689633e-27;
const double cgs_h = llnlunits::constants::h.value_as(g*cm*cm/s);

/// speed of light in cgs
// const double cgs_c = 299792458.e2;
const double cgs_c = llnlunits::constants::c.value_as(cm/s);
#endif

t_constant k_B;
t_constant m_u;
t_constant h;
t_constant c;
t_constant sigma;
t_constant _G;
t_constant _R;

#ifndef NDEBUG
t_constant &G = _G;
t_constant &R = _R;
#endif

t_constant::t_constant()
{
    m_symbol = NULL;
    m_code_value = 1.0;
    m_cgs_value = 1.0;
    m_cgs_unit_symbol = NULL;
}

t_constant::~t_constant()
{
    delete[] m_symbol;
    delete[] m_cgs_unit_symbol;
}

void t_constant::set_symbol(const char *symbol)
{
    // delete old symbol
    delete[] m_symbol;

    // aquire memory for symbol
    const unsigned int length = strlen(symbol) + 1;
    m_symbol = new char[length];

    // copy symbol
    strncpy(m_symbol, symbol, length);
}

void t_constant::set_code_value(double value) { m_code_value = value; }

void t_constant::set_cgs_value(double value) { m_cgs_value = value; }

void t_constant::set_cgs_unit_symbol(const char *symbol)
{
    // delete old symbol
    delete[] m_cgs_unit_symbol;

    // aquire memory for symbol
    const unsigned int length = strlen(symbol) + 1;
    m_cgs_unit_symbol = new char[length];

    // copy symbol
    strncpy(m_cgs_unit_symbol, symbol, length);
}

const char *t_constant::get_symbol(void) const { return m_symbol; }

double t_constant::get_code_value() const { return m_code_value; }

double t_constant::get_cgs_value() const { return m_cgs_value; }

const char *t_constant::get_cgs_unit_symbol(void) const
{
    return m_cgs_unit_symbol;
}

/**
	Initialize constants objects with symbol, cgs value and cgs unit symbol
*/
void initialize_constants()
{
    _G.set_symbol("G");
    _G.set_cgs_value(cgs_G);
    _G.set_cgs_unit_symbol("cm^3 g^-1 s^-2");

    k_B.set_symbol("k_B");
    k_B.set_cgs_value(cgs_k_B);
    k_B.set_cgs_unit_symbol("erg K^-1");

    m_u.set_symbol("m_u");
    m_u.set_cgs_value(cgs_m_u);
    m_u.set_cgs_unit_symbol("g");

    h.set_symbol("h");
    h.set_cgs_value(cgs_h);
    h.set_cgs_unit_symbol("erg s");

    c.set_symbol("c");
    c.set_cgs_value(cgs_c);
    c.set_cgs_unit_symbol("cm s^-1");

    _R.set_symbol("R");
    _R.set_cgs_value(k_B.get_cgs_value() / (m_u.get_cgs_value()));
    _R.set_cgs_unit_symbol("erg K^-1 g^-1");

    sigma.set_symbol("sigma");
#ifdef DO_USE_PLUTO_UNITS
	sigma.set_cgs_value(5.67051e-5);
#else
    sigma.set_cgs_value(
	2. * pow(M_PI, 5) * pow(k_B.get_cgs_value(), 4) /
	(15. * pow(h.get_cgs_value(), 3) * pow(c.get_cgs_value(), 2)));
#endif
    sigma.set_cgs_unit_symbol("erg cm^-2 s^-1 K^-4");
}

/**
	Calculate constants in code units from cgs units. Must be called AFTER
   all units have been initalized properly.
*/
void calculate_constants_in_code_units()
{
#ifndef NDEBUG
    G.set_code_value(G.get_cgs_value() /
		     (units::length * units::length * units::length /
		      (units::mass * units::time * units::time)));
#endif
    k_B.set_code_value(k_B.get_cgs_value() /
		       (units::energy / units::temperature));
    m_u.set_code_value(m_u.get_cgs_value() / (units::mass));
    h.set_code_value(h.get_cgs_value() / (units::energy * units::time));
    c.set_code_value(c.get_cgs_value() / (units::length / units::time));
#ifndef NDEBUG
    R.set_code_value(R.get_cgs_value() /
		     (units::energy / (units::temperature * units::mass)));
#endif
    sigma.set_code_value(
	sigma.get_cgs_value() /
	(units::energy /
	 (units::length * units::length * units::time * units::temperature *
	  units::temperature * units::temperature * units::temperature)));
}

void print_constants()
{
#ifdef DO_USE_PLUTO_UNITS
logging::print_master(LOG_INFO "Using PLUTO like units\n");
#endif
    logging::print_master(LOG_VERBOSE "Code constants:\n");
#ifndef NDEBUG
    logging::print_master(LOG_VERBOSE
			  "     gravitational constant: %8s = %15g = %15g %s\n",
			  G.get_symbol(), G.get_code_value(), G.get_cgs_value(),
			  G.get_cgs_unit_symbol());
#else
    logging::print_master(
	LOG_VERBOSE
	"     gravitational constant: 1 (hardcoded, compile without NDEBUG to calculate it dynamically)\n");
#endif
    logging::print_master(LOG_VERBOSE
			  "         Boltzmann constant: %8s = %15g = %15g %s\n",
			  k_B.get_symbol(), k_B.get_code_value(),
			  k_B.get_cgs_value(), k_B.get_cgs_unit_symbol());
    logging::print_master(LOG_VERBOSE
			  "             molecular mass: %8s = %15g = %15g %s\n",
			  m_u.get_symbol(), m_u.get_code_value(),
			  m_u.get_cgs_value(), m_u.get_cgs_unit_symbol());
    logging::print_master(LOG_VERBOSE
			  "            Planck constant: %8s = %15g = %15g %s\n",
			  h.get_symbol(), h.get_code_value(), h.get_cgs_value(),
			  h.get_cgs_unit_symbol());
    logging::print_master(LOG_VERBOSE
			  "             speed of light: %8s = %15g = %15g %s\n",
			  c.get_symbol(), c.get_code_value(), c.get_cgs_value(),
			  c.get_cgs_unit_symbol());
#ifndef NDEBUG
    logging::print_master(LOG_VERBOSE
			  "      specific gas constant: %8s = %15g = %15g %s\n",
			  R.get_symbol(), R.get_code_value(), R.get_cgs_value(),
			  R.get_cgs_unit_symbol());
#else
    logging::print_master(
	LOG_VERBOSE
	"      specific gas constant: 1 (hardcoded, compile without NDEBUG to calculate it dynamically)\n");
#endif
    logging::print_master(LOG_VERBOSE
			  "  Stefan-Boltzmann constant: %8s = %15g = %15g %s\n",
			  sigma.get_symbol(), sigma.get_code_value(),
			  sigma.get_cgs_value(), sigma.get_cgs_unit_symbol());
}

} // namespace constants
