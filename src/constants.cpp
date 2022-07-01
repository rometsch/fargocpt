/**
	\file constants.cpp
	\author Tobias Mueller <Tobias_Mueller@twam.info>

	This file handles all kinds physical constants used in this code
*/

#include "constants.h"
#include "LowTasks.h"
#include "logging.h"
#include "units.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

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
/// gravitational constant in cgs source: https://ssd.jpl.nasa.gov/astro_par.html
const double cgs_G = 6.67430e-8;

/// molecular mass (dalton) in cgs
/// Source: https://physics.nist.gov/cgi-bin/cuu/Value?u
const double cgs_m_u = 1.66053906660e-24;

/// Source: CODATA recommended values of the fundamental physical constants: 2018
/// https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.93.025010#fulltext
/// Boltzmann constant in cgs
const double cgs_k_B = 1.380649e-16;
/// Planck constant in cgs
const double cgs_h = 6.62607015e-27;
/// speed of light in cgs
const double cgs_c = 299792458.e2;
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
    m_symbol = new char[strlen(symbol) + 1];

    // copy symbol
    strcpy(m_symbol, symbol);
}

void t_constant::set_code_value(double value) { m_code_value = value; }

void t_constant::set_cgs_value(double value) { m_cgs_value = value; }

void t_constant::set_cgs_unit_symbol(const char *symbol)
{
    // delete old symbol
    delete[] m_cgs_unit_symbol;

    // aquire memory for symbol
    m_cgs_unit_symbol = new char[strlen(symbol) + 1];

    // copy symbol
    strcpy(m_cgs_unit_symbol, symbol);
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

	/// Stefan Boltzmann constant
	sigma.set_symbol("sigma");
	sigma.set_cgs_unit_symbol("erg cm^-2 s^-1 K^-4");
#ifdef DO_USE_PLUTO_UNITS
	sigma.set_cgs_value(5.67051e-5);
#else
	sigma.set_cgs_value(
	2. * pow(M_PI, 5) * pow(k_B.get_cgs_value(), 4) /
	(15. * pow(h.get_cgs_value(), 3) * pow(c.get_cgs_value(), 2)));
#endif
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
