#pragma once

#include <string>
#include <type_traits>

#define UNITS_DEFAULT_DOMAIN llnlunits::domains::astronomy
#include "units/units.hpp"


namespace units
{

typedef llnlunits::precise_unit precise_unit;
namespace precise = llnlunits::precise;

#ifdef DO_USE_PLUTO_UNITS
/// Source: pluto.h
const double cgs_AU = 1.49597892e13;
const double cgs_Msol = 2.e33;
const double cgs_Year = 365.25636 * 24.0 * 60.0 * 60.0; // pluto.h does not have year defined
#else
/// astronomical unit in cgs
/// source: https://ssd.jpl.nasa.gov/astro_par.html
// const double cgs_AU = 149.60e11;
const double cgs_AU = 1.495978707e13;

/// solar mass in cgs
const double cgs_Msol = 1.98892e33;

/// seconds of a year
const double cgs_Year = 31556925.261; // 365.*24.*60.*60.;
#endif

const double solar_radius_in_au = 0.00465047;

class t_unit
{
  private:
    /// cgs conversion factor
    double m_cgs_factor;
    double m_inverse_cgs_factor;

    /// cgs unit symbol
    char *m_cgs_symbol;

  public:
    t_unit();
    ~t_unit();

    // setter
    void set_cgs_factor(double);
    void set_cgs_symbol(const char *);

    // getter
    const char *get_cgs_symbol(void) const;
    /// get conversion factor to cgs system
    inline double get_cgs_factor(void) const { return m_cgs_factor; }
    /// get conversion factor from cgs systems
    inline double get_inverse_cgs_factor(void) const
    {
	return m_inverse_cgs_factor;
    }

    // operator
    inline operator const double &() const { return m_cgs_factor; }
    /* inline operator double&() { return m_cgs_factor; } */

    std::string get_cgs_factor_symbol();
};


extern llnlunits::precise_unit L0;
extern llnlunits::precise_unit M0;
extern llnlunits::precise_unit T0;
extern llnlunits::precise_unit Temp0;

extern llnlunits::precise_unit solMass;
extern llnlunits::precise_unit au;

extern t_unit not_yet_defined;
extern t_unit unitless;
extern t_unit length;
extern t_unit mass;
extern t_unit time;
extern t_unit temperature;
extern t_unit energy;
extern t_unit energy_density;
extern t_unit density;
extern t_unit surface_density;
extern t_unit opacity;
extern t_unit potential;
extern t_unit energy_flux;
extern t_unit velocity;
extern t_unit angular_momentum;
extern t_unit kinematic_viscosity;
extern t_unit dynamic_viscosity;
extern t_unit acceleration;
extern t_unit stress;
extern t_unit pressure;
extern t_unit power;
extern t_unit torque;
extern t_unit force;
extern t_unit mass_accretion_rate;

bool has_unit(const std::string &val);

void set_baseunits(
  const std::string &l0, 
  const std::string &m0,
  const std::string &t0,
  const std::string &temp0);

template <typename T> T parse_units(const std::string &val);
template <typename T> T parse_units(const std::string &val, const precise_unit &baseunit);

void calculate_unit_factors();
void print_code_units();
void write_code_unit_file();

} // namespace units
