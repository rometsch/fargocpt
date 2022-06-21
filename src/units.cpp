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
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string.h>

namespace units
{

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

t_unit::t_unit()
{
    m_cgs_symbol = NULL;
    m_cgs_factor = 1;
}

t_unit::~t_unit() { delete[] m_cgs_symbol; }

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
	get unit symbol in cgs system
*/
const char *t_unit::get_cgs_symbol(void) const { return m_cgs_symbol; }

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

#include <stdio.h>

void calculate_unit_factors()
{

    length.set_cgs_factor(parameters::L0 * cgs_AU);
    length.set_cgs_symbol("cm");

    mass.set_cgs_factor(parameters::M0 * cgs_Msol);
    mass.set_cgs_symbol("g");

    time.set_cgs_factor(sqrt((length * length * length) /
			     (mass * constants::_G.get_cgs_value())));
    time.set_cgs_symbol("s");

    energy.set_cgs_factor(length * length * mass / (time * time));
    energy.set_cgs_symbol("erg");

    energy_density.set_cgs_factor(mass / (time * time));
    energy_density.set_cgs_symbol("erg cm^-2");

    temperature.set_cgs_factor((constants::_G.get_cgs_value() * mass) /
			       (constants::_R.get_cgs_value() * length));
    temperature.set_cgs_symbol("K");

    density.set_cgs_factor(mass / (length * length * length));
    density.set_cgs_symbol("g cm^-3");

    surface_density.set_cgs_factor(mass / (length * length));
    surface_density.set_cgs_symbol("g cm^-2");

    opacity.set_cgs_factor(length * length / mass);
    opacity.set_cgs_symbol("g^-1 cm^2");

    energy_flux.set_cgs_factor(energy / (length * length * time));
    energy_flux.set_cgs_symbol("erg cm^-2 s^-1");

    velocity.set_cgs_factor(length / time);
    velocity.set_cgs_symbol("cm s^-1");

    acceleration.set_cgs_factor(length / (time * time));
    acceleration.set_cgs_symbol("cm s^-2");

    angular_momentum.set_cgs_factor(length * mass * velocity);
    angular_momentum.set_cgs_symbol("cm^2 g s^-1");

    kinematic_viscosity.set_cgs_factor(length * length / time);
    kinematic_viscosity.set_cgs_symbol("cm^2 s^-1");

    dynamic_viscosity.set_cgs_factor(mass / (length * time));
    dynamic_viscosity.set_cgs_symbol("P");

    stress.set_cgs_factor(mass / (time * time));
    stress.set_cgs_symbol("g s^-2");

    pressure.set_cgs_factor(mass / (time * time));
    pressure.set_cgs_symbol("dyn cm^-1");

    power.set_cgs_factor(mass * length * length / (time * time * time));
    power.set_cgs_symbol("erg/s");

    potential.set_cgs_factor(length * length / (time * time));
    potential.set_cgs_symbol("erg/g");

    torque.set_cgs_factor(length * length * mass / (time * time));
    torque.set_cgs_symbol("erg");

    force.set_cgs_factor(mass * length / (time * time));
    force.set_cgs_symbol("dyn");

    mass_accretion_rate.set_cgs_factor(mass / time);
    mass_accretion_rate.set_cgs_symbol("g s^-1");

    // after all units have calculated, calculate constants in code units
    constants::calculate_constants_in_code_units();
}

void print_code_units()
{
    logging::print_master(LOG_VERBOSE "Code units:\n");
    logging::print_master(LOG_VERBOSE
			  "                     length:       l0 = %.16g %s\n",
			  length.get_cgs_factor(), length.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                       mass:       m0 = %.16g %s\n",
			  mass.get_cgs_factor(), mass.get_cgs_symbol());
    logging::print_master(
	LOG_VERBOSE
	"                       time:       t0 = %.16g %s\t\t= %15.10g a\n",
	time.get_cgs_factor(), time.get_cgs_symbol(),
	time.get_cgs_factor() / (24 * 60 * 60 * 365.2425));
    logging::print_master(
	LOG_VERBOSE "                temperature:       T0 = %.16g %s\n",
	temperature.get_cgs_factor(), temperature.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                     energy:       E0 = %.16g %s\n",
			  energy.get_cgs_factor(), energy.get_cgs_symbol());
    logging::print_master(
	LOG_VERBOSE "             energy density:       e0 = %.16g %s\n",
	energy_density.get_cgs_factor(), energy_density.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                    density:     rho0 = %.16g %s\n",
			  density.get_cgs_factor(), density.get_cgs_symbol());
    logging::print_master(
	LOG_VERBOSE "            surface density:   Sigma0 = %.16g %s\n",
	surface_density.get_cgs_factor(), surface_density.get_cgs_symbol());
    logging::print_master(
	LOG_VERBOSE "                energy flux:       S0 = %.16g %s\n",
	energy_flux.get_cgs_factor(), energy_flux.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                   velocity:       v0 = %.16g %s\n",
			  velocity.get_cgs_factor(), velocity.get_cgs_symbol());
    logging::print_master(
	LOG_VERBOSE "               acceleration:       a0 = %.16g %s\n",
	acceleration.get_cgs_factor(), acceleration.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "        kinematic viscosity:      nu0 = %.16g %s\n",
			  kinematic_viscosity.get_cgs_factor(),
			  kinematic_viscosity.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                     stress:   sigma0 = %.16g %s\n",
			  stress.get_cgs_factor(), stress.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                   pressure:       p0 = %.16g %s\n",
			  pressure.get_cgs_factor(), pressure.get_cgs_symbol());
    logging::print_master(
	LOG_VERBOSE "           angular momentum:      L0 = %.16g %s\n",
	angular_momentum.get_cgs_factor(), angular_momentum.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                      power:       P0 = %.16g %s\n",
			  power.get_cgs_factor(), power.get_cgs_symbol());

    logging::print_master(
	LOG_VERBOSE "                      potential:       V0 = %.16g %s\n",
	potential.get_cgs_factor(), potential.get_cgs_symbol());

    logging::print_master(LOG_VERBOSE
			  "                     torque:     tau0 = %.16g %s\n",
			  torque.get_cgs_factor(), torque.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                      force:       F0 = %.16g %s\n",
			  force.get_cgs_factor(), force.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "        mass accretion rate:    Mdot0 = %.16g %s\n",
			  mass_accretion_rate.get_cgs_factor(),
			  mass_accretion_rate.get_cgs_symbol());
    logging::print_master(LOG_VERBOSE
			  "                    opacity:   kappa0 = %.16g %s\n",
			  opacity.get_cgs_factor(), opacity.get_cgs_symbol());
}

void write_code_unit_file()
{
    /* Write a file containing the base units to the output folder. */

    FILE *fd = 0;
    char *fd_filename;

    if (CPU_Master) {
	if (asprintf(&fd_filename, "%s%s", OUTPUTDIR, "units.dat") == -1) {
	    logging::print_master(LOG_ERROR
				  "Not enough memory for string buffer.\n");
	    PersonalExit(1);
	}
	fd = fopen(fd_filename, "w");
	if (fd == NULL) {
	    logging::print_master(LOG_ERROR
				  "Can't write 'units.dat' file. Aborting.\n");
	    PersonalExit(1);
	}

	free(fd_filename);

	fprintf(fd, "# units-file : 1.0\n");
	fprintf(fd, "# log output of units:\n");
	fprintf(fd, "#                    length:       l0 = %.16g %s\n",
		length.get_cgs_factor(), length.get_cgs_symbol());
	fprintf(fd, "#                      mass:       m0 = %.16g %s\n",
		mass.get_cgs_factor(), mass.get_cgs_symbol());
	fprintf(
	    fd,
	    "#                      time:       t0 = %.16g %s\t\t= %15.10g a\n",
	    time.get_cgs_factor(), time.get_cgs_symbol(),
	    time.get_cgs_factor() / (24 * 60 * 60 * 365.2425));
	fprintf(fd, "#               temperature:       T0 = %.16g %s\n",
		temperature.get_cgs_factor(), temperature.get_cgs_symbol());
	fprintf(fd, "#                    energy:       E0 = %.16g %s\n",
		energy.get_cgs_factor(), energy.get_cgs_symbol());
	fprintf(fd, "#            energy density:       e0 = %.16g %s\n",
		energy_density.get_cgs_factor(),
		energy_density.get_cgs_symbol());
	fprintf(fd, "#                   density:     rho0 = %.16g %s\n",
		density.get_cgs_factor(), density.get_cgs_symbol());
	fprintf(fd, "#           surface density:   Sigma0 = %.16g %s\n",
		surface_density.get_cgs_factor(),
		surface_density.get_cgs_symbol());
	fprintf(fd, "#               energy flux:       S0 = %.16g %s\n",
		energy_flux.get_cgs_factor(), energy_flux.get_cgs_symbol());
	fprintf(fd, "#                  velocity:       v0 = %.16g %s\n",
		velocity.get_cgs_factor(), velocity.get_cgs_symbol());
	fprintf(fd, "#              acceleration:       a0 = %.16g %s\n",
		acceleration.get_cgs_factor(), acceleration.get_cgs_symbol());
	fprintf(fd, "#       kinematic viscosity:      nu0 = %.16g %s\n",
		kinematic_viscosity.get_cgs_factor(),
		kinematic_viscosity.get_cgs_symbol());
	fprintf(fd, "#                    stress:   sigma0 = %.16g %s\n",
		stress.get_cgs_factor(), stress.get_cgs_symbol());
	fprintf(fd, "#                  pressure:       p0 = %.16g %s\n",
		pressure.get_cgs_factor(), pressure.get_cgs_symbol());

	fprintf(fd, "#          angular momentum:       L0 = %.16g %s\n",
		angular_momentum.get_cgs_factor(),
		angular_momentum.get_cgs_symbol());
	fprintf(fd, "#                     power:       P0 = %.16g %s\n",
		power.get_cgs_factor(), power.get_cgs_symbol());
	fprintf(fd, "#                     potential:	V0 = %.16g %s\n",
		potential.get_cgs_factor(), potential.get_cgs_symbol());
	fprintf(fd, "#                    torque:     tau0 = %.16g %s\n",
		torque.get_cgs_factor(), torque.get_cgs_symbol());
	fprintf(fd, "#                     force:       F0 = %.16g %s\n",
		force.get_cgs_factor(), force.get_cgs_symbol());
	fprintf(fd, "#       mass_accretion_rate:    Mdot0 = %.16g %s\n",
		mass_accretion_rate.get_cgs_factor(),
		mass_accretion_rate.get_cgs_symbol());
	fprintf(fd, "#                   opacity:   kappa0 = %.16g %s\n",
		opacity.get_cgs_factor(), opacity.get_cgs_symbol());

	fprintf(fd, "# Syntax: base unit <tab> value <tab> unit name\n");
	fprintf(fd, "length\t%.16e\t%s\n", length.get_cgs_factor(),
		length.get_cgs_symbol());
	fprintf(fd, "mass\t%.16e\t%s\n", mass.get_cgs_factor(),
		mass.get_cgs_symbol());
	fprintf(fd, "time\t%.16e\t%s\n", time.get_cgs_factor(),
		time.get_cgs_symbol());
	fprintf(fd, "current\t \t \n");
	fprintf(fd, "temperature\t%.16e\t%s\n", temperature.get_cgs_factor(),
		temperature.get_cgs_symbol());
	fprintf(fd, "amount\t \t \n");
	fprintf(fd, "intensity\t \t ");
	fclose(fd);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace units
