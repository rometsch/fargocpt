/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"

#include "../global.h"
#include "../logging.h"
#include "../config.h"
#include "../logging.h"

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


	std::string sigma_inner_name;
	std::string sigma_outer_name;
	std::string energy_inner_name;
	std::string energy_outer_name;
	std::string vrad_inner_name;
	std::string vrad_outer_name;
	std::string vaz_inner_name;
	std::string vaz_outer_name;



	static void sigma_inner() {
		const std::string str = config::cfg.get_lowercase("BoundarySigmaInner", "zerogradient");
		sigma_inner_name = str;
		if (str == "zerogradient") {
			sigma_inner_func = zero_gradient_inner;
		} else if (str == "diskmodel") {
		} else if (str == "reference") {
		} else {
			throw std::runtime_error("Unknown boundary condition for sigma inner: " + str);
		}
		logging::print_master(LOG_INFO "BC: Sigma inner = %s\n", str.c_str());
	}

	static void sigma_outer() {
		const std::string str = config::cfg.get_lowercase("BoundarySigmaOuter", "zerogradient");
		sigma_outer_name = str;
		if (str == "zerogradient") {
			sigma_outer_func = zero_gradient_outer;
		} else if (str == "diskmodel") {
		} else if (str == "reference") {
		} else {
			throw std::runtime_error("Unknown boundary condition for sigma outer: " + str);
		}
		logging::print_master(LOG_INFO "BC: Sigma outer = %s\n", str.c_str());
	}

	static void energy_inner() {
		const std::string str = config::cfg.get_lowercase("BoundaryEnergyInner", "zerogradient");
		energy_inner_name = str;
		if (str == "zerogradient") {
			energy_inner_func = zero_gradient_inner;
		} else if (str == "diskmodel") {
		} else if (str == "reference") {
		} else {
			throw std::runtime_error("Unknown boundary condition for energy inner: " + str);
		}
		logging::print_master(LOG_INFO "BC: Energy inner = %s\n", str.c_str());
	}

	static void energy_outer() {
		const std::string str = config::cfg.get_lowercase("BoundaryEnergyOuter", "zerogradient");
		energy_outer_name = str;
		if (str == "zerogradient") {
			energy_outer_func = zero_gradient_outer;
		} else if (str == "diskmodel") {
		} else if (str == "reference") {
		} else {
			throw std::runtime_error("Unknown boundary condition for energy outer: " + str);
		}
		logging::print_master(LOG_INFO "BC: Energy outer = %s\n", str.c_str());
	}

	static void vrad_inner() {
		const std::string str = config::cfg.get_lowercase("BoundaryVradInner", "zerogradient");
		if (str == "zerogradient") {
			vrad_inner_func = zero_gradient_inner;
		} else if (str == "diskmodel") {
		} else if (str == "reference") {
		} else if (str == "reflective") {
		} else if (str == "outflow") {
			vrad_inner_func = outflow_inner;
		} else if (str == "viscous") {
		} else if (str == "keplerian") {
		} else {
			throw std::runtime_error("Unknown boundary condition for vrad inner: " + str);
		}
		logging::print_master(LOG_INFO "BC: Vrad inner = %s\n", str.c_str());
	}

	static void vrad_outer() {
		const std::string str = config::cfg.get_lowercase("BoundaryVradOuter", "zerogradient");
		if (str == "zerogradient") {
			vrad_outer_func = zero_gradient_outer;
		} else if (str == "diskmodel") {
		} else if (str == "reference") {
		} else if (str == "reflective") {
		} else if (str == "outflow") {
			vrad_outer_func = outflow_outer;
		} else if (str == "viscous") {
		} else if (str == "keplerian") {
		} else {
			throw std::runtime_error("Unknown boundary condition for vrad outer: " + str);
		}
		logging::print_master(LOG_INFO "BC: Vrad outer = %s\n", str.c_str());
	}

	static void vaz_inner() {
		const std::string str = config::cfg.get_lowercase("BoundaryVazInner", "zerogradient");
		if (str == "zerogradient") {
			vaz_inner_func = zero_gradient_inner;
		} else if (str == "diskmodel") {
		} else if (str == "reference") {
		} else if (str == "zeroshear") {
		} else if (str == "balanced") {
		} else if (str == "keplerian") {
		} else {
			throw std::runtime_error("Unknown boundary condition for vaz inner: " + str);
		}
		logging::print_master(LOG_INFO "BC: Vaz inner = %s\n", str.c_str());
	}

	static void vaz_outer() {
		const std::string str = config::cfg.get_lowercase("BoundaryVazOuter", "zerogradient");
		if (str == "zerogradient") {
			vaz_outer_func = zero_gradient_outer;
		} else if (str == "diskmodel") {
		} else if (str == "reference") {
		} else if (str == "zeroshear") {
		} else if (str == "balanced") {
		} else if (str == "keplerian") {
		} else {
			throw std::runtime_error("Unknown boundary condition for vaz outer: " + str);
		}
		logging::print_master(LOG_INFO "BC: Vaz outer = %s\n", str.c_str());
	}

	static void special() {
		const std::string str = config::cfg.get_lowercase("BoundarySpecial", "none");
		if (str == "com") {
		} else if (str == "massoverflow") {
		} else if (str == "none") {
		} else {
			throw std::runtime_error("Unknown boundary condition for special: " + str);
		}
		logging::print_master(LOG_INFO "BC: Special = %s\n", str.c_str());
	}


	void parse_config() {

		sigma_inner();
		sigma_outer();

		energy_inner();
		energy_outer();

		vrad_inner();
		vrad_outer();

		vaz_inner();
		vaz_outer();

		special();
	}

} // namespace boundary_conditions
