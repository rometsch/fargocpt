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
	t_bc_scalar bc_sigma_inner;
	t_bc_scalar bc_sigma_outer;
	t_bc_scalar bc_energy_inner;
	t_bc_scalar bc_energy_outer;
	t_bc_vrad bc_vrad_inner;
	t_bc_vrad bc_vrad_outer;
	t_bc_vaz bc_vaz_inner;
	t_bc_vaz bc_vaz_outer;
	t_bc_special bc_special;

	static void sigma_inner() {
		const std::string str = config::cfg.get_lowercase("BoundarySigmaInner", "zerogradient");
		if (str == "zerogradient") {
			bc_sigma_inner = t_bc_scalar::bcs_zero_gradient;
		} else if (str == "diskmodel") {
			bc_sigma_inner = t_bc_scalar::bcs_disk_model;
		} else if (str == "reference") {
			bc_sigma_inner = t_bc_scalar::bcs_reference;
		} else {
			throw std::runtime_error("Unknown boundary condition for sigma inner: " + str);
		}
		logging::print_master(LOG_INFO "BC: Sigma inner = %s\n", str.c_str());
	}

	static void sigma_outer() {
		const std::string str = config::cfg.get_lowercase("BoundarySigmaOuter", "zerogradient");
		if (str == "zerogradient") {
			bc_sigma_outer = t_bc_scalar::bcs_zero_gradient;
		} else if (str == "diskmodel") {
			bc_sigma_outer = t_bc_scalar::bcs_disk_model;
		} else if (str == "reference") {
			bc_sigma_outer = t_bc_scalar::bcs_reference;
		} else {
			throw std::runtime_error("Unknown boundary condition for sigma outer: " + str);
		}
		logging::print_master(LOG_INFO "BC: Sigma outer = %s\n", str.c_str());
	}

	static void energy_inner() {
		const std::string str = config::cfg.get_lowercase("BoundaryEnergyInner", "zerogradient");
		if (str == "zerogradient") {
			bc_energy_inner = t_bc_scalar::bcs_zero_gradient;
		} else if (str == "diskmodel") {
			bc_energy_inner = t_bc_scalar::bcs_disk_model;
		} else if (str == "reference") {
			bc_energy_inner = t_bc_scalar::bcs_reference;
		} else {
			throw std::runtime_error("Unknown boundary condition for energy inner: " + str);
		}
		logging::print_master(LOG_INFO "BC: Energy inner = %s\n", str.c_str());
	}

	static void energy_outer() {
		const std::string str = config::cfg.get_lowercase("BoundaryEnergyOuter", "zerogradient");
		if (str == "zerogradient") {
			bc_energy_outer = t_bc_scalar::bcs_zero_gradient;
		} else if (str == "diskmodel") {
			bc_energy_outer = t_bc_scalar::bcs_disk_model;
		} else if (str == "reference") {
			bc_energy_outer = t_bc_scalar::bcs_reference;
		} else {
			throw std::runtime_error("Unknown boundary condition for energy outer: " + str);
		}
		logging::print_master(LOG_INFO "BC: Energy outer = %s\n", str.c_str());
	}

	static void vrad_inner() {
		const std::string str = config::cfg.get_lowercase("BoundaryVradInner", "zerogradient");
		if (str == "zerogradient") {
			bc_vrad_inner = t_bc_vrad::bcr_zero_gradient;
		} else if (str == "diskmodel") {
			bc_vrad_inner = t_bc_vrad::bcr_disk_model;
		} else if (str == "reference") {
			bc_vrad_inner = t_bc_vrad::bcr_reference;
		} else if (str == "reflective") {
			bc_vrad_inner = t_bc_vrad::bcr_reflective;
		} else if (str == "outflow") {
			bc_vrad_inner = t_bc_vrad::bcr_outflow;
		} else if (str == "viscous") {
			bc_vrad_inner = t_bc_vrad::bcr_viscous;
		} else if (str == "keplerian") {
			bc_vrad_inner = t_bc_vrad::bcr_Keplerian;
		} else {
			throw std::runtime_error("Unknown boundary condition for vrad inner: " + str);
		}
		logging::print_master(LOG_INFO "BC: Vrad inner = %s\n", str.c_str());
	}

	static void vrad_outer() {
		const std::string str = config::cfg.get_lowercase("BoundaryVradOuter", "zerogradient");
		if (str == "zerogradient") {
			bc_vrad_outer = t_bc_vrad::bcr_zero_gradient;
		} else if (str == "diskmodel") {
			bc_vrad_outer = t_bc_vrad::bcr_disk_model;
		} else if (str == "reference") {
			bc_vrad_outer = t_bc_vrad::bcr_reference;
		} else if (str == "reflective") {
			bc_vrad_outer = t_bc_vrad::bcr_reflective;
		} else if (str == "outflow") {
			bc_vrad_outer = t_bc_vrad::bcr_outflow;
		} else if (str == "viscous") {
			bc_vrad_outer = t_bc_vrad::bcr_viscous;
		} else if (str == "keplerian") {
			bc_vrad_outer = t_bc_vrad::bcr_Keplerian;
		} else {
			throw std::runtime_error("Unknown boundary condition for vrad outer: " + str);
		}
		logging::print_master(LOG_INFO "BC: Vrad outer = %s\n", str.c_str());
	}

	static void vaz_inner() {
		const std::string str = config::cfg.get_lowercase("BoundaryVazInner", "zerogradient");
		if (str == "zerogradient") {
			bc_vaz_inner = t_bc_vaz::bca_zero_gradient;
		} else if (str == "diskmodel") {
			bc_vaz_inner = t_bc_vaz::bca_disk_model;
		} else if (str == "reference") {
			bc_vaz_inner = t_bc_vaz::bca_reference;
		} else if (str == "zeroshear") {
			bc_vaz_inner = t_bc_vaz::bca_zero_shear;
		} else if (str == "balanced") {
			bc_vaz_inner = t_bc_vaz::bca_balanced;
		} else if (str == "keplerian") {
			bc_vaz_inner = t_bc_vaz::bca_Keplerian;
		} else {
			throw std::runtime_error("Unknown boundary condition for vaz inner: " + str);
		}
		logging::print_master(LOG_INFO "BC: Vaz inner = %s\n", str.c_str());
	}

	static void vaz_outer() {
		const std::string str = config::cfg.get_lowercase("BoundaryVazOuter", "zerogradient");
		if (str == "zerogradient") {
			bc_vaz_outer = t_bc_vaz::bca_zero_gradient;
		} else if (str == "diskmodel") {
			bc_vaz_outer = t_bc_vaz::bca_disk_model;
		} else if (str == "reference") {
			bc_vaz_outer = t_bc_vaz::bca_reference;
		} else if (str == "zeroshear") {
			bc_vaz_outer = t_bc_vaz::bca_zero_shear;
		} else if (str == "balanced") {
			bc_vaz_outer = t_bc_vaz::bca_balanced;
		} else if (str == "keplerian") {
			bc_vaz_outer = t_bc_vaz::bca_Keplerian;
		} else {
			throw std::runtime_error("Unknown boundary condition for vaz outer: " + str);
		}
		logging::print_master(LOG_INFO "BC: Vaz outer = %s\n", str.c_str());
	}

	static void special() {
		const std::string str = config::cfg.get_lowercase("BoundarySpecial", "none");
		if (str == "com") {
			bc_special = t_bc_special::bcsp_com;
		} else if (str == "massoverflow") {
			bc_special = t_bc_special::bcsp_mass_overflow;
		} else if (str == "none") {
			bc_special = t_bc_special::bcsp_none;
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
