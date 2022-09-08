/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"
#include "Theo.h"
#include "find_cell_id.h"
#include "global.h"
#include "logging.h"
#include "parameters.h"
#include "util.h"
#include "quantities.h"
#include <algorithm>
#include <cstring>
#include <cmath>
#include <vector>

#include "constants.h"
#include <cassert>
#include <experimental/filesystem>
#include <fstream>
#include <iostream>

// temporary
#include "LowTasks.h"
#include "SideEuler.h"
extern boolean OuterSourceMass;
extern std::vector<parameters::t_DampingType> parameters::damping_vector;

namespace boundary_conditions
{

void init_prescribed_time_variable_boundaries(t_data &data)
{

    // delete old data
    if (data[t_data::PRESCRIBED_DENSITY_OUTER].Field != nullptr) {
	delete[] data[t_data::PRESCRIBED_DENSITY_OUTER].Field;
	data[t_data::PRESCRIBED_DENSITY_OUTER].Field = nullptr;
    }
    if (data[t_data::PRESCRIBED_ENERGY_OUTER].Field != nullptr) {
	delete[] data[t_data::PRESCRIBED_ENERGY_OUTER].Field;
	data[t_data::PRESCRIBED_ENERGY_OUTER].Field = nullptr;
    }
    if (data[t_data::PRESCRIBED_V_RADIAL_OUTER].Field != nullptr) {
	delete[] data[t_data::PRESCRIBED_V_RADIAL_OUTER].Field;
	data[t_data::PRESCRIBED_V_RADIAL_OUTER].Field = nullptr;
    }
    if (data[t_data::PRESCRIBED_V_AZIMUTHAL_OUTER].Field != nullptr) {
	delete[] data[t_data::PRESCRIBED_V_AZIMUTHAL_OUTER].Field;
	data[t_data::PRESCRIBED_V_AZIMUTHAL_OUTER].Field = nullptr;
    }

    if (CPU_Rank != 0 && CPU_Rank != CPU_Highest)
	return;

    if (CPU_Rank == 0) {
	if (parameters::boundary_inner ==
	    parameters::boundary_condition_precribed_time_variable) {
	    die("Inner precribed time variable boundary condition is not implemented yet!\n");
	}
    }

    if (CPU_Rank == CPU_Highest) {
	if (parameters::boundary_outer ==
	    parameters::boundary_condition_precribed_time_variable) {
	    if (PRESCRIBED_BOUNDARY_OUTER_FILE == NULL) {
		die("Outer prescribed time variable boundary condition is enabled but the supplied file folder is not found!\n");
	    } else {

		// TODO: naming convention might need adjustment
		const int Nphi =
		    data[t_data::PRESCRIBED_DENSITY_OUTER].get_size_azimuthal();
		char *file_name_body_char;
		asprintf(&file_name_body_char, "%s/%dshift",
			 PRESCRIBED_BOUNDARY_OUTER_FILE, Nphi);
		std::string file_name_body{file_name_body_char};

		char *file_name_test_char;
		asprintf(&file_name_test_char, "%s0.dat", file_name_body_char);
		if (!std::experimental::filesystem::exists(
			file_name_test_char)) {
		    die("Prescribed boundary file %s does not exist!\n",
			file_name_test_char);
		}

		// get number of files
		int num_files = 0;
		const std::experimental::filesystem::path File_Folder{
		    PRESCRIBED_BOUNDARY_OUTER_FILE};
		for (auto const &dir_entry :
		     std::experimental::filesystem::directory_iterator{
			 File_Folder}) {
		    std::string path_string{dir_entry.path()};
		    if (path_string.find(file_name_body) != std::string::npos) {
			num_files++;
		    }
		}

		PRESCRIBED_TIME_SEGMENT_NUMBER = num_files;

		// load data from files
		int num_cells = num_files * Nphi;

		data[t_data::PRESCRIBED_DENSITY_OUTER].Nrad = num_files;
		data[t_data::PRESCRIBED_ENERGY_OUTER].Nrad = num_files;
		data[t_data::PRESCRIBED_V_RADIAL_OUTER].Nrad = num_files;
		data[t_data::PRESCRIBED_V_AZIMUTHAL_OUTER].Nrad = num_files;

		// assign new memory
		data[t_data::PRESCRIBED_DENSITY_OUTER].Field =
		    new double[num_cells];
		data[t_data::PRESCRIBED_ENERGY_OUTER].Field =
		    new double[num_cells];
		data[t_data::PRESCRIBED_V_RADIAL_OUTER].Field =
		    new double[num_cells];
		data[t_data::PRESCRIBED_V_AZIMUTHAL_OUTER].Field =
		    new double[num_cells];

		// set to 0
		data[t_data::PRESCRIBED_DENSITY_OUTER].clear();
		data[t_data::PRESCRIBED_ENERGY_OUTER].clear();
		data[t_data::PRESCRIBED_V_RADIAL_OUTER].clear();
		data[t_data::PRESCRIBED_V_AZIMUTHAL_OUTER].clear();

		std::vector<int> count_files;

		// read data
		for (auto const &dir_entry :
		     std::experimental::filesystem::directory_iterator{
			 File_Folder}) {
		    std::string path_string{dir_entry.path()};
		    if (path_string.find(file_name_body) != std::string::npos) {
			int file_id;
			std::string filename = file_name_body + "%d.dat";
			std::sscanf(path_string.c_str(), filename.c_str(),
				    &file_id);

			count_files.push_back(file_id);

			if (file_id >= num_files) {
			    die("file_id %d is out of range for %d number of files!\n",
				file_id, num_files);
			}

			std::fstream infile(path_string, std::ios_base::in);
			std::string line;
			int n_azimuthal = 0;
			while (std::getline(infile, line)) {

			    // TODO for now, hardcoded units for the loaded data
			    const double pluto_m0 = 1.8 * units::cgs_Msol;
			    const double pluto_l0 = 20.0 * units::cgs_AU;
			    const double pluto_sigma =
				pluto_m0 / pluto_l0 / pluto_l0;
			    const double pluto_t0 = std::sqrt(
				pluto_l0 * pluto_l0 * pluto_l0 /
				(pluto_m0 * constants::_G.get_cgs_value()));
			    const double pluto_v0 = pluto_l0 / pluto_t0;
			    const double pluto_aspect_ratio = 0.1;

			    double phi;
			    double sigma;
			    double vr;
			    double vphi;

			    sscanf(line.c_str(), "%lf	%lf	%lf	%lf\n",
				   &phi, &sigma, &vr, &vphi);

			    /*
			    double phi_calc = 2*M_PI / (double)(Nphi-0.5) *
			    (double(n_azimuthal) - 1.0); if(file_id == 0)
			    printf("phi = (%.5e	%.5e)	%.5e	%.5e	%.5e
			    %.5e\n", phi, phi_calc, (phi-phi_calc)/phi, sigma,
			    vr,vphi);
				       */

			    sigma = std::max(sigma * pluto_sigma *
						 units::surface_density
						     .get_inverse_cgs_factor(),
					     parameters::sigma_floor *
						 parameters::sigma0);
			    vr = vr * pluto_v0 *
				 units::velocity.get_inverse_cgs_factor();
			    vphi = vphi * pluto_v0 *
				   units::velocity.get_inverse_cgs_factor();

				const double m_bin = data.get_planetary_system().get_mass(2);
			    const double Cs =
				pluto_aspect_ratio *
				sqrt(constants::G * m_bin / RMAX) *
				std::pow(RMAX, FLARINGINDEX);
			    const double P = sigma * Cs * Cs;
			    const double T =
				parameters::MU / constants::R * P / sigma;
			    const double energy = T * sigma / parameters::MU *
						  constants::R /
						  (ADIABATICINDEX - 1.0);

			    data[t_data::PRESCRIBED_DENSITY_OUTER](
				file_id, n_azimuthal) = sigma;
			    data[t_data::PRESCRIBED_ENERGY_OUTER](
				file_id, n_azimuthal) = energy;
			    data[t_data::PRESCRIBED_V_RADIAL_OUTER](
				file_id, n_azimuthal) = vr;
			    data[t_data::PRESCRIBED_V_AZIMUTHAL_OUTER](
				file_id, n_azimuthal) = vphi;

			    n_azimuthal++;
			    /*
			    // debugging output of data
			    if(file_id == 0)
			    printf("Nphi = %d	Sigma = %.5e	T = %.5e
			    vr = %.5e	vphi = %.5e\n", n_azimuthal,
			    sigma*units::surface_density.get_cgs_factor(),
			    T*units::temperature.get_cgs_factor(), vr, vphi);
				       */
			}
			if (n_azimuthal !=
			    data[t_data::PRESCRIBED_DENSITY_OUTER]
				.get_size_azimuthal()) {
			    die("Could not read full ring from file %s, only %d / %d lines found\n",
				path_string.c_str(), n_azimuthal,
				data[t_data::PRESCRIBED_DENSITY_OUTER]
				    .get_max_azimuthal());
			}
		    }
		}

		if (count_files.size() != (unsigned int)num_files) {
		    die("Only prescribed boundary files %d / %d files were loaded!\n");
		}

		std::sort(count_files.begin(), count_files.end());
		for (unsigned int i = 0; i < count_files.size(); ++i) {
		    if (i != (unsigned int)count_files[i]) {
			die("Prescribed boundary file number %d was not loaded correctly!\n");
		    }
		}
	    }
	}
    }
}

void apply_boundary_condition(t_data &data, const double current_time, const double dt, const bool final)
{
    // if this is the final boundary condition and damping is enable, do it
    if (final && parameters::damping) {
	damping(data, dt);

	if(ECC_GROWTH_MONITOR){
		quantities::calculate_disk_delta_ecc_peri(data, delta_ecc_damp, delta_peri_damp);
	}
    }

    // inner boundary
    switch (parameters::boundary_inner) {
    case parameters::boundary_condition_open:
	open_boundary_inner(data);
	break;
    case parameters::boundary_condition_reflecting:
	reflecting_boundary_inner(data);
	break;
    case parameters::boundary_condition_zero_gradient:
	zero_gradient_boundary_inner(data);
	break;
    case parameters::boundary_condition_boundary_layer:
	boundary_layer_inner_boundary(data);
	break;
    case parameters::boundary_condition_nonreflecting:
	NonReflectingBoundary_inner(data, &data[t_data::V_RADIAL],
				    &data[t_data::DENSITY],
				    &data[t_data::ENERGY]);
	break;
    case parameters::boundary_condition_viscous_outflow:
	viscous_outflow_boundary_inner(data);
	break;
    case parameters::boundary_condition_jibin_spreading_ring:
	jibin_boundary_inner(data);
	break;
    case parameters::boundary_condition_evanescent: // evanescent works only for
						    // inner and outer together
						    // until now
	if (parameters::boundary_outer ==
	    parameters::boundary_condition_evanescent) {
	    EvanescentBoundary(data, dt);
	} else {
	    logging::print_master(
		LOG_ERROR
		"Different EvanescentBoundary Parameters. Old .par File?\n");
	    die("inner/outer evanescent boundary not implemented yet");
	}
	break;
    case parameters::boundary_condition_keplerian:
	keplerian2d_boundary_inner(data);
	break;
    case parameters::boundary_condition_precribed_time_variable:
	die("Inner precribed time variable boundary condition is not implemented yet!\n");
	break;
    case parameters::boundary_condition_center_of_mass_initial:
	die("Inner initial center of mass boundary is not implemented yet!\n");
	break;
    }

    // outer boundary
    switch (parameters::boundary_outer) {
    case parameters::boundary_condition_open:
	open_boundary_outer(data);
	break;
    case parameters::boundary_condition_reflecting:
	reflecting_boundary_outer(data);
	break;
    case parameters::boundary_condition_center_of_mass_initial: {
	if (final) {
	    damping_initial_center_of_mass_outer(data, dt);
	}
	initial_center_of_mass_boundary(data);
	break;
    }
    case parameters::boundary_condition_zero_gradient:
	zero_gradient_boundary_outer(data);
	break;
    case parameters::boundary_condition_boundary_layer:
	boundary_layer_outer_boundary(data);
	break;
    case parameters::boundary_condition_nonreflecting:
	NonReflectingBoundary_outer(data, &data[t_data::V_RADIAL],
				    &data[t_data::DENSITY],
				    &data[t_data::ENERGY]);
	break;
    case parameters::boundary_condition_jibin_spreading_ring:
	jibin_boundary_outer(data);
	break;
    case parameters::boundary_condition_precribed_time_variable: {
	boundary_condition_precribed_time_variable_outer(
		data, &data[t_data::DENSITY], current_time);
    } break;
    case parameters::boundary_condition_viscous_outflow:
	die("outer viscous outflow boundary not implemented");
	break;
    case parameters::boundary_condition_evanescent:
	// EvanescentBoundary_outer(VRadial, VAzimuthal, Density, Energy, dt,
	// sys);
	break;
    case parameters::boundary_condition_keplerian:
	keplerian2d_boundary_outer(data);
	break;
    }

    /// d Omega / dr = 0 has really bad effects on the massflow test
    /// not recommended
    if (parameters::domegadr_zero) {
	/// Vphi_i / r_i - Vphi_(i-1) / r_(i-1) = 0
	if (CPU_Rank == CPU_Highest) {
		#pragma omp parallel for
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::V_AZIMUTHAL].get_max_azimuthal();
		 ++n_azimuthal) {
		// this is a work around as long as V_AZIMUTHAL is defined as a
		// vector
		const double R_N =
		    Rmed[data[t_data::V_AZIMUTHAL].get_max_radial()];
		const double R_Nm1 =
		    Rmed[data[t_data::V_AZIMUTHAL].get_max_radial() - 1];

		data[t_data::V_AZIMUTHAL](
		    data[t_data::V_AZIMUTHAL].get_max_radial(), n_azimuthal) =
		    R_N / R_Nm1 *
		    data[t_data::V_AZIMUTHAL](
			data[t_data::V_AZIMUTHAL].get_max_radial() - 1,
			n_azimuthal);
	    }
	}

	if (CPU_Rank == 0) {
		#pragma omp parallel for
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::V_AZIMUTHAL].get_max_azimuthal();
		 ++n_azimuthal) {
		// this is a work around as long as V_AZIMUTHAL is defined as a
		// vector
		const double R1 = Rmed[1];
		const double R0 = Rmed[0];

		data[t_data::V_AZIMUTHAL](0, n_azimuthal) =
		    R0 / R1 * data[t_data::V_AZIMUTHAL](1, n_azimuthal);
	    }
	}
    }

    if (OuterSourceMass) {
	ApplyOuterSourceMass(&data[t_data::DENSITY], &data[t_data::V_RADIAL]);
    }

    if (parameters::massoverflow) {
	boundary_conditions::mass_overflow_willy(data, nullptr, false);
    }
}

/**
	inner open boundary condition
*/
void open_boundary_inner(t_data &data)
{
    if (CPU_Rank != 0)
	return;

	#pragma omp parallel for
    for (unsigned int n_azimuthal = 0;
	 n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
	 ++n_azimuthal) {
	// copy first ring into ghost ring
	data[t_data::DENSITY](0, n_azimuthal) =
	    data[t_data::DENSITY](1, n_azimuthal);
	data[t_data::ENERGY](0, n_azimuthal) =
	    data[t_data::ENERGY](1, n_azimuthal);

	// set velocity to min(v[+1],0) (allow only outflow)
	if (data[t_data::V_RADIAL](2, n_azimuthal) > 0.0) {
	    data[t_data::V_RADIAL](1, n_azimuthal) = 0.0;
	    data[t_data::V_RADIAL](0, n_azimuthal) = 0.0;
	} else {
	    data[t_data::V_RADIAL](1, n_azimuthal) =
		data[t_data::V_RADIAL](2, n_azimuthal);
	    data[t_data::V_RADIAL](0, n_azimuthal) =
		data[t_data::V_RADIAL](2, n_azimuthal);
	}
    }
}

/**
	outer open boundary condition
 */
void open_boundary_outer(t_data &data)
{
    if (CPU_Rank != CPU_Highest)
	return;

	#pragma omp parallel for
    for (unsigned int n_azimuthal = 0;
	 n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
	 ++n_azimuthal) {
	// copy last ring into ghost ring
	data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial(),
			      n_azimuthal) =
	    data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial() - 1,
				  n_azimuthal);
	data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial(),
			     n_azimuthal) =
	    data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial() - 1,
				 n_azimuthal);

	// set velocity to min(v[+1],0) (allow only outflow)
	if (data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial() - 2,
				   n_azimuthal) < 0.0) {
	    data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial() - 1,
				   n_azimuthal) = 0.0;
	    data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial(),
				   n_azimuthal) = 0.0;
	} else {
	    data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial() - 1,
				   n_azimuthal) =
		data[t_data::V_RADIAL](
		    data[t_data::V_RADIAL].get_max_radial() - 2, n_azimuthal);
	    data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial(),
				   n_azimuthal) =
		data[t_data::V_RADIAL](
		    data[t_data::V_RADIAL].get_max_radial() - 2, n_azimuthal);
	}
    }
}

void zero_gradient_boundary_inner(t_data &data)
{
    if (CPU_Rank != 0)
	return;

	#pragma omp parallel for
    for (unsigned int n_azimuthal = 0;
	 n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
	 ++n_azimuthal) {
	// copy first ring into ghost ring
	data[t_data::DENSITY](0, n_azimuthal) =
	    data[t_data::DENSITY](1, n_azimuthal);
	data[t_data::ENERGY](0, n_azimuthal) =
	    data[t_data::ENERGY](1, n_azimuthal);

	data[t_data::V_RADIAL](1, n_azimuthal) =
	    data[t_data::V_RADIAL](2, n_azimuthal);

	data[t_data::V_RADIAL](0, n_azimuthal) =
	    data[t_data::V_RADIAL](2, n_azimuthal);
    }
}

void zero_gradient_boundary_outer(t_data &data)
{
    if (CPU_Rank != CPU_Highest)
	return;

	#pragma omp parallel for
    for (unsigned int n_azimuthal = 0;
	 n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
	 ++n_azimuthal) {
	// copy last ring into ghost ring
	data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial(),
			      n_azimuthal) =
	    data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial() - 1,
				  n_azimuthal);
	data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial(),
			     n_azimuthal) =
	    data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial() - 1,
				 n_azimuthal);

	data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial() - 1,
			       n_azimuthal) =
	    data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial() - 2,
				   n_azimuthal);
	data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial(),
			       n_azimuthal) =
	    data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial() - 2,
				   n_azimuthal);

	const double R_outer = Rmed[data[t_data::V_AZIMUTHAL].get_max_radial()];
	const double R_inner =
	    Rmed[data[t_data::V_AZIMUTHAL].get_max_radial() - 1];

	data[t_data::V_AZIMUTHAL](data[t_data::V_AZIMUTHAL].get_max_radial(),
				  n_azimuthal) =
	    std::sqrt(R_outer / R_inner) *
	    data[t_data::V_AZIMUTHAL](
		data[t_data::V_AZIMUTHAL].get_max_radial() - 1, n_azimuthal);
    }
}

/**
	inner reflecting boundary condition
*/
void reflecting_boundary_inner(t_data &data)
{
    if (CPU_Rank == 0) {
	#pragma omp parallel for
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
	     ++n_azimuthal) {
	    // copy first ring into ghost ring
	    data[t_data::DENSITY](0, n_azimuthal) =
		data[t_data::DENSITY](1, n_azimuthal);
	    data[t_data::ENERGY](0, n_azimuthal) =
		data[t_data::ENERGY](1, n_azimuthal);

	    data[t_data::V_RADIAL](1, n_azimuthal) = 0.0;
	    data[t_data::V_RADIAL](0, n_azimuthal) =
		-data[t_data::V_RADIAL](2, n_azimuthal);
	}
    }
}

/**
	outer boundary_condition_precribed_time_variable_outer
*/
void boundary_condition_precribed_time_variable_outer(t_data &data,
							  t_polargrid *densitystar, const double current_time)
{
    if (CPU_Rank == CPU_Highest) {
	const int n_radial = data[t_data::DENSITY].get_max_radial();
	const double T_bin =
	    data.get_planetary_system().get_planet(1).get_orbital_period();

	const double step_size = T_bin / (double)PRESCRIBED_TIME_SEGMENT_NUMBER;
	const double real_time = current_time / step_size;
	const int integer_time = (int)std::floor(real_time);
	const int time_id = integer_time % PRESCRIBED_TIME_SEGMENT_NUMBER;
	const int time_id_next = (time_id + 1) % PRESCRIBED_TIME_SEGMENT_NUMBER;
	const double percent_to_next_timestep =
	    real_time - double(integer_time);

	#pragma omp parallel for
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
	     ++n_azimuthal) {

	    const double vr =
		data[t_data::PRESCRIBED_V_RADIAL_OUTER](time_id, n_azimuthal);
	    const double vr_next = data[t_data::PRESCRIBED_V_RADIAL_OUTER](
		time_id_next, n_azimuthal);
	    const double vr_cell =
		vr + percent_to_next_timestep * (vr_next - vr);

	    if (vr < 0.0) // allow inflow
	    {
		const double sigma = data[t_data::PRESCRIBED_DENSITY_OUTER](
		    time_id, n_azimuthal);
		const double sigma_next =
		    data[t_data::PRESCRIBED_DENSITY_OUTER](time_id_next,
							   n_azimuthal);
		const double sigma_cell =
		    sigma + percent_to_next_timestep * (sigma_next - sigma);

		const double energy =
		    data[t_data::PRESCRIBED_ENERGY_OUTER](time_id, n_azimuthal);
		const double energy_next =
		    data[t_data::PRESCRIBED_ENERGY_OUTER](time_id_next,
							  n_azimuthal);
		const double energy_cell =
		    energy + percent_to_next_timestep * (energy_next - energy);

		const double vphi = data[t_data::PRESCRIBED_V_AZIMUTHAL_OUTER](
		    time_id, n_azimuthal);
		const double vphi_next =
		    data[t_data::PRESCRIBED_V_AZIMUTHAL_OUTER](time_id_next,
							       n_azimuthal);
		const double vphi_cell =
		    vphi + percent_to_next_timestep * (vphi_next - vphi);

		// copy interpolated values into outer ghost ring
		(*densitystar)(n_radial, n_azimuthal) = sigma_cell;
		data[t_data::ENERGY](n_radial, n_azimuthal) = energy_cell;
		data[t_data::V_RADIAL](n_radial, n_azimuthal) = vr_cell;
		data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) = vphi_cell;

		/*
		const double r = Rmed[n_radial];
		const double v_kep = sqrt(constants::G * hydro_center_mass / r);


		printf("Nphi = %d	dens = (%.3e	%.3e)	vr = (%.3e
		%.3e) vphi = (%.3e	%.3e	%.3e)\n", n_azimuthal,
			   data[t_data::PRESCRIBED_DENSITY_OUTER](time_id,
		n_azimuthal)*units::surface_density.get_cgs_factor(),
				data[t_data::DENSITY](n_radial-5,
		n_azimuthal)*units::surface_density.get_cgs_factor(),
				data[t_data::PRESCRIBED_V_RADIAL_OUTER](time_id,
		n_azimuthal), data[t_data::V_RADIAL](n_radial-5, n_azimuthal),
				data[t_data::PRESCRIBED_V_AZIMUTHAL_OUTER](time_id,
		n_azimuthal), data[t_data::V_AZIMUTHAL](n_radial-5,
		n_azimuthal), v_kep);*/
	    } else { // normal outflow
		// copy last ring into ghost ring
		(*densitystar)(data[t_data::ENERGY].get_max_radial(),
			       n_azimuthal) =
		    (*densitystar)(data[t_data::ENERGY].get_max_radial() - 1,
				   n_azimuthal);

		data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial(),
				     n_azimuthal) =
		    data[t_data::ENERGY](
			data[t_data::ENERGY].get_max_radial() - 1, n_azimuthal);

		data[t_data::V_RADIAL](
		    data[t_data::V_RADIAL].get_max_radial() - 1, n_azimuthal) =
		    data[t_data::V_RADIAL](
			data[t_data::V_RADIAL].get_max_radial() - 2,
			n_azimuthal);

		data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial(),
				       n_azimuthal) =
		    data[t_data::V_RADIAL](
			data[t_data::V_RADIAL].get_max_radial() - 2,
			n_azimuthal);
	    }
	}
    }
}

/**
	outer reflecting boundary condition
*/
void reflecting_boundary_outer(t_data &data)
{
    if (CPU_Rank == CPU_Highest) {
	#pragma omp parallel for
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
	     ++n_azimuthal) {
	    // copy last ring into ghost ring
	    data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial(),
				  n_azimuthal) =
		data[t_data::DENSITY](
		    data[t_data::DENSITY].get_max_radial() - 1, n_azimuthal);
	    data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial(),
				 n_azimuthal) =
		data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial() - 1,
				     n_azimuthal);

	    data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial() - 1,
				   n_azimuthal) = 0.0;
	    data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial(),
				   n_azimuthal) =
		-data[t_data::V_RADIAL](
		    data[t_data::V_RADIAL].get_max_radial() - 2, n_azimuthal);
	}
    }
}

void viscous_outflow_boundary_inner(t_data &data)
{
    if (CPU_Rank == 0) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::VISCOSITY].get_max_azimuthal();
	     ++n_azimuthal) {
	    data[t_data::DENSITY](0, n_azimuthal) =
		data[t_data::DENSITY](1, n_azimuthal);
	    data[t_data::ENERGY](0, n_azimuthal) =
		data[t_data::ENERGY](1, n_azimuthal);
	}

	const double s = parameters::viscous_outflow_speed;

	#pragma omp parallel for
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::V_RADIAL].get_max_azimuthal();
	     ++n_azimuthal) {

	    const double Nu0 = data[t_data::VISCOSITY](0, n_azimuthal);
	    const double Nu1 = data[t_data::VISCOSITY](1, n_azimuthal);
	    const double Nu = 0.5 * (Nu0 + Nu1);

	    // V_rad =  - 1.5 / r * Nu (Kley, Papaloizou and Ogilvie, 2008)
	    data[t_data::V_RADIAL](1, n_azimuthal) = -1.5 * s / Rinf[1] * Nu;
	    data[t_data::V_RADIAL](0, n_azimuthal) = -1.5 * s / Rinf[0] * Nu;
	}
    }
}

void damping_single_inner(t_polargrid &quantity, t_polargrid &quantity0,
			  double dt)
{
    // use the correct radius array corresponding to quantity
    t_radialarray &radius = quantity.is_scalar() ? Rb : Ra;
    double delta;
    double X, X0, Xnew;
    bool is_density = strcmp(quantity.get_name(), "dens") == 0;

    // is this CPU in the inner damping domain?
    if ((parameters::damping_inner_limit > 1.0) &&
	(radius[0] < RMIN * parameters::damping_inner_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMIN * parameters::damping_inner_limit), false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMIN * parameters::damping_inner_limit), true);
	}

	double tau = parameters::damping_time_factor * 2.0 * M_PI /
		     calculate_omega_kepler(RMIN);

	// Needed for OpenMP to work
	double &InnerWaveDampingPositive = MassDelta.InnerWaveDampingPositive;
	double &InnerWaveDampingNegative = MassDelta.InnerWaveDampingNegative;

	#pragma omp parallel for reduction (+ : InnerWaveDampingPositive, InnerWaveDampingNegative)
	for (unsigned int n_radial = 0; n_radial <= limit; ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMIN * parameters::damping_inner_limit) /
		    (RMIN - RMIN * parameters::damping_inner_limit),
		2);
	    double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
		X = quantity(n_radial, n_azimuthal);
		X0 = quantity0(n_radial, n_azimuthal);
		Xnew = (X - X0) * exp_factor + X0;
		quantity(n_radial, n_azimuthal) = Xnew;
		delta = Xnew - X;
		if (is_density) {
		    if (delta > 0) {
			sum_without_ghost_cells(InnerWaveDampingPositive,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(InnerWaveDampingNegative,
						delta * Surf[n_radial],
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
    double delta;
    double X, X0, Xnew;
    bool is_density = strcmp(quantity.get_name(), "dens") == 0;

    // is this CPU in the outer damping domain?
    if ((parameters::damping_outer_limit < 1.0) &&
	(radius[quantity.get_max_radial()] >
	 RMAX * parameters::damping_outer_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMAX * parameters::damping_outer_limit) + 1, false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMAX * parameters::damping_outer_limit) + 1, true);
	}

	double tau = parameters::damping_time_factor * 2.0 * M_PI /
		     calculate_omega_kepler(RMAX);

	// Needed for OpenMP to work
	double &OuterWaveDampingPositive = MassDelta.OuterWaveDampingPositive;
	double &OuterWaveDampingNegative = MassDelta.OuterWaveDampingNegative;

	#pragma omp parallel for reduction (+ : OuterWaveDampingPositive, OuterWaveDampingNegative)
	for (unsigned int n_radial = limit;
	     n_radial < quantity.get_size_radial(); ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMAX * parameters::damping_outer_limit) /
		    (RMAX - RMAX * parameters::damping_outer_limit),
		2);
	    double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
		X = quantity(n_radial, n_azimuthal);
		X0 = quantity0(n_radial, n_azimuthal);
		Xnew = (X - X0) * exp_factor + X0;
		quantity(n_radial, n_azimuthal) = Xnew;
		delta = Xnew - X;
		if (is_density) {
		    if (delta > 0) {
			sum_without_ghost_cells(OuterWaveDampingPositive,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(OuterWaveDampingNegative,
						delta * Surf[n_radial],
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
    bool is_density = strcmp(quantity.get_name(), "dens") == 0;

    // is this CPU in the inner damping domain?
    if ((parameters::damping_inner_limit > 1.0) &&
	(radius[0] < RMIN * parameters::damping_inner_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMIN * parameters::damping_inner_limit), false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMIN * parameters::damping_inner_limit), true);
	}

	double tau = parameters::damping_time_factor * 2.0 * M_PI /
		     calculate_omega_kepler(RMIN);

	// Needed for OpenMP to work
	double &InnerWaveDampingPositive = MassDelta.InnerWaveDampingPositive;
	double &InnerWaveDampingNegative = MassDelta.InnerWaveDampingNegative;

	#pragma omp parallel for reduction (+ : InnerWaveDampingPositive, InnerWaveDampingNegative)
	for (unsigned int n_radial = 0; n_radial <= limit; ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMIN * parameters::damping_inner_limit) /
		    (RMIN - RMIN * parameters::damping_inner_limit),
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
			sum_without_ghost_cells(InnerWaveDampingPositive,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(InnerWaveDampingNegative,
						delta * Surf[n_radial],
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

    bool is_density = strcmp(quantity.get_name(), "dens") == 0;

    // is this CPU in the outer damping domain?
    if ((parameters::damping_outer_limit < 1.0) &&
	(radius[quantity.get_max_radial()] >
	 RMAX * parameters::damping_outer_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMAX * parameters::damping_outer_limit) + 1, false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMAX * parameters::damping_outer_limit) + 1, true);
	}

	double tau = parameters::damping_time_factor * 2.0 * M_PI /
		     calculate_omega_kepler(RMAX);

	// Needed for OpenMP to work
	double &OuterWaveDampingPositive = MassDelta.OuterWaveDampingPositive;
	double &OuterWaveDampingNegative = MassDelta.OuterWaveDampingNegative;

	#pragma omp parallel for reduction (+ : OuterWaveDampingPositive, OuterWaveDampingNegative)
	for (unsigned int n_radial = limit;
	     n_radial < quantity.get_size_radial(); ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMAX * parameters::damping_outer_limit) /
		    (RMAX - RMAX * parameters::damping_outer_limit),
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
			sum_without_ghost_cells(OuterWaveDampingPositive,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(OuterWaveDampingNegative,
						delta * Surf[n_radial],
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
    double delta;
    double X, X0, Xnew;
    bool is_density = strcmp(quantity.get_name(), "dens") == 0;

    // is this CPU in the inner damping domain?
    if ((parameters::damping_inner_limit > 1.0) &&
	(radius[0] < RMIN * parameters::damping_inner_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMIN * parameters::damping_inner_limit), false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMIN * parameters::damping_inner_limit), true);
	}

	double tau = parameters::damping_time_factor * 2.0 * M_PI /
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

	// Needed for OpenMP to work
	double &InnerWaveDampingPositive = MassDelta.InnerWaveDampingPositive;
	double &InnerWaveDampingNegative = MassDelta.InnerWaveDampingNegative;

	#pragma omp parallel for reduction (+ : InnerWaveDampingPositive, InnerWaveDampingNegative)
	for (unsigned int n_radial = 0; n_radial <= limit; ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMIN * parameters::damping_inner_limit) /
		    (RMIN - RMIN * parameters::damping_inner_limit),
		2);
	    double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
		X = quantity(n_radial, n_azimuthal);
		X0 = quantity0(n_radial, 0);
		Xnew = (X - X0) * exp_factor + X0;
		quantity(n_radial, n_azimuthal) = Xnew;
		delta = Xnew - X;
		if (is_density) {
		    if (delta > 0) {
			sum_without_ghost_cells(InnerWaveDampingPositive,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(InnerWaveDampingNegative,
						delta * Surf[n_radial],
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
    if ((parameters::damping_inner_limit > 1.0) &&
	(rinf[0] < RMIN * parameters::damping_inner_limit)) {
	// find range

	const static unsigned int limit =
	    get_rinf_id(RMIN * parameters::damping_inner_limit);

	const double tau = parameters::damping_time_factor * 2.0 * M_PI /
			   calculate_omega_kepler(RMIN);

	const double s = parameters::viscous_outflow_speed;

	#pragma omp parallel for
	for (unsigned int n_radial = Zero_no_ghost; n_radial <= limit;
	     ++n_radial) {
	    double factor = std::pow(
		(rinf[n_radial] - RMIN * parameters::damping_inner_limit) /
		    (RMIN - RMIN * parameters::damping_inner_limit),
		2);
	    double exp_factor = std::exp(-dt * factor / tau);

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
    double delta;
    double X, X0, Xnew;
    bool is_density = strcmp(quantity.get_name(), "dens") == 0;

    // is this CPU in the outer damping domain?
    if ((parameters::damping_outer_limit < 1.0) &&
	(radius[quantity.get_max_radial()] >
	 RMAX * parameters::damping_outer_limit)) {
	// find range
	unsigned int limit;
	if (quantity.is_scalar()) {
	    limit = clamp_r_id_to_rmed_grid(
		get_rmed_id(RMAX * parameters::damping_outer_limit) + 1, false);
	} else {
	    limit = clamp_r_id_to_radii_grid(
		get_rinf_id(RMAX * parameters::damping_outer_limit) + 1, true);
	}

	double tau = parameters::damping_time_factor * 2.0 * M_PI /
		     calculate_omega_kepler(RMAX);

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
	double &OuterWaveDampingPositive = MassDelta.OuterWaveDampingPositive;
	double &OuterWaveDampingNegative = MassDelta.OuterWaveDampingNegative;

	#pragma omp parallel for reduction (+ : OuterWaveDampingPositive, OuterWaveDampingNegative)
	for (unsigned int n_radial = limit;
	     n_radial < quantity.get_size_radial(); ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMAX * parameters::damping_outer_limit) /
		    (RMAX - RMAX * parameters::damping_outer_limit),
		2);
	    double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < quantity.get_size_azimuthal(); ++n_azimuthal) {
		X = quantity(n_radial, n_azimuthal);
		X0 = quantity0(n_radial, 0);
		Xnew = (X - X0) * exp_factor + X0;
		quantity(n_radial, n_azimuthal) = Xnew;
		delta = Xnew - X;
		if (is_density) {
		    if (delta > 0) {
			sum_without_ghost_cells(OuterWaveDampingPositive,
						delta * Surf[n_radial],
						n_radial);
		    } else {
			sum_without_ghost_cells(OuterWaveDampingNegative,
						delta * Surf[n_radial],
						n_radial);
		    }
		}
	    }
	}
    }
}

void damping_initial_center_of_mass_outer(t_data &data, double dt)
{

    // use the correct radius array corresponding to quantity
    const t_radialarray &radius = Rinf;
    t_polargrid &vrad_arr = data[t_data::V_RADIAL];
    t_polargrid &vphi_arr = data[t_data::V_AZIMUTHAL];

    const unsigned int np = data.get_planetary_system().get_number_of_planets();
    const Pair com_pos = data.get_planetary_system().get_center_of_mass(np);
    const Pair com_vel =
	data.get_planetary_system().get_center_of_mass_velocity(np);
    const double com_mass = data.get_planetary_system().get_mass(np);

    // is this CPU in the outer damping domain?
    if ((parameters::damping_outer_limit < 1.0) &&
	(radius[vrad_arr.get_max_radial()] >
	 RMAX * parameters::damping_outer_limit)) {

	const unsigned int clamped_vrad_id = clamp_r_id_to_radii_grid(
	    get_rinf_id(RMAX * parameters::damping_outer_limit) + 1,
	    vrad_arr.is_vector());

	double tau = parameters::damping_time_factor * 2.0 * M_PI /
		     calculate_omega_kepler(RMAX);

	#pragma omp parallel for
	for (unsigned int n_radial = clamped_vrad_id;
	     n_radial < vrad_arr.get_size_radial(); ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMAX * parameters::damping_outer_limit) /
		    (RMAX - RMAX * parameters::damping_outer_limit),
		2);
	    double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < vrad_arr.get_size_azimuthal(); ++n_azimuthal) {

		const double phi = (double)n_azimuthal * dphi;
		const double rinf = radius[n_radial];

		const double cell_x = rinf * std::cos(phi);
		const double cell_y = rinf * std::sin(phi);

		// Position in center of mass frame
		const double x_com = cell_x - com_pos.x;
		const double y_com = cell_y - com_pos.y;
		const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

		// pressure support correction
		double corr;
		double vr_init;
		if (parameters::initialize_pure_keplerian) {
		    corr = 1.0;
		    vr_init = 0.0;
		} else {
		    corr = std::sqrt(
			1.0 - std::pow(ASPECTRATIO_REF, 2) *
				  std::pow(r_com, 2.0 * FLARINGINDEX) *
				  (1. + SIGMASLOPE - 2.0 * FLARINGINDEX));

		    const double v_k =
			std::sqrt(constants::G * com_mass / r_com);
		    const double h =
			ASPECTRATIO_REF * std::pow(r_com, FLARINGINDEX);
		    const double cs_iso = h * v_k;
		    const double H = h * r_com;
		    const double nu = ALPHAVISCOSITY * cs_iso * H;
		    vr_init = -3.0 * nu / r_com *
			      (-SIGMASLOPE + 2.0 * FLARINGINDEX + 1.0);
		}

		// Velocity in center of mass frame
		const double cell_vphi_com =
			std::sqrt(constants::G * com_mass / r_com) * corr;
		const double cell_vr_com = vr_init;

		const double cell_vx_com =
		    (cell_vr_com * x_com - cell_vphi_com * y_com) / r_com;
		const double cell_vy_com =
		    (cell_vr_com * y_com + cell_vphi_com * x_com) / r_com;

		// shift velocity from center of mass frame to primary frame
		const double cell_vx = cell_vx_com + com_vel.x;
		const double cell_vy = cell_vy_com + com_vel.y;

		const double vr0 = (cell_x * cell_vx + cell_y * cell_vy) / rinf;

		const double vr = vrad_arr(n_radial, n_azimuthal);
		const double vr_new = (vr - vr0) * exp_factor + vr0;
		vrad_arr(n_radial, n_azimuthal) = vr_new;
	    }
	}

	const unsigned int clamped_vphi_id = clamp_r_id_to_radii_grid(
	    get_rinf_id(RMAX * parameters::damping_outer_limit) + 1,
	    vphi_arr.is_vector());

	#pragma omp parallel for
	for (unsigned int n_radial = clamped_vphi_id;
	     n_radial < vphi_arr.get_size_radial(); ++n_radial) {
	    double factor = std::pow(
		(radius[n_radial] - RMAX * parameters::damping_outer_limit) /
		    (RMAX - RMAX * parameters::damping_outer_limit),
		2);
	    double exp_factor = std::exp(-dt * factor / tau);

	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < vphi_arr.get_size_azimuthal(); ++n_azimuthal) {

		const double phi = ((double)n_azimuthal - 0.5) * dphi;
		const double rmed = Rmed[n_radial];

		const double cell_x = rmed * std::cos(phi);
		const double cell_y = rmed * std::sin(phi);

		// Position in center of mass frame
		const double x_com = cell_x - com_pos.x;
		const double y_com = cell_y - com_pos.y;
		const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

		// pressure support correction
		double corr;
		double vr0;
		if (parameters::initialize_pure_keplerian) {
		    corr = 1.0;
		    vr0 = 0.0;
		} else {
		    corr = std::sqrt(
			1.0 - std::pow(ASPECTRATIO_REF, 2) *
				  std::pow(r_com, 2.0 * FLARINGINDEX) *
				  (1. + SIGMASLOPE - 2.0 * FLARINGINDEX));

		    const double v_k =
			std::sqrt(constants::G * com_mass / r_com);
		    const double h =
			ASPECTRATIO_REF * std::pow(r_com, FLARINGINDEX);
		    const double cs_iso = h * v_k;
		    const double H = h * r_com;
		    const double nu = ALPHAVISCOSITY * cs_iso * H;
		    vr0 = -3.0 * nu / r_com *
			  (-SIGMASLOPE + 2.0 * FLARINGINDEX + 1.0);
		}

		// Velocity in center of mass frame
		const double cell_vphi_com =
			std::sqrt(constants::G * com_mass / r_com) * corr;
		const double cell_vr_com = vr0;

		const double cell_vx_com =
		    (cell_vr_com * x_com - cell_vphi_com * y_com) / r_com;
		const double cell_vy_com =
		    (cell_vr_com * y_com + cell_vphi_com * x_com) / r_com;

		// shift velocity from center of mass frame to primary frame
		const double cell_vx = cell_vx_com + com_vel.x;
		const double cell_vy = cell_vy_com + com_vel.y;

		const double vp0 = (cell_x * cell_vy - cell_vx * cell_y) / rmed -
				OmegaFrame * rmed;

		const double vp = vphi_arr(n_radial, n_azimuthal);
		const double vp_new = (vp - vp0) * exp_factor + vp0;
		vphi_arr(n_radial, n_azimuthal) = vp_new;
	    }
	}

	if (parameters::Adiabatic){
	t_polargrid &energy = data[t_data::ENERGY];
	//t_polargrid &sigma = data[t_data::DENSITY];
	#pragma omp parallel for
	for (unsigned int nr = clamped_vphi_id;
		 nr < energy.get_size_radial(); ++nr) {
		double factor = std::pow(
		(radius[nr] - RMAX * parameters::damping_outer_limit) /
			(RMAX - RMAX * parameters::damping_outer_limit),
		2);
		double exp_factor = std::exp(-dt * factor / tau);

		for (unsigned int naz = 0;
		 naz < energy.get_size_azimuthal(); ++naz) {
	const double cell_x = (*CellCenterX)(nr, naz);
	const double cell_y = (*CellCenterY)(nr, naz);

	// Position in center of mass frame
	const double x_com = cell_x - com_pos.x;
	const double y_com = cell_y - com_pos.y;
	const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

	/// Initial profile temperature
	const double cell_energy_profile =
	1.0 / (ADIABATICINDEX - 1.0) * parameters::sigma0 *
	std::pow(ASPECTRATIO_REF, 2) *
	std::pow(r_com, -SIGMASLOPE - 1.0 + 2.0 * FLARINGINDEX) *
	constants::G * com_mass;
	/*
	const double cell_sigma = sigma(nr, naz);
	const double temperature_floor =
	parameters::minimum_temperature *
	units::temperature.get_inverse_cgs_factor();

	const double energy_floor = temperature_floor * cell_sigma /
				parameters::MU * constants::R /
				(ADIABATICINDEX - 1.0);
	const double cell_energy0 = std::max(cell_energy_profile, energy_floor);
				*/
	const double cell_energy0 = cell_energy_profile;

	const double cell_energy = energy(nr, naz);
	const double energy_new = (cell_energy - cell_energy0) * exp_factor + cell_energy0;

	energy(nr, naz)  = energy_new;
		}
	}
	}


	/*

	t_polargrid &sigma = data[t_data::DENSITY];
	for (unsigned int nr = clamped_vphi_id;
		 nr < sigma.get_size_radial(); ++nr) {
		double factor = std::pow(
		(radius[nr] - RMAX * parameters::damping_outer_limit) /
			(RMAX - RMAX * parameters::damping_outer_limit),
		2);
		double exp_factor = std::exp(-dt * factor / tau);

		for (unsigned int naz = 0;
		 naz < sigma.get_size_azimuthal(); ++naz) {
	const double cell_x = (*CellCenterX)(nr, naz);
	const double cell_y = (*CellCenterY)(nr, naz);

	// Position in center of mass frame
	const double x_com = cell_x - com_pos.x;
	const double y_com = cell_y - com_pos.y;
	const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

	const double cell_sigma0 =
	parameters::sigma0 *
	std::pow(r_com,
		 -SIGMASLOPE); // we assume the floor is not reached.

	const double cell_sigma = sigma(nr, naz);
	const double sigma_new = (cell_sigma - cell_sigma0) * exp_factor + cell_sigma0;

	sigma(nr, naz)  = sigma_new;
		}
	}*/

    }
}

/**
	damping of all selected quantities
*/
void damping(t_data &data, double dt)
{
    if (parameters::damping) {
	const unsigned int number_of_quantities_to_damp =
	    parameters::damping_vector.size();

	for (unsigned int i = 0; i < number_of_quantities_to_damp; ++i) {
	    const parameters::t_DampingType *damper =
		&parameters::damping_vector[i];
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

void mass_overflow(t_data &data, const double current_time)
{
    if (CPU_Rank != CPU_Highest) {
	return;
    }

    static double last_PhysicalTime = 0.0;

	double dt = current_time - last_PhysicalTime;

    if (dt == 0.0) {
	return;
    }

	last_PhysicalTime = current_time;

    // get location of binary star
    if (parameters::mof_planet + 1 >
	data.get_planetary_system().get_number_of_planets()) {
	logging::print_master(
	    LOG_ERROR
	    "Wrong Planet/Star for Mass Overflow specified! Old .par File?\n");
	die("Wrong Planet/Star for Mass Overflow specified! Old .par File?");
    }
    const t_planet &planet =
	data.get_planetary_system().get_planet(parameters::mof_planet);
    const double xplanet = planet.get_x();
    const double yplanet = planet.get_y();

    const double omega_planet = planet.get_omega();
    // get grid cell where binary star is nearest
    // atan2(y,x) is from -PI to PI
    double angle = std::atan2(yplanet, xplanet) * 0.5 * M_1_PI;
    if (angle < 0) {
	angle += 1.0;
    }

    const unsigned int Nrad = data[t_data::DENSITY].get_max_radial();
    const unsigned int Nphi = data[t_data::DENSITY].get_size_azimuthal();
    const double r_cell = Rmed[Nrad];

    const double vr_fraction = 0.002;
    const double vr_stream =
	-omega_planet * r_cell *
	vr_fraction; // radial inflow velocity is a fraction of v_Kepler, (e.g.
		     // vr_fraction = 0.002)
    const double vphi_stream =
	(omega_planet - OmegaFrame) *
	r_cell; // set angular velocity to zero (in co-rotating frame)

    const double mdot_cgs =
	parameters::mof_value * units::cgs_Msol / units::cgs_Year;

    const double mdot_code = mdot_cgs * units::time.get_cgs_factor() *
			     units::mass.get_inverse_cgs_factor();

    const double Sigma_stream = mdot_code * dt / Surf[Nrad - 2];

    int nearest_grid_cell =
	(int)((double)Nphi * angle + 0.5); // nearest gridcell to binary star
    nearest_grid_cell = nearest_grid_cell % Nphi;

    // Calculate sigma from temperature according to
    // Meyer & Meyer-Hofmeister (1983a) equation 17
    const double Porb =
	2.0 * M_PI / omega_planet * units::time.get_cgs_factor() / 3600.0;
    // cross section Q
    const double Q = 2.4e13 * parameters::mof_temperature * Porb * Porb;
    // stream radius W
    const double W = std::sqrt(Q / M_PI);
    // circumference circ
    const double circ = 2.0 * M_PI * r_cell * units::length.get_cgs_factor();

    double sigma = 2 * W / circ;
    int number_of_cells =
	(double)Nphi * 3.0 *
	sigma; // walk through 3 sigma, so we get 99.7% accuracy
    double sigmabar = Nphi * sigma;

    double check = 0.0;
    for (int i = -number_of_cells; i <= number_of_cells; i++) {

	const double t_ramp =
	    parameters::mof_rampingtime * planet.get_orbital_period();

	double ramp_factor;
	if (PhysicalTime < t_ramp) {
	    ramp_factor = std::pow(std::sin(PhysicalTime * M_PI_2 / t_ramp), 4);
	} else {
	    ramp_factor = 1.0;
	}

	// adapt gauss profile
	double weight_factor;
	int gridcell;
	if (number_of_cells == 0) {
	    gridcell = nearest_grid_cell;
	    weight_factor = 1.0;
	} else {
	    gridcell = (nearest_grid_cell + i + Nphi) % Nphi;
	    weight_factor = 1.0 / (sigmabar * std::sqrt(2.0 * M_PI)) *
			    std::exp(-1.0 / 2.0 * std::pow(i / sigmabar, 2));
	}
	check += weight_factor;

	int gridcell_r = gridcell + 1;
	gridcell_r = gridcell_r % Nphi;

	double dens = ramp_factor * weight_factor * Sigma_stream;
	if (dens < parameters::sigma0 * parameters::sigma_floor) {
	    dens = parameters::sigma0 * parameters::sigma_floor;
	}

	if (parameters::Adiabatic) {
	    const double T_stream =
		parameters::mof_temperature; // Stream Temperature in Kelvin
	    const double e_stream =
		T_stream * units::temperature.get_inverse_cgs_factor() * dens /
		parameters::MU * constants::R /
		(ADIABATICINDEX - 1.0); // energy density equivalent to T_stream

	    data[t_data::ENERGY](Nrad - 2, gridcell) = e_stream;
	}

	data[t_data::DENSITY](Nrad - 2, gridcell) += dens;

	data[t_data::V_RADIAL](Nrad - 2, gridcell) = vr_stream;
	data[t_data::V_RADIAL](Nrad - 1, gridcell) = vr_stream;
	data[t_data::V_AZIMUTHAL](Nrad - 2, gridcell) = vphi_stream;
	data[t_data::V_AZIMUTHAL](Nrad - 2, gridcell_r) = vphi_stream;

#ifndef NDEBUG
	logging::print(
	    LOG_VERBOSE
	    "dens %lE, WF %lE , angle %lf, nearest_grid_cell %i, mass_stream %lE, gridcell %i, Nphi %i , noc %i , i %i \n",
	    data[t_data::DENSITY](Nrad, gridcell), 1.0, angle, gridcell,
	    Sigma_stream, gridcell, Nrad, Nphi * 2 + 1, 0);
#endif
    }

    // check if mass overflow is constant
    if (!(check > 0.99) || !(check < 1.01)) {
	logging::print_master(LOG_ERROR "weight_factor %lf should be 0.997 \n",
			      check);
    }
}

/**
 *	Mass overflow function copied from RH2D
 *	added gaussian profile to stream
 */
void mass_overflow_willy(t_data &data, t_polargrid *densitystar, bool transport)
{
    if (CPU_Rank != CPU_Highest) {
	return;
    }

    // get location of binary star
    if (parameters::mof_planet + 1 >
	data.get_planetary_system().get_number_of_planets()) {
	logging::print_master(
	    LOG_ERROR
	    "Wrong Planet/Star for Mass Overflow specified! Old .par File?\n");
	die("Wrong Planet/Star for Mass Overflow specified! Old .par File?");
    }

    const t_planet &planet =
	data.get_planetary_system().get_planet(parameters::mof_planet);
    const double xplanet = planet.get_x();
    const double yplanet = planet.get_y();

    const double omega_planet = planet.get_omega();
    // get grid cell where binary star is nearest
    // atan2(y,x) is from -PI to PI
    double angle = std::atan2(yplanet, xplanet) * 0.5 * M_1_PI;
    if (angle < 0) {
	angle += 1.0;
    }

    const unsigned int Nrad = data[t_data::DENSITY].get_max_radial();
    const unsigned int Nphi = data[t_data::DENSITY].get_size_azimuthal();
    const double r_cell = Rmed[Nrad];

    const double vr_fraction = 0.002;
    const double vr_stream =
	-omega_planet * r_cell *
	vr_fraction; // radial inflow velocity is a fraction of v_Kepler, (e.g.
		     // vr_fraction = 0.002)
    const double vphi_stream =
	(omega_planet - OmegaFrame) *
	r_cell; // set angular velocity to zero (in co-rotating frame)

    const double stream_mdot_cgs =
	parameters::mof_value * units::cgs_Msol / units::cgs_Year;
    const double stream_mdot_code = stream_mdot_cgs *
				    units::mass.get_inverse_cgs_factor() *
				    units::time.get_cgs_factor();
    const double Sigma_stream =
	fabs(stream_mdot_code / (dphi * Rinf[Nrad] * vr_stream));

    const int nearest_grid_cell = ((int)((double)Nphi * angle + 0.5)) %
				  Nphi; // nearest gridcell to binary star

    // Calculate sigma from temperature according to
    // Meyer & Meyer-Hofmeister (1983a) equation 17
    const double Porb =
	2.0 * M_PI / omega_planet * units::time.get_cgs_factor() / 3600.0;
    // cross section Q
    const double Q = 2.4e13 * parameters::mof_temperature * Porb * Porb;
    // stream radius W
    const double W = std::sqrt(Q * M_1_PI);
    // circumference circ
    const double circ = 2.0 * M_PI * r_cell * units::length.get_cgs_factor();

    double sigma = 2.0 * W / circ;
    int number_of_cells =
	(double)Nphi * 3.0 *
	sigma; // walk through 3 sigma, so we get 99.7% accuracy
    double sigmabar = Nphi * sigma;

    const double t_ramp =
	parameters::mof_rampingtime * planet.get_orbital_period();

    double ramp_factor;
    if (PhysicalTime < t_ramp) {
	ramp_factor = std::pow(std::sin(PhysicalTime * M_PI_2 / t_ramp), 6);
    } else {
	ramp_factor = 1.0;
    }

    double check = 0.0;
    for (int i = -number_of_cells; i <= number_of_cells; i++) {

	// adapt gauss profile
	double weight_factor;
	int gridcell;
	if (number_of_cells == 0) {
	    gridcell = nearest_grid_cell;
	    weight_factor = 1.0;
	} else {
	    gridcell = (nearest_grid_cell + i + Nphi) % Nphi;
	    weight_factor = 1.0 / (sigmabar * std::sqrt(2.0 * M_PI)) *
			    std::exp(-1.0 / 2.0 * std::pow(i / sigmabar, 2));
	}
	check += weight_factor;

	int gridcell_r = gridcell + 1;
	gridcell_r = gridcell_r % Nphi;

	double dens = ramp_factor * weight_factor * Sigma_stream;
	if (dens < parameters::sigma0 * parameters::sigma_floor) {
	    dens = parameters::sigma0 * parameters::sigma_floor;
	}

	if (parameters::Adiabatic) {
	    const double T_stream =
		parameters::mof_temperature; // Stream Temperature in Kelvin
	    const double e_stream =
		T_stream * units::temperature.get_inverse_cgs_factor() * dens /
		parameters::MU * constants::R /
		(ADIABATICINDEX - 1.0); // energy density equivalent to T_stream

	    data[t_data::ENERGY](Nrad, gridcell) = e_stream;
	}

	if (transport && (densitystar != nullptr)) {
	    (*densitystar)(Nrad, gridcell) = dens;
	} else {
	    data[t_data::DENSITY](Nrad, gridcell) = dens;
	}

	data[t_data::V_RADIAL](Nrad, gridcell) = vr_stream;
	data[t_data::V_RADIAL](Nrad + 1, gridcell) = vr_stream;
	data[t_data::V_AZIMUTHAL](Nrad, gridcell) = vphi_stream;
	data[t_data::V_AZIMUTHAL](Nrad, gridcell_r) = vphi_stream;

#ifndef NDEBUG
	logging::print(
	    LOG_VERBOSE
	    "dens %lE, WF %lE , angle %lf, nearest_grid_cell %i, mass_stream %lE, gridcell %i, Nphi %i , noc %i , i %i \n",
	    data[t_data::DENSITY](Nrad, gridcell), 1.0, angle, gridcell,
	    Sigma_stream, gridcell, Nrad, Nphi * 2 + 1, 0);
#endif
    }

    // check if mass overflow is constant
    if (!(check > 0.99) || !(check < 1.01)) {
	logging::print_master(LOG_ERROR "weight_factor %lf should be 0.997 \n",
			      check);
    }
}

void apply_boundary_condition_temperature(t_data &data)
{
    if (CPU_Rank == 0) {
	#pragma omp parallel for
	for (unsigned int n_radial = 0; n_radial <= GHOSTCELLS_B; ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal();
		 ++n_azimuthal) {
		data[t_data::ENERGY](n_radial, n_azimuthal) =
		    data[t_data::TEMPERATURE](n_radial, n_azimuthal) *
		    data[t_data::DENSITY](n_radial, n_azimuthal) /
		    (ADIABATICINDEX - 1.0) / parameters::MU * constants::R;
	    }
	}

	switch (parameters::boundary_inner) {
	case parameters::boundary_condition_open:
	    open_boundary_inner(data);
	    break;
	case parameters::boundary_condition_reflecting:
	    reflecting_boundary_inner(data);
	    break;
	default:
	    logging::print_master(
		LOG_ERROR
		"Boundary condition at inner boundary not supported for temperature!\n");
	    break;
	}

	#pragma omp parallel for
	for (unsigned int n_radial = 0; n_radial <= GHOSTCELLS_B; ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal();
		 ++n_azimuthal) {
		data[t_data::TEMPERATURE](n_radial, n_azimuthal) =
		    data[t_data::ENERGY](n_radial, n_azimuthal) /
		    data[t_data::DENSITY](n_radial, n_azimuthal) *
		    (ADIABATICINDEX - 1.0) * parameters::MU / constants::R;
	    }
	}
    }

    if (CPU_Rank == CPU_Highest) {
	#pragma omp parallel for
	for (unsigned int n_radial =
		 data[t_data::ENERGY].get_max_radial() - GHOSTCELLS_B;
	     n_radial <= data[t_data::ENERGY].get_max_radial(); ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal();
		 ++n_azimuthal) {
		data[t_data::ENERGY](n_radial, n_azimuthal) =
		    data[t_data::TEMPERATURE](n_radial, n_azimuthal) *
		    data[t_data::DENSITY](n_radial, n_azimuthal) /
		    (ADIABATICINDEX - 1.0) / parameters::MU * constants::R;
	    }
	}

	// outer boundary
	switch (parameters::boundary_outer) {
	case parameters::boundary_condition_open:
	    open_boundary_outer(data);
	    break;
	case parameters::boundary_condition_reflecting:
	    reflecting_boundary_outer(data);
	    break;
	default:
	    logging::print_master(
		LOG_ERROR
		"Boundary condition at outer boundary not supported for temperature!\n");
	    break;
	}

	#pragma omp parallel for
	for (unsigned int n_radial =
		 data[t_data::ENERGY].get_max_radial() - GHOSTCELLS_B;
	     n_radial <= data[t_data::ENERGY].get_max_radial(); ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal();
		 ++n_azimuthal) {
		data[t_data::TEMPERATURE](n_radial, n_azimuthal) =
		    data[t_data::ENERGY](n_radial, n_azimuthal) /
		    data[t_data::DENSITY](n_radial, n_azimuthal) *
		    (ADIABATICINDEX - 1.0) * parameters::MU / constants::R;
	    }
	}
    }
}

/**
	Boundary conditions for calculation of the boundary layer (BL) starting
   here:
	TODO: Stellar radiative flux into disk via implicit routine!
*/

/**
	Inner boundary: zero gradient & fixed velocities
*/

void boundary_layer_inner_boundary(t_data &data)
{

    if (CPU_Rank != 0)
	return;

	#pragma omp parallel for
    for (unsigned int n_azimuthal = 0;
	 n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
	 ++n_azimuthal) {
	// zero gradient
	data[t_data::DENSITY](0, n_azimuthal) =
	    data[t_data::DENSITY](1, n_azimuthal);
	data[t_data::ENERGY](0, n_azimuthal) =
	    data[t_data::ENERGY](1, n_azimuthal);

	// set vrad to fraction of Keplerian velocity
	data[t_data::V_RADIAL](1, n_azimuthal) =
	    -1. * parameters::vrad_fraction_of_kepler * std::sqrt(1. / Ra[1]);
	data[t_data::V_RADIAL](0, n_azimuthal) =
	    data[t_data::V_RADIAL](1, n_azimuthal);

	// set vphi to stellar rotation rate
	data[t_data::V_AZIMUTHAL](0, n_azimuthal) =
	    parameters::stellar_rotation_rate * std::sqrt(1. / Rb[0]);
    }
}

/**
	Outer boundary: floating boundary conditions & pressure correction for
   Omega
*/

void boundary_layer_outer_boundary(t_data &data)
{

    if (CPU_Rank != CPU_Highest)
	return;

	#pragma omp parallel for
    for (unsigned int n_azimuthal = 0;
	 n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
	 ++n_azimuthal) {
	// floating BCs
	data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial(),
			      n_azimuthal) =
	    data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial() - 1,
				  n_azimuthal) *
	    std::sqrt(Ra[data[t_data::DENSITY].get_max_radial() - 1] /
		      Ra[data[t_data::DENSITY].get_max_radial()]);

	data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial(),
			     n_azimuthal) =
	    data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial() - 1,
				 n_azimuthal) *
	    std::pow(Ra[data[t_data::ENERGY].get_max_radial() - 1] /
			 Ra[data[t_data::ENERGY].get_max_radial()],
		     1.25);

	data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial() - 1,
			       n_azimuthal) =
	    -1. *
	    fabs(data[t_data::V_RADIAL](
		data[t_data::V_RADIAL].get_max_radial() - 2, n_azimuthal)) *
	    std::sqrt(Ra[data[t_data::V_RADIAL].get_max_radial() - 2] /
		      Ra[data[t_data::V_RADIAL].get_max_radial() - 1]);

	data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial(),
			       n_azimuthal) =
	    -1. *
	    fabs(data[t_data::V_RADIAL](
		data[t_data::V_RADIAL].get_max_radial() - 2, n_azimuthal)) *
	    std::sqrt(Ra[data[t_data::V_RADIAL].get_max_radial() - 2] /
		      Ra[data[t_data::V_RADIAL].get_max_radial()]);

	// Omega at outer boundary equals calculate_omega_kepler (plus leading
	// order pressure correction)
	data[t_data::V_AZIMUTHAL](data[t_data::V_AZIMUTHAL].get_max_radial(),
				  n_azimuthal) =
	    1. / std::sqrt(Rb[data[t_data::DENSITY].get_max_radial()]);
	// TODO: Include pressure correction, like in uphi[*jN]
	// = 1./sqrt(Rb[*jN]) +
	// 0.5/Sigma[*jN]*sqrt(pow3(Rb[*jN])*pow2(Rb[*jN]))*.5/Rb[*jN]*(P[*jN+1]-P[*jN-1])/DeltaRa[*jN+1];
    }
}

/**
	End BL Boundary Conditions!
*/

void keplerian2d_boundary_inner(t_data &data)
{
    for (unsigned int n_azimuthal = 0;
	 n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal();
	 ++n_azimuthal) {
	data[t_data::DENSITY](1, n_azimuthal) =
	    parameters::sigma0 * std::pow(Rmed[1], -SIGMASLOPE);
	data[t_data::DENSITY](0, n_azimuthal) =
	    parameters::sigma0 * std::pow(Rmed[0], -SIGMASLOPE);
	data[t_data::ENERGY](1, n_azimuthal) =
	    1.0 / (ADIABATICINDEX - 1.0) * parameters::sigma0 *
	    std::pow(ASPECTRATIO_REF, 2) *
	    std::pow(Rmed[1], -SIGMASLOPE - 1.0 + 2.0 * FLARINGINDEX);
	data[t_data::ENERGY](0, n_azimuthal) =
	    1.0 / (ADIABATICINDEX - 1.0) * parameters::sigma0 *
	    std::pow(ASPECTRATIO_REF, 2) *
	    std::pow(Rmed[0], -SIGMASLOPE - 1.0 + 2.0 * FLARINGINDEX);
	data[t_data::TEMPERATURE](1, n_azimuthal) =
	    data[t_data::ENERGY](1, n_azimuthal) /
	    data[t_data::DENSITY](1, n_azimuthal) * (ADIABATICINDEX - 1.0) *
	    parameters::MU * constants::R;
	data[t_data::TEMPERATURE](0, n_azimuthal) =
	    data[t_data::ENERGY](0, n_azimuthal) /
	    data[t_data::DENSITY](0, n_azimuthal) * (ADIABATICINDEX - 1.0) *
	    parameters::MU * constants::R;
	data[t_data::V_RADIAL](1, n_azimuthal) = 0.0;
	data[t_data::V_RADIAL](0, n_azimuthal) =
	    -data[t_data::V_RADIAL](2, n_azimuthal);
	data[t_data::V_AZIMUTHAL](1, n_azimuthal) =
	    std::sqrt(constants::G * hydro_center_mass / Rmed[1]);
	data[t_data::V_AZIMUTHAL](0, n_azimuthal) =
	    std::sqrt(constants::G * hydro_center_mass / Rmed[0]);
    }
}

void keplerian2d_boundary_outer(t_data &data)
{
	#pragma omp parallel for
    for (unsigned int n_azimuthal = 0;
	 n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal();
	 ++n_azimuthal) {
	data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial(),
			      n_azimuthal) =
	    parameters::sigma0 *
	    std::pow(Rmed[data[t_data::DENSITY].get_max_radial()], -SIGMASLOPE);
	data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial() - 1,
			      n_azimuthal) =
	    parameters::sigma0 *
	    std::pow(Rmed[data[t_data::DENSITY].get_max_radial() - 1],
		     -SIGMASLOPE);
	data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial(),
			     n_azimuthal) =
	    1.0 / (ADIABATICINDEX - 1.0) * parameters::sigma0 *
	    std::pow(ASPECTRATIO_REF, 2) *
	    std::pow(Rmed[data[t_data::DENSITY].get_max_radial()],
		     -SIGMASLOPE - 1.0 + 2.0 * FLARINGINDEX);
	data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial() - 1,
			     n_azimuthal) =
	    1.0 / (ADIABATICINDEX - 1.0) * parameters::sigma0 *
	    std::pow(ASPECTRATIO_REF, 2) *
	    std::pow(Rmed[data[t_data::DENSITY].get_max_radial() - 1],
		     -SIGMASLOPE - 1.0 + 2.0 * FLARINGINDEX);
	data[t_data::TEMPERATURE](data[t_data::TEMPERATURE].get_max_radial(),
				  n_azimuthal) =
	    data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial(),
				 n_azimuthal) /
	    data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial(),
				  n_azimuthal) *
	    (ADIABATICINDEX - 1.0) * parameters::MU * constants::R;
	data[t_data::TEMPERATURE](
	    data[t_data::TEMPERATURE].get_max_radial() - 1, n_azimuthal) =
	    data[t_data::ENERGY](data[t_data::ENERGY].get_max_radial() - 1,
				 n_azimuthal) /
	    data[t_data::DENSITY](data[t_data::DENSITY].get_max_radial() - 1,
				  n_azimuthal) *
	    (ADIABATICINDEX - 1.0) * parameters::MU * constants::R;
	data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial(),
			       n_azimuthal) =
	    -data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial() - 2,
				    n_azimuthal);
	data[t_data::V_RADIAL](data[t_data::V_RADIAL].get_max_radial() - 1,
			       n_azimuthal) = 0.0;
	data[t_data::V_AZIMUTHAL](data[t_data::V_AZIMUTHAL].get_max_radial(),
				  n_azimuthal) =
	    std::sqrt(constants::G * hydro_center_mass /
		      Rmed[data[t_data::DENSITY].get_max_radial()]);
	data[t_data::V_AZIMUTHAL](
	    data[t_data::V_AZIMUTHAL].get_max_radial() - 1, n_azimuthal) =
	    std::sqrt(constants::G * hydro_center_mass /
		      Rmed[data[t_data::DENSITY].get_max_radial() - 1]);
    }
}

/**
 * @brief initial_center_of_mass_boundary: sets the outer boundary
 *  to the initial profile in the center of mass and then shifts it to the
 * primary center. Lucas thinks this is important when simulating a circumbinary
 * disk when the coordination system is centered on the primary.
 * @param data
 */
void initial_center_of_mass_boundary(t_data &data)
{

    if (CPU_Rank != CPU_Highest)
	return;

    const unsigned int np = data.get_planetary_system().get_number_of_planets();
    const Pair com_pos = data.get_planetary_system().get_center_of_mass(np);
    const Pair com_vel =
	data.get_planetary_system().get_center_of_mass_velocity(np);
    const double com_mass = data.get_planetary_system().get_mass(np);

    auto &sigma = data[t_data::DENSITY];
    auto &energy = data[t_data::ENERGY];
    auto &vrad = data[t_data::V_RADIAL];
    auto &vaz = data[t_data::V_AZIMUTHAL];

    const unsigned int nr = data[t_data::DENSITY].get_max_radial();
	#pragma omp parallel for
    for (unsigned int naz = 0; naz <= data[t_data::DENSITY].get_max_azimuthal();
	 ++naz) {

	{ /// V_PHI
	    const double phi = ((double)naz - 0.5) * dphi;
	    const double rmed = Rmed[nr];

	    const double cell_x = rmed * std::cos(phi);
	    const double cell_y = rmed * std::sin(phi);

	    // Position in center of mass frame
	    const double x_com = cell_x - com_pos.x;
	    const double y_com = cell_y - com_pos.y;
	    const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

	    // pressure support correction
	    double corr;
	    double vr0;
	    if (parameters::initialize_pure_keplerian) {
		corr = 1.0;
		vr0 = 0.0;
	    } else {
		corr =
			std::sqrt(1.0 - std::pow(ASPECTRATIO_REF, 2) *
					std::pow(r_com, 2.0 * FLARINGINDEX) *
					(1. + SIGMASLOPE - 2.0 * FLARINGINDEX));

		const double v_k = std::sqrt(constants::G * com_mass / r_com);
		const double h =
		    ASPECTRATIO_REF * std::pow(r_com, FLARINGINDEX);
		const double cs_iso = h * v_k;
		const double H = h * r_com;
		const double nu = ALPHAVISCOSITY * cs_iso * H;
		vr0 = -3.0 * nu / r_com *
		      (-SIGMASLOPE + 2.0 * FLARINGINDEX + 1.0);
	    }

	    // Velocity in center of mass frame
	    const double cell_vphi_com =
		std::sqrt(constants::G * com_mass / r_com) * corr;
	    const double cell_vr_com = vr0;

	    const double cell_vx_com =
		(cell_vr_com * x_com - cell_vphi_com * y_com) / r_com;
	    const double cell_vy_com =
		(cell_vr_com * y_com + cell_vphi_com * x_com) / r_com;

	    // shift velocity from center of mass frame to primary frame
	    const double cell_vx = cell_vx_com + com_vel.x;
	    const double cell_vy = cell_vy_com + com_vel.y;

	    const double cell_vphi =
		(cell_x * cell_vy - cell_vx * cell_y) / rmed;
		vaz(nr, naz) = cell_vphi - OmegaFrame * rmed;
	} /// END V_PHI

	{ /// V_R
	    const double phi = (double)naz * dphi;
	    const double rinf = Rinf[nr];

	    const double cell_x = rinf * std::cos(phi);
	    const double cell_y = rinf * std::sin(phi);

	    // Position in center of mass frame
	    const double x_com = cell_x - com_pos.x;
	    const double y_com = cell_y - com_pos.y;
	    const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

	    // pressure support correction
	    double corr;
	    double vr0;
	    if (parameters::initialize_pure_keplerian) {
		corr = 1.0;
		vr0 = 0.0;
	    } else {
		corr =
			std::sqrt(1.0 - std::pow(ASPECTRATIO_REF, 2) *
					std::pow(r_com, 2.0 * FLARINGINDEX) *
					(1. + SIGMASLOPE - 2.0 * FLARINGINDEX));

		const double v_k = std::sqrt(constants::G * com_mass / r_com);
		const double h =
		    ASPECTRATIO_REF * std::pow(r_com, FLARINGINDEX);
		const double cs_iso = h * v_k;
		const double H = h * r_com;
		const double nu = ALPHAVISCOSITY * cs_iso * H;
		vr0 = -3.0 * nu / r_com *
		      (-SIGMASLOPE + 2.0 * FLARINGINDEX + 1.0);
	    }

	    // Velocity in center of mass frame
	    const double cell_vphi_com =
		std::sqrt(constants::G * com_mass / r_com) * corr;
	    const double cell_vr_com = vr0;

	    const double cell_vx_com =
		(cell_vr_com * x_com - cell_vphi_com * y_com) / r_com;
	    const double cell_vy_com =
		(cell_vr_com * y_com + cell_vphi_com * x_com) / r_com;

	    // shift velocity from center of mass frame to primary frame
	    const double cell_vx = cell_vx_com + com_vel.x;
	    const double cell_vy = cell_vy_com + com_vel.y;

	    const double cell_vr = (cell_x * cell_vx + cell_y * cell_vy) / rinf;
	    vrad(nr, naz) = cell_vr;
	} /// END V_R

	{ /// V_R GHOST CELL
	    const double phi = (double)naz * dphi;
	    const double rinf = Rsup[nr];

	    const double cell_x = rinf * std::cos(phi);
	    const double cell_y = rinf * std::sin(phi);

	    // Position in center of mass frame
	    const double x_com = cell_x - com_pos.x;
	    const double y_com = cell_y - com_pos.y;
	    const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

	    // pressure support correction
	    double corr;
	    double vr0;
	    if (parameters::initialize_pure_keplerian) {
		corr = 1.0;
		vr0 = 0.0;
	    } else {
		corr =
			std::sqrt(1.0 - std::pow(ASPECTRATIO_REF, 2) *
					std::pow(r_com, 2.0 * FLARINGINDEX) *
					(1. + SIGMASLOPE - 2.0 * FLARINGINDEX));

		const double v_k = std::sqrt(constants::G * com_mass / r_com);
		const double h =
		    ASPECTRATIO_REF * std::pow(r_com, FLARINGINDEX);
		const double cs_iso = h * v_k;
		const double H = h * r_com;
		const double nu = ALPHAVISCOSITY * cs_iso * H;
		vr0 = -3.0 * nu / r_com *
		      (-SIGMASLOPE + 2.0 * FLARINGINDEX + 1.0);
	    }

	    // Velocity in center of mass frame
	    const double cell_vphi_com =
		std::sqrt(constants::G * com_mass / r_com) * corr;
	    const double cell_vr_com = vr0;

	    const double cell_vx_com =
		(cell_vr_com * x_com - cell_vphi_com * y_com) / r_com;
	    const double cell_vy_com =
		(cell_vr_com * y_com + cell_vphi_com * x_com) / r_com;

	    // shift velocity from center of mass frame to primary frame
	    const double cell_vx = cell_vx_com + com_vel.x;
	    const double cell_vy = cell_vy_com + com_vel.y;

	    const double cell_vr = (cell_x * cell_vx + cell_y * cell_vy) / rinf;
	    vrad(nr + 1, naz) = cell_vr;
	} /// END V_R GHOST CELL

	{ /// DENSITY and ENERGY
	    const double cell_x = (*CellCenterX)(nr, naz);
	    const double cell_y = (*CellCenterY)(nr, naz);

	    // Position in center of mass frame
	    const double x_com = cell_x - com_pos.x;
	    const double y_com = cell_y - com_pos.y;
	    const double r_com = std::sqrt(x_com * x_com + y_com * y_com);

	    const double cell_sigma =
		parameters::sigma0 *
		std::pow(r_com,
			 -SIGMASLOPE); // we assume the floor is not reached.
	    sigma(nr, naz) = cell_sigma;

	    /// Initial profile temperature
	const double cell_energy =
	    1.0 / (ADIABATICINDEX - 1.0) * parameters::sigma0 *
	    std::pow(ASPECTRATIO_REF, 2) *
	    std::pow(r_com, -SIGMASLOPE - 1.0 + 2.0 * FLARINGINDEX) *
	    constants::G * com_mass;

	const double temperature_floor =
	    parameters::minimum_temperature *
	    units::temperature.get_inverse_cgs_factor();

	const double energy_floor = temperature_floor * cell_sigma /
				    parameters::MU * constants::R /
				    (ADIABATICINDEX - 1.0);

	energy(nr, naz) = std::max(cell_energy, energy_floor);

	    /// dP / dr = 0
		/// energy(nr, naz) = energy(nr - 1, naz);
	} /// END DENSITY and ENERGY
    }
}

/**
	for viscous spreading ring comparison simulations for Jibin
*/
void jibin_boundary_inner(t_data &data)
{
    if (CPU_Rank != 0)
	return;

    const double h = ASPECTRATIO_REF;
    const double p = SIGMASLOPE;
    const double q = 2.0 * FLARINGINDEX - 1.0;
    const double R = Rmed[0];
    const double OmegaK = 1.0 / (R * std::sqrt(R));
    const double corr = std::sqrt(1.0 + (p + q) * h * h);
    const double vaz = R * OmegaK * corr - R * OmegaFrame;

	#pragma omp parallel for
    for (unsigned int n_azimuthal = 0;
	 n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
	 ++n_azimuthal) {

	// copy first ring into ghost ring
	data[t_data::DENSITY](0, n_azimuthal) =
	    data[t_data::DENSITY](1, n_azimuthal);

	data[t_data::V_AZIMUTHAL](0, n_azimuthal) = vaz;

	if (data[t_data::V_RADIAL](2, n_azimuthal) <= 0.0) { // outflow
	    data[t_data::V_RADIAL](1, n_azimuthal) =
		data[t_data::V_RADIAL](2, n_azimuthal);
	    data[t_data::V_RADIAL](0, n_azimuthal) =
		data[t_data::V_RADIAL](2, n_azimuthal);
	} else { // reflective
	    data[t_data::V_RADIAL](1, n_azimuthal) =
		-data[t_data::V_RADIAL](2, n_azimuthal);
	    data[t_data::V_RADIAL](0, n_azimuthal) =
		-data[t_data::V_RADIAL](2, n_azimuthal);
	}
    }
}

/**
	for viscous spreading ring comparison simulations for Jibin
*/
void jibin_boundary_outer(t_data &data)
{
    if (CPU_Rank == CPU_Highest) {

	const double h = ASPECTRATIO_REF;
	const double p = SIGMASLOPE;
	const double q = 2.0 * FLARINGINDEX - 1.0;
	const double R = Rmed[data[t_data::V_AZIMUTHAL].get_max_radial()];
	const double OmegaK = 1.0 / (R * std::sqrt(R));
	const double corr = std::sqrt(1.0 + (p + q) * h * h);
	const double vaz = R * OmegaK * corr - R * OmegaFrame;

	#pragma omp parallel for
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal();
	     ++n_azimuthal) {
	    data[t_data::V_AZIMUTHAL](
		data[t_data::V_AZIMUTHAL].get_max_radial(), n_azimuthal) = vaz;
	}
    }
}

} // namespace boundary_conditions
