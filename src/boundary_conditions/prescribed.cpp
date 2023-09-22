/**
	\file boundary_conditions.cpp

	All boundary conditions are handled in this file.
*/
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boundary_conditions.h"
#include "../global.h"
#include "../parameters.h"
#include "../constants.h"

#include <filesystem>

namespace boundary_conditions
{

int PRESCRIBED_TIME_SEGMENT_NUMBER;

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
	    if (parameters::PRESCRIBED_BOUNDARY_OUTER_FILE == "") {
		die("Outer prescribed time variable boundary condition is enabled but the supplied file folder is not found!\n");
	    } else {

		// TODO: naming convention might need adjustment
		const int Nphi =
		    data[t_data::PRESCRIBED_DENSITY_OUTER].get_size_azimuthal();
		std::string file_name_body = parameters::PRESCRIBED_BOUNDARY_OUTER_FILE + "/" + std::to_string(Nphi) + "shift";

		std::string file_name_test = file_name_body + "0.dat";
		if (!std::filesystem::exists(
			file_name_test)) {
		    die("Prescribed boundary file %s does not exist!\n",
			file_name_test.c_str());
		}

		// get number of files
		int num_files = 0;
		const std::filesystem::path File_Folder{
		    parameters::PRESCRIBED_BOUNDARY_OUTER_FILE};
		for (auto const &dir_entry :
		     std::filesystem::directory_iterator{
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
		     std::filesystem::directory_iterator{
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
				std::pow(RMAX, parameters::FLARINGINDEX);
			    const double P = sigma * Cs * Cs;
			    const double T =
				parameters::MU / constants::R * P / sigma;
			    const double energy = T * sigma / parameters::MU *
						  constants::R /
						  (parameters::ADIABATICINDEX - 1.0);

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

/**
	outer boundary_condition_precribed_time_variable_outer
*/
void boundary_condition_precribed_time_variable_outer(t_data &data,
							  t_polargrid *densitystar, const double current_time)
{
    if (CPU_Rank == CPU_Highest) {
	const int n_radial = data[t_data::SIGMA].get_max_radial();
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
	     n_azimuthal <= data[t_data::SIGMA].get_max_azimuthal();
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

} // namespace boundary_conditions
