#include "quantities.h"

#include "Theo.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "parameters.h"
#include "stress.h"
#include "util.h"
#include "frame_of_reference.h"
#include <math.h>
#include <mpi.h>
#include <vector>
#include "SourceEuler.h"
#include "pvte_law.h"
#include "simulation.h"
#include "gas_torques.h"
#include "viscosity/viscosity.h"


namespace quantities
{

void fill_alpha_array(t_data &data, unsigned int timestep,
				   bool force_update){

	static int last_timestep_calculated = -1;

	if (!force_update) {
	if (last_timestep_calculated == (int)timestep) {
		return;
	} else {
		last_timestep_calculated = timestep;
	}
	}

	const unsigned int Nr = data[t_data::SIGMA].get_size_radial();
	const unsigned int Nphi = data[t_data::SIGMA].get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const double alpha = viscosity::get_alpha(nr, naz, data);
		data[t_data::ALPHA](nr, naz) = alpha;
	}
	}

}

/**
	Calculates total gas mass.
*/
double gas_total_mass(t_data &data, const double quantitiy_radius)
{
    double local_mass = 0.0;
    double global_mass = 0.0;

    // calculate mass of this process' cells
	#pragma omp parallel for reduction(+ : local_mass)
    for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size; ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::SIGMA].get_size_azimuthal();
	     ++n_azimuthal) {
	    if (Rmed[n_radial] <= quantitiy_radius) {
		local_mass +=
		    Surf[n_radial] * data[t_data::SIGMA](n_radial, n_azimuthal);
	    }
	}
    }

	MPI_Allreduce(&local_mass, &global_mass, 1, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);

    return global_mass;
}

/**
 * @brief gas_quantity_reduce
 Warning: Reduce, only CPU 0 has the correct value
 */
double gas_quantity_reduce(const t_polargrid& arr, const double quantitiy_radius)
{

	double global_reduced_quantity = 0.0;
	double local_reduced_quantity = 0.0;

	// Loop thru all cells excluding GHOSTCELLS & CPUOVERLAP cells (otherwise
	// they would be included twice!)
	#pragma omp parallel for reduction(+ : local_reduced_quantity)
	for (unsigned int nr = radial_first_active; nr < radial_active_size; ++nr) {
	for (unsigned int naz = 0; naz < arr.get_size_azimuthal(); ++naz) {
		// eccentricity and semi major axis weighted with cellmass
		if (Rmed[nr] <= quantitiy_radius) {
		local_reduced_quantity += arr(nr, naz) ;
		}
	}
	}

	// synchronize threads
	MPI_Reduce(&local_reduced_quantity, &global_reduced_quantity, 1, MPI_DOUBLE, MPI_SUM, 0,
		  MPI_COMM_WORLD);

	return global_reduced_quantity;
}



double gas_allreduce_mass_average(t_data &data, const t_polargrid& arr, const double quantitiy_radius)
{

	const t_polargrid& sigma = data[t_data::SIGMA];

	double local_mass = 0.0;
	double global_mass = 0.0;

	double global_reduced_quantity = 0.0;
	double local_reduced_quantity = 0.0;

	// Loop thru all cells excluding GHOSTCELLS & CPUOVERLAP cells (otherwise
	// they would be included twice!)
	#pragma omp parallel for reduction(+ : local_reduced_quantity, local_mass)
	for (unsigned int nr = radial_first_active; nr < radial_active_size; ++nr) {
	for (unsigned int naz = 0; naz < arr.get_size_azimuthal(); ++naz) {
		// eccentricity and semi major axis weighted with cellmass
		if (Rmed[nr] <= quantitiy_radius) {
		const double cell_mass = sigma(nr, naz) * Surf[nr];
		local_mass += cell_mass;
		local_reduced_quantity += arr(nr, naz) * cell_mass;
		}
	}
	}

	MPI_Allreduce(&local_mass, &global_mass, 1, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);

	// synchronize threads
	MPI_Allreduce(&local_reduced_quantity, &global_reduced_quantity, 1, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);

	global_reduced_quantity /= global_mass;
	return global_reduced_quantity;

}


double gas_reduce_mass_average(t_data &data, const t_polargrid& arr, const double quantitiy_radius)
{

	const t_polargrid& sigma = data[t_data::SIGMA];

	double local_mass = 0.0;
	double global_mass = 0.0;

	double local_reduced_quantity = 0.0;
	double global_reduced_quantity = 0.0;

	// Loop thru all cells excluding GHOSTCELLS & CPUOVERLAP cells (otherwise
	// they would be included twice!)
	#pragma omp parallel for reduction(+ : local_reduced_quantity, local_mass)
	for (unsigned int nr = radial_first_active; nr < radial_active_size; ++nr) {
	for (unsigned int naz = 0; naz < arr.get_size_azimuthal(); ++naz) {
		// eccentricity and semi major axis weighted with cellmass
		if (Rmed[nr] <= quantitiy_radius) {
		const double cell_mass = sigma(nr, naz) * Surf[nr];
		local_mass += cell_mass;
		local_reduced_quantity += arr(nr, naz) * cell_mass;
		}
	}
	}

	MPI_Reduce(&local_mass, &global_mass, 1, MPI_DOUBLE, MPI_SUM, 0,
		  MPI_COMM_WORLD);

	// synchronize threads
	MPI_Reduce(&local_reduced_quantity, &global_reduced_quantity, 1, MPI_DOUBLE, MPI_SUM,
		  0, MPI_COMM_WORLD);

	if(CPU_Master && global_mass > 0.0){
	global_reduced_quantity /= global_mass;
	return global_reduced_quantity;
	} else {
		return 0.0;
	}
}


/**
 * @brief gas_disk_radius
 * @param data
 * @return
 */
double gas_disk_radius(t_data &data, const double total_mass)
{

    const unsigned int local_array_start = Zero_or_active;
    const unsigned int local_array_end = Max_or_active;
	static const unsigned int send_size = local_array_end - local_array_start;

	static std::vector<double> local_mass(send_size);

	#pragma omp parallel for
    for (unsigned int n_radial = local_array_start; n_radial < local_array_end;
	 ++n_radial) {
	local_mass[n_radial - local_array_start] = 0.0;
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::SIGMA].get_size_azimuthal();
	     ++n_azimuthal) {
	    local_mass[n_radial - local_array_start] +=
		Surf[n_radial] * data[t_data::SIGMA](n_radial, n_azimuthal);
	}
    }

    double radius = 0.0;
    double current_mass = 0.0;

    MPI_Gatherv(&local_mass[0], send_size, MPI_DOUBLE, GLOBAL_bufarray,
		RootNradialLocalSizes, RootNradialDisplacements, MPI_DOUBLE, 0,
		MPI_COMM_WORLD);

    if (CPU_Master) {
	int j = 0;
	for (int rank = 0; rank < CPU_Number; ++rank) {
	    int id = RootRanksOrdered[rank];
		for (int i = RootIMIN[id]; i <= RootIMAX[id]; ++i) {
		++j;
		current_mass += GLOBAL_bufarray[i];
		if (current_mass > parameters::disk_radius_mass_fraction * total_mass) {
		    radius = GlobalRmed[j];
		    goto found_radius; // break out of nested loop
		}
	    }
	}
    found_radius:
	void();
	}

    return radius;
}

/**
	Calculates angular gas momentum.
*/
double gas_angular_momentum(t_data &data, const double quantitiy_radius)
{
    double local_angular_momentum = 0.0;
    double global_angular_momentum = 0.0;

	#pragma omp parallel for reduction(+ : local_angular_momentum)
    for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size; ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::SIGMA].get_size_azimuthal();
	     ++n_azimuthal) {
	    if (Rmed[n_radial] <= quantitiy_radius) {
		local_angular_momentum +=
		    Surf[n_radial] * 0.5 *
		    (data[t_data::SIGMA](n_radial, n_azimuthal) +
		     data[t_data::SIGMA](
			 n_radial, n_azimuthal == 0
				       ? data[t_data::SIGMA].get_max_azimuthal()
				       : n_azimuthal - 1)) *
		    Rmed[n_radial] *
		    (data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) +
		     refframe::OmegaFrame * Rmed[n_radial]);
	    }
	    // local_angular_momentum +=
	    // Surf[n_radial]*data[t_data::DENSITY](n_radial,n_azimuthal)*Rmed[n_radial]*(0.5*(data[t_data::V_AZIMUTHAL](n_radial,n_azimuthal)+data[t_data::V_AZIMUTHAL](n_radial,n_azimuthal
	    // == data[t_data::V_AZIMUTHAL].get_max_azimuthal() ? 0 :
	    // n_azimuthal+1))+frame_of_reference::OmegaFrame*Rmed[n_radial]);
	}
    }

    MPI_Allreduce(&local_angular_momentum, &global_angular_momentum, 1,
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return global_angular_momentum;
}

/**
	Calculates gas internal energy
*/
double gas_internal_energy(t_data &data, const double quantitiy_radius)
{
    double local_internal_energy = 0.0;
    double global_internal_energy = 0.0;

	#pragma omp parallel for reduction(+ : local_internal_energy)
    for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size; ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::ENERGY].get_size_azimuthal();
	     ++n_azimuthal) {
	    if (Rmed[n_radial] <= quantitiy_radius) {
		local_internal_energy +=
		    Surf[n_radial] *
		    data[t_data::ENERGY](n_radial, n_azimuthal);
	    }
	}
    }

    MPI_Reduce(&local_internal_energy, &global_internal_energy, 1, MPI_DOUBLE,
	       MPI_SUM, 0, MPI_COMM_WORLD);

    return global_internal_energy;
}

double gas_viscous_dissipation(t_data &data, const double quantitiy_radius)
{
    double local_qplus = 0.0;
    double global_qplus = 0.0;

	#pragma omp parallel for reduction(+ : local_qplus)
    for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size; ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::QPLUS].get_size_azimuthal();
	     ++n_azimuthal) {
	    if (Rmed[n_radial] <= quantitiy_radius) {
		local_qplus +=
		    Surf[n_radial] * data[t_data::QPLUS](n_radial, n_azimuthal);
	    }
	}
    }

    MPI_Reduce(&local_qplus, &global_qplus, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);

    return global_qplus;
}

double gas_luminosity(t_data &data, const double quantitiy_radius)
{
    double local_qminus = 0.0;
    double global_qminus = 0.0;

	#pragma omp parallel for reduction(+ : local_qminus)
    for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size; ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::QMINUS].get_size_azimuthal();
	     ++n_azimuthal) {
	    if (Rmed[n_radial] <= quantitiy_radius) {
		local_qminus += Surf[n_radial] *
				data[t_data::QMINUS](n_radial, n_azimuthal);
	    }
	}
    }

    MPI_Reduce(&local_qminus, &global_qminus, 1, MPI_DOUBLE, MPI_SUM, 0,
	       MPI_COMM_WORLD);

    return global_qminus;
}

/**
	Calculates gas kinematic energy
*/
double gas_kinematic_energy(t_data &data, const double quantitiy_radius)
{
    double local_kinematic_energy = 0.0;
    double global_kinematic_energy = 0.0;

	#pragma omp parallel for reduction(+ : local_kinematic_energy)
    for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size; ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::SIGMA].get_size_azimuthal();
	     ++n_azimuthal) {
	    // centered-in-cell radial velocity
	    if (Rmed[n_radial] <= quantitiy_radius) {
		double v_radial_center =
		    (Rmed[n_radial] - Rinf[n_radial]) *
			data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) +
		    (Rsup[n_radial] - Rmed[n_radial]) *
			data[t_data::V_RADIAL](n_radial, n_azimuthal);
		v_radial_center /= (Rsup[n_radial] - Rinf[n_radial]);

		// centered-in-cell azimuthal velocity
		double v_azimuthal_center =
		    0.5 *
			(data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) +
			 data[t_data::V_AZIMUTHAL](
			     n_radial, n_azimuthal == data[t_data::V_AZIMUTHAL]
							  .get_max_azimuthal()
					   ? 0
					   : n_azimuthal + 1)) +
		    Rmed[n_radial] * refframe::OmegaFrame;

		local_kinematic_energy +=
		    0.5 * Surf[n_radial] *
		    data[t_data::SIGMA](n_radial, n_azimuthal) *
		    (std::pow(v_radial_center, 2) +
		     std::pow(v_azimuthal_center, 2));
	    }
	}
    }

    MPI_Reduce(&local_kinematic_energy, &global_kinematic_energy, 1, MPI_DOUBLE,
	       MPI_SUM, 0, MPI_COMM_WORLD);

    return global_kinematic_energy;
}

/**
	Calculates gas kinematic energy
*/
double gas_radial_kinematic_energy(t_data &data, const double quantitiy_radius)
{
    double local_kinematic_energy = 0.0;
    double global_kinematic_energy = 0.0;

	#pragma omp parallel for reduction(+ : local_kinematic_energy)
    for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size; ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::SIGMA].get_size_azimuthal();
	     ++n_azimuthal) {
	    if (Rmed[n_radial] <= quantitiy_radius) {
		// centered-in-cell radial velocity
		double v_radial_center =
		    (Rmed[n_radial] - Rinf[n_radial]) *
			data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) +
		    (Rsup[n_radial] - Rmed[n_radial]) *
			data[t_data::V_RADIAL](n_radial, n_azimuthal);
		v_radial_center /= (Rsup[n_radial] - Rinf[n_radial]);

		local_kinematic_energy +=
		    0.5 * Surf[n_radial] *
		    data[t_data::SIGMA](n_radial, n_azimuthal) *
		    std::pow(v_radial_center, 2);
	    }
	}
    }

    MPI_Reduce(&local_kinematic_energy, &global_kinematic_energy, 1, MPI_DOUBLE,
	       MPI_SUM, 0, MPI_COMM_WORLD);

    return global_kinematic_energy;
}

/**
	Calculates gas kinematic energy
*/
double gas_azimuthal_kinematic_energy(t_data &data,
				      const double quantitiy_radius)
{
    double local_kinematic_energy = 0.0;
    double global_kinematic_energy = 0.0;

	#pragma omp parallel for reduction(+ : local_kinematic_energy)
    for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size; ++n_radial) {
	for (unsigned int n_azimuthal = 0;
	     n_azimuthal < data[t_data::SIGMA].get_size_azimuthal();
	     ++n_azimuthal) {
	    if (Rmed[n_radial] <= quantitiy_radius) {
		// centered-in-cell azimuthal velocity
		double v_azimuthal_center =
		    0.5 *
			(data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) +
			 data[t_data::V_AZIMUTHAL](
			     n_radial, n_azimuthal == data[t_data::V_AZIMUTHAL]
							  .get_max_azimuthal()
					   ? 0
					   : n_azimuthal + 1)) +
		    Rmed[n_radial] * refframe::OmegaFrame;

		local_kinematic_energy +=
		    0.5 * Surf[n_radial] *
		    data[t_data::SIGMA](n_radial, n_azimuthal) *
		    std::pow(v_azimuthal_center, 2);
	    }
	}
    }

    MPI_Reduce(&local_kinematic_energy, &global_kinematic_energy, 1, MPI_DOUBLE,
	       MPI_SUM, 0, MPI_COMM_WORLD);

    return global_kinematic_energy;
}

static void calculate_disk_ecc_peri_nbody_center(t_data &data, unsigned int timestep,
				   bool force_update)
{
	static int last_timestep_calculated = -1;

	const double num_nbody = data.get_planetary_system().get_number_of_planets();
	const double cms_mass = data.get_planetary_system().get_mass(num_nbody);
	const Pair cms_pos = data.get_planetary_system().get_center_of_mass(num_nbody);
	const Pair cms_vel = data.get_planetary_system().get_center_of_mass_velocity(num_nbody);

	if (!force_update) {
	if (last_timestep_calculated == (int)timestep) {
		return;
	} else {
		last_timestep_calculated = timestep;
	}
	}
	// calculations outside the loop for speedup
	double sinFrameAngle = std::sin(refframe::FrameAngle);
	double cosFrameAngle = std::cos(refframe::FrameAngle);

	const unsigned int Nr = data[t_data::SIGMA].get_size_radial();
	const unsigned int Nphi = data[t_data::SIGMA].get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const double total_mass =
		cms_mass + data[t_data::SIGMA](nr, naz) * Surf[nr];

		// location of the cell
		const double angle = (double)naz * dphi;
		const double r_x = Rmed[nr] * std::cos(angle) - cms_pos.x;
		const double r_y = Rmed[nr] * std::sin(angle) - cms_pos.y;
		const double dist = std::sqrt(r_x*r_x + r_y*r_y);

		// averaged velocities
		const double v_xmed =
		std::cos(angle) * 0.5 *
			(data[t_data::V_RADIAL](nr, naz) +
			 data[t_data::V_RADIAL](nr + 1, naz)) -
		std::sin(angle) *
			(0.5 *
			 (data[t_data::V_AZIMUTHAL](nr, naz) +
			  data[t_data::V_AZIMUTHAL](
				  nr, naz == data[t_data::V_AZIMUTHAL]
							   .get_max_azimuthal()
						? 0
						: naz + 1)) +
			 refframe::OmegaFrame * Rmed[nr]) - cms_vel.x;
		const double v_ymed =
		std::sin(angle) * 0.5 *
			(data[t_data::V_RADIAL](nr, naz) +
			 data[t_data::V_RADIAL](nr + 1, naz)) +
		std::cos(angle) *
			(0.5 *
			 (data[t_data::V_AZIMUTHAL](nr, naz) +
			  data[t_data::V_AZIMUTHAL](
				  nr, naz == data[t_data::V_AZIMUTHAL]
							   .get_max_azimuthal()
						? 0
						: naz + 1)) +
			 refframe::OmegaFrame * Rmed[nr]) - cms_vel.y;

		// specific angular momentum for each cell j = j*e_z
		const double j = r_x * v_ymed - r_y * v_xmed;
		// Runge-Lenz vector Ax = x*vy*vy-y*vx*vy-G*m*x/d;
		const double e_x =
		j * v_ymed / (constants::G * total_mass) - r_x / dist;
		const double e_y = -1.0 * j * v_xmed / (constants::G * total_mass) -
		  r_y / dist;

		data[t_data::ECCENTRICITY](nr, naz) =
		std::sqrt(std::pow(e_x, 2) + std::pow(e_y, 2));

		if (refframe::FrameAngle != 0.0) {
		// periastron grid is rotated to non-rotating coordinate system
		// to prevent phase jumps of atan2 in later transformations like
		// you would have had if you back-transform the output
		// periastron values
		data[t_data::PERIASTRON](nr, naz) =
			std::atan2(e_y * cosFrameAngle + e_x * sinFrameAngle,
				   e_x * cosFrameAngle - e_y * sinFrameAngle);
		} else {
		data[t_data::PERIASTRON](nr, naz) =
			std::atan2(e_y, e_x);
		}
	}
	}
}

static void calculate_disk_ecc_peri_hydro_center(t_data &data, unsigned int timestep,
			       bool force_update)
{
    static int last_timestep_calculated = -1;
	const double cms_mass = hydro_center_mass;
	const Pair cms_pos = data.get_planetary_system().get_hydro_frame_center_position();
	const Pair cms_vel = data.get_planetary_system().get_hydro_frame_center_velocity();

    if (!force_update) {
	if (last_timestep_calculated == (int)timestep) {
	    return;
	} else {
	    last_timestep_calculated = timestep;
	}
    }

    // calculations outside the loop for speedup
	const double sinFrameAngle = std::sin(refframe::FrameAngle);
	const double cosFrameAngle = std::cos(refframe::FrameAngle);

	const unsigned int Nr = data[t_data::SIGMA].get_size_radial();
	const unsigned int Nphi = data[t_data::SIGMA].get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const double total_mass =
		cms_mass +	data[t_data::SIGMA](nr, naz) * Surf[nr];

	    // location of the cell
		const double angle = (double)naz * dphi;
		const double r_x = Rmed[nr] * std::cos(angle) - cms_pos.x;
		const double r_y = Rmed[nr] * std::sin(angle) - cms_pos.y;
		const double dist = std::sqrt(r_x*r_x + r_y*r_y);


	    // averaged velocities
		const double v_xmed =
		std::cos(angle) * 0.5 *
			(data[t_data::V_RADIAL](nr, naz) +
			 data[t_data::V_RADIAL](nr + 1, naz)) -
		std::sin(angle) *
		    (0.5 *
			 (data[t_data::V_AZIMUTHAL](nr, naz) +
			  data[t_data::V_AZIMUTHAL](
				  nr, naz == data[t_data::V_AZIMUTHAL]
							   .get_max_azimuthal()
					    ? 0
						: naz + 1)) +
			 refframe::OmegaFrame * Rmed[nr]) - cms_vel.x;
		const double v_ymed =
		std::sin(angle) * 0.5 *
			(data[t_data::V_RADIAL](nr, naz) +
			 data[t_data::V_RADIAL](nr + 1, naz)) +
		std::cos(angle) *
		    (0.5 *
			 (data[t_data::V_AZIMUTHAL](nr, naz) +
			  data[t_data::V_AZIMUTHAL](
				  nr, 
				  naz == data[t_data::V_AZIMUTHAL].get_max_azimuthal() ? 0: naz + 1)) +
			 refframe::OmegaFrame * Rmed[nr]) - cms_vel.y;

	    // specific angular momentum for each cell j = j*e_z
		const double j = r_x * v_ymed - r_y * v_xmed;
	    // Runge-Lenz vector Ax = x*vy*vy-y*vx*vy-G*m*x/d;
		const double e_x =
		j * v_ymed / (constants::G * total_mass) - r_x / dist;
		const double e_y = -1.0 * j * v_xmed / (constants::G * total_mass) -
		  r_y / dist;

		data[t_data::ECCENTRICITY](nr, naz) =
		std::sqrt(std::pow(e_x, 2) + std::pow(e_y, 2));

	    if (refframe::FrameAngle != 0.0) {
		// periastron grid is rotated to non-rotating coordinate system
		// to prevent phase jumps of atan2 in later transformations like
		// you would have had if you back-transform the output
		// periastron values
		data[t_data::PERIASTRON](nr, naz) =
		    std::atan2(e_y * cosFrameAngle + e_x * sinFrameAngle,
			       e_x * cosFrameAngle - e_y * sinFrameAngle);
	    } else {
		data[t_data::PERIASTRON](nr, naz) =
		    std::atan2(e_y, e_x);
	    }
	}
    }
}

void calculate_disk_ecc_peri(t_data &data, unsigned int timestep,
				   bool force_update){
	if(parameters::n_bodies_for_hydroframe_center == 1){
		if(data.get_planetary_system().get_number_of_planets() > 1){
			// Binary has effects out to ~ 15 abin, if that is not inside the domain, compute ecc around primary
			if(data.get_planetary_system().get_planet(1).get_semi_major_axis()*15.0 > RMAX){
				calculate_disk_ecc_peri_hydro_center(data, timestep, force_update);
			} else {
				// We are looking at a circumbinary (or more Nbodies) disk
				calculate_disk_ecc_peri_nbody_center(data, timestep, force_update);
			}
		} else {
			// We only have a star, compute ecc around primary
			calculate_disk_ecc_peri_hydro_center(data, timestep, force_update);
		}
	} else {
		// If we have multiple objects as hydro center, always compute eccentricity around hydro center
		calculate_disk_ecc_peri_hydro_center(data, timestep, force_update);
	}
}


static void calculate_disk_ecc_vector_worker(t_data &data, const unsigned int num_planets_for_center)
{

	const double cms_mass = data.get_planetary_system().get_mass(num_planets_for_center);
	const Pair cms_pos = data.get_planetary_system().get_center_of_mass(num_planets_for_center);
	const Pair cms_vel = data.get_planetary_system().get_center_of_mass_velocity(num_planets_for_center);

	const t_polargrid &density = data[t_data::SIGMA];
	const t_polargrid &vr = data[t_data::V_RADIAL];
	const t_polargrid &vphi = data[t_data::V_AZIMUTHAL];

	const double sinFrameAngle = std::sin(refframe::FrameAngle);
	const double cosFrameAngle = std::cos(refframe::FrameAngle);

	const unsigned int Nr = density.get_size_radial();
	const unsigned int Naz = density.get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Naz; ++naz) {

		const unsigned int naz_next = (naz == Naz - 1 ? 0 : naz + 1);
		const double total_mass = cms_mass + density(nr, naz) * Surf[nr];

		// location of the cell
		const double angle = (double)naz * dphi;
		const double r_x = Rmed[nr] * std::cos(angle) - cms_pos.x;
		const double r_y = Rmed[nr] * std::sin(angle) - cms_pos.y;
		const double dist = std::sqrt(r_x*r_x + r_y*r_y);

		// averaged velocities
		const double v_xmed =
		std::cos(angle) * 0.5 *	(vr(nr, naz) + vr(nr + 1, naz))
		- std::sin(angle) *	(0.5 * (vphi(nr, naz) + vphi(nr, naz_next))
							 + refframe::OmegaFrame * Rmed[nr]) - cms_vel.x;
		const double v_ymed =
		std::sin(angle) * 0.5 *	(vr(nr, naz) + vr(nr + 1, naz)) +
		std::cos(angle) * (0.5 * (vphi(nr, naz) + vphi(nr, naz_next)) +
			 refframe::OmegaFrame * Rmed[nr]) - cms_vel.y;

		// specific angular momentum for each cell j = j*e_z
		const double j = r_x * v_ymed - r_y * v_xmed;
		// Runge-Lenz vector Ax = x*vy*vy-y*vx*vy-G*m*x/d;
		const double e_x =
		j * v_ymed / (constants::G * total_mass) - r_x / dist;
		const double e_y = -1.0 * j * v_xmed / (constants::G * total_mass) -
		  r_y / dist;

		data[t_data::ECCENTRICITY](nr, naz) = std::sqrt(std::pow(e_x, 2) + std::pow(e_y, 2));

		if (refframe::FrameAngle != 0.0) {
		// periastron grid is rotated to non-rotating coordinate system
		// to prevent phase jumps of atan2 in later transformations like
		// you would have had if you back-transform the output
		// periastron values
		data[t_data::PERIASTRON](nr, naz) =
			std::atan2(e_y * cosFrameAngle + e_x * sinFrameAngle,
				   e_x * cosFrameAngle - e_y * sinFrameAngle);
		} else {
		data[t_data::PERIASTRON](nr, naz) =
			std::atan2(e_y, e_x);
		}
	}

	}
}

void calculate_disk_ecc_vector(t_data &data){
	int n_bodies_for_cms;
	if(parameters::n_bodies_for_hydroframe_center == 1){
		if(data.get_planetary_system().get_number_of_planets() > 1){
			// Binary has effects out to ~ 15 abin, if that is not inside the domain, compute ecc around primary
			if(data.get_planetary_system().get_planet(1).get_semi_major_axis()*15.0 > RMAX){
				n_bodies_for_cms = parameters::n_bodies_for_hydroframe_center;
			} else {
				// We are looking at a circumbinary (or more Nbodies) disk
				n_bodies_for_cms = data.get_planetary_system().get_number_of_planets();
			}
		} else {
			// We only have a star, compute ecc around primary
			n_bodies_for_cms = parameters::n_bodies_for_hydroframe_center;
		}
	} else {
		// If we have multiple objects as hydro center, always compute eccentricity around hydro center
		n_bodies_for_cms = parameters::n_bodies_for_hydroframe_center;
	}
	calculate_disk_ecc_vector_worker(data, n_bodies_for_cms);

}

void state_disk_ecc_peri_calculation_center(t_data &data){
	if(parameters::n_bodies_for_hydroframe_center == 1){
		if(data.get_planetary_system().get_number_of_planets() > 1){
			// Binary has effects out to ~ 15 abin, if that is not inside the domain, compute ecc around primary
			if(data.get_planetary_system().get_planet(1).get_semi_major_axis()*15.0 > RMAX){
				logging::print_master(LOG_INFO "Computing eccentricity / pericenter with respect to the hydro frame (primary) center!\n");
			} else {
				logging::print_master(LOG_INFO "Computing eccentricity / pericenter with respect to the center of mass of the Nbody system!\n");
			}
		} else {
			// We only have a star, compute ecc around primary
			logging::print_master(LOG_INFO "Computing eccentricity / pericenter with respect to the hydro frame (primary) center!\n");
		}
	} else {
		// If we have multiple objects as hydro center, always compute eccentricity around hydro center
		logging::print_master(LOG_INFO "Computing eccentricity / pericenter with respect to the hydro frame (Nbody) center!\n");
	}
}

void calculate_disk_delta_ecc_peri(t_data &data, double &dEcc, double &dPer)
{
	// ecc holds the current eccentricity
	t_polargrid &ecc_new = data[t_data::ECCENTRICITY];
	t_polargrid &P_new = data[t_data::PERIASTRON];

	// compute new eccentricity into ecc
	calculate_disk_ecc_vector(data);

	const double e_new = quantities::gas_reduce_mass_average(data, ecc_new, parameters::quantities_radius_limit);
	const double peri_new = quantities::gas_reduce_mass_average(data, P_new, parameters::quantities_radius_limit);

	const double de = e_new - ecc_old;
	double dp = peri_new - peri_old;

	if(dp < -M_PI){
		dp += 2.0*M_PI;
	}
	if(dp > M_PI){
		dp -= 2.0*M_PI;
	}

	ecc_old = e_new;
	peri_old = peri_new;

	dEcc += de;
	dPer += dp;
	return;
}


/**
	compute alpha gravitational

	alpha(R) = |d ln Omega/d ln R|^-1 (Tgrav)/(Sigma cs^2)
*/
void calculate_alpha_grav(t_data &data, unsigned int timestep,
			  bool force_update)
{
    static int last_timestep_calculated = -1;

    if (parameters::self_gravity != true)
	return;

    if (!force_update) {
	if (last_timestep_calculated == (int)timestep) {
	    return;
	} else {
	    last_timestep_calculated = timestep;
	}
    }

    stress::calculate_gravitational_stress(data);

	const unsigned int Nr = data[t_data::ALPHA_GRAV].get_size_radial();
	const unsigned int Nphi = data[t_data::ALPHA_GRAV].get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
	    /*data[t_data::ALPHA_GRAV](n_radial, n_azimuthal) = -2.0/3.0 *
	     * (data[t_data::T_GRAVITATIONAL](n_radial,n_azimuthal)+data[t_data::T_REYNOLDS](n_radial,
	     * n_azimuthal))/(data[t_data::DENSITY](n_radial,n_azimuthal)*pow2(data[t_data::SOUNDSPEED](n_radial,
	     * n_azimuthal)));*/
		data[t_data::ALPHA_GRAV](nr, naz) =
		2.0 / 3.0 *
		data[t_data::T_GRAVITATIONAL](nr, naz) /
		(data[t_data::SIGMA](nr, naz) *
		 std::pow(data[t_data::SOUNDSPEED](nr, naz), 2));
	}
    }
}

void calculate_alpha_grav_mean_sumup(t_data &data, unsigned int timestep,
				     double dt)
{
    calculate_alpha_grav(data, timestep, true);

	const unsigned int Nr = data[t_data::ALPHA_GRAV_MEAN].get_size_radial();
	const unsigned int Nphi = data[t_data::ALPHA_GRAV_MEAN].get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		data[t_data::ALPHA_GRAV_MEAN](nr, naz) +=
		data[t_data::ALPHA_GRAV](nr, naz) * dt;
	}
    }
}

/**
	compute alpha Reynolds

	alpha(R) = |d ln Omega/d ln R|^-1 (Trey)/(Sigma cs^2)
*/
void calculate_alpha_reynolds(t_data &data, unsigned int timestep,
			      bool force_update)
{
    static int last_timestep_calculated = -1;

    stress::calculate_Reynolds_stress(data);

    if (!force_update) {
	if (last_timestep_calculated == (int)timestep) {
	    return;
	} else {
	    last_timestep_calculated = timestep;
	}
    }

	const unsigned int Nr = data[t_data::ALPHA_REYNOLDS].get_size_radial();
	const unsigned int Nphi = data[t_data::ALPHA_REYNOLDS].get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		data[t_data::ALPHA_REYNOLDS](nr, naz) =
		2.0 / 3.0 * data[t_data::T_REYNOLDS](nr, naz) /
		(data[t_data::SIGMA](nr, naz) *
		 std::pow(data[t_data::SOUNDSPEED](nr, naz), 2));
	}
    }
}

void calculate_alpha_reynolds_mean_sumup(t_data &data, unsigned int timestep,
					 double dt)
{
    calculate_alpha_reynolds(data, timestep, true);

	const unsigned int Nr = data[t_data::ALPHA_REYNOLDS_MEAN].get_size_radial();
	const unsigned int Nphi = data[t_data::ALPHA_REYNOLDS_MEAN].get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		data[t_data::ALPHA_REYNOLDS_MEAN](nr, naz) +=
		data[t_data::ALPHA_REYNOLDS](nr, naz) * dt;
	}
    }
}

/**
	Calculates Toomre Q parameter
*/
void calculate_toomre(t_data &data, unsigned int /* timestep */,
		      bool /* force_update */)
{
    double kappa;

	const unsigned int Nr = data[t_data::TOOMRE].get_size_radial();
	const unsigned int Nphi = data[t_data::TOOMRE].get_size_azimuthal();

	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 1; nr < Nr; ++nr) {
	for (unsigned int naz = 0; naz < Nphi; ++naz) {
	    // kappa^2 = 1/r^3 d((r^2 Omega)^2)/dr = 1/r^3 d((r*v_phi)^2)/dr
	    kappa = std::sqrt(std::fabs(
		std::pow(InvRmed[nr], 3) *
		(std::pow(data[t_data::V_AZIMUTHAL](nr, naz) * Rmed[nr],  2) -
		 std::pow(data[t_data::V_AZIMUTHAL](nr - 1, naz) * Rmed[nr - 1], 2)) *
		InvDiffRmed[nr]));

	    // Q = (c_s kappa) / (Pi G Sigma)
	    // data[t_data::TOOMRE](n_radial, n_azimuthal) =
	    // data[t_data::SOUNDSPEED](n_radial,
	    // n_azimuthal)*calculate_omega_kepler(Rmed[n_radial])/(PI*G*data[t_data::DENSITY](n_radial,
	    // n_azimuthal));
		data[t_data::TOOMRE](nr, naz) = data[t_data::SOUNDSPEED](nr, naz) * kappa / (M_PI * constants::G * data[t_data::SIGMA](nr, naz));
	}
    }
}

void calculate_radial_luminosity(t_data &data, unsigned int timestep,
				 bool force_update)
{
    static int last_timestep_calculated = -1;

    if ((!force_update) && (last_timestep_calculated == (int)timestep)) {
	return;
    }

    last_timestep_calculated = timestep;

	const unsigned int Nr = data[t_data::LUMINOSITY_1D].get_size_radial();
	const unsigned int Nphi = data[t_data::QMINUS].get_size_azimuthal();

	#pragma omp parallel for
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	// L = int( int(sigma T^4 r ,phi) ,r);
	data[t_data::LUMINOSITY_1D](nr) = 0.0;

	for (unsigned int naz = 0; naz < Nphi; ++naz) {
		const double dr = (Rsup[nr] - Rinf[nr]);
		data[t_data::LUMINOSITY_1D](nr) +=
		data[t_data::QMINUS](nr, naz) * Rmed[nr] *
		dr * dphi;
	}
    }
}

void calculate_radial_dissipation(t_data &data, unsigned int timestep,
				  bool force_update)
{
    static int last_timestep_calculated = -1;

    if ((!force_update) && (last_timestep_calculated == (int)timestep)) {
	return;
    }

    last_timestep_calculated = timestep;

	const unsigned int Nr = data[t_data::DISSIPATION_1D].get_size_radial();
	const unsigned int Nphi = data[t_data::QPLUS].get_size_azimuthal();

	#pragma omp parallel for
	for (unsigned int nr = 0; nr < Nr; ++nr) {
	data[t_data::DISSIPATION_1D](nr) = 0.0;

	for (unsigned int naz = 0;  naz < Nphi; ++naz) {
		double dr = (Rsup[nr] - Rinf[nr]);

		data[t_data::DISSIPATION_1D](nr) +=
		data[t_data::QPLUS](nr, naz) * Rmed[nr] *
		dr * dphi;
	}
    }
}

void calculate_massflow(t_data &data, unsigned int timestep, bool force_update)
{
    (void)timestep;
    (void)force_update;

	const double denom = parameters::NINTERM * parameters::DT;

    // divide the data in massflow by the large timestep DT before writing out
    // to obtain the massflow from the mass difference
	 data[t_data::MASSFLOW] /= denom;
}


void compute_aspectratio(t_data &data, unsigned int timestep, bool force_update)
{
    static int last_timestep_calculated = -1;

    if ((!force_update) && (last_timestep_calculated == (int)timestep)) {
	return;
    }

	const unsigned int Nr = data[t_data::SCALE_HEIGHT].get_size_radial();
	const unsigned int Nphi = data[t_data::SCALE_HEIGHT].get_size_azimuthal();

    switch (parameters::ASPECTRATIO_MODE) {
    case 0: {
			#pragma omp parallel for collapse(2)
			for (unsigned int nr = 0; nr < Nr; ++nr) {
				for (unsigned int naz = 0; naz < Nphi; ++naz) {
					const double h = data[t_data::SCALE_HEIGHT](nr, naz) / Rb[nr];
					data[t_data::ASPECTRATIO](nr, naz) = h;
				}
			}

	break;
    }
    case 1: {

			static const unsigned int N_planets =
					data.get_planetary_system().get_number_of_planets();

			// setup planet data
			for (unsigned int k = 0; k < N_planets; k++) {
				const t_planet &planet = data.get_planetary_system().get_planet(k);
				g_mpl[k] = planet.get_rampup_mass(sim::PhysicalTime);
				g_xpl[k] = planet.get_x();
				g_ypl[k] = planet.get_y();
				g_rpl[k] = planet.get_planet_radial_extend();
			}

			// h = H/r
			// H = = c_s,iso / (GM/r^3) = c_s/sqrt(gamma) / / (GM/r^3)
			// for an Nbody system, H^-2 = sum_n (H_n)^-2
			// See G??nter & Kley 2003 Eq. 8, but beware of wrong extra square.
			// Better see Thun et al. 2017 Eq. 8 instead.
			#pragma omp parallel for collapse(2)
			for (unsigned int nr = 0; nr < Nr; ++nr) {
			for (unsigned int naz = 0; naz < Nphi; ++naz) {

					const int cell = get_cell_id(nr, naz);
					const double x = CellCenterX->Field[cell];
					const double y = CellCenterY->Field[cell];
					const double cs2 =
							std::pow(data[t_data::SOUNDSPEED](nr, naz), 2);

					double inv_h2 = 0.0; // inverse aspectratio squared

					for (unsigned int k = 0; k < N_planets; k++) {

						/// since the mass is distributed homogeniously distributed
						/// inside the cell, we assume that the planet is always at
						/// least cell_size / 2 plus planet radius away from the gas
						/// this is an rough estimate without explanation
						/// alternatively you can think about it yourself
						const double min_dist =
								0.5 * std::max(Rsup[nr] - Rinf[nr],
											   Rmed[nr] * dphi) +
								g_rpl[k];

						const double dx = x - g_xpl[k];
						const double dy = y - g_ypl[k];

						const double dist = std::max(
									std::sqrt(std::pow(dx, 2) + std::pow(dy, 2)), min_dist);

						// H^2 = (GM / dist^3 / Cs_iso^2)^-1
						if (parameters::Adiabatic || parameters::Polytropic) {
							const double gamma1 = pvte::get_gamma1(data, nr, naz);

								const double tmp_inv_h2 =
										constants::G * g_mpl[k] * gamma1 / (dist * cs2);
								inv_h2 += tmp_inv_h2;

						} else {

								const double tmp_inv_h2 =
										constants::G * g_mpl[k] / (dist * cs2);
								inv_h2 += tmp_inv_h2;
						}
					}

						const double h = std::sqrt(1.0 / inv_h2);
						data[t_data::ASPECTRATIO](nr, naz) = h;
				}
			}

	break;
    }
    case 2: {

		const Pair r_cm = data.get_planetary_system().get_center_of_mass();
		const double m_cm = data.get_planetary_system().get_mass();

		#pragma omp parallel for collapse(2)
		for (unsigned int nr = 0; nr < Nr; ++nr) {
		for (unsigned int n_az = 0; n_az < Nphi; ++n_az) {

			const int cell = get_cell_id(nr, n_az);
			const double x = CellCenterX->Field[cell];
			const double y = CellCenterY->Field[cell];
			const double cs = data[t_data::SOUNDSPEED](nr, n_az);

			// const double min_dist =
			//	0.5 * std::max(Rsup[n_rad] - Rinf[n_rad],
			//		   Rmed[n_rad] * dphi);

			const double dx = x - r_cm.x;
			const double dy = y - r_cm.y;

			// const double dist = std::max(
			//	std::sqrt(std::pow(dx, 2) + std::pow(dy, 2)), min_dist);
			const double dist = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));

			// h^2 = Cs_iso / vk = (Cs_iso^2 / (GM / dist))
			// H^2 = Cs_iso / Omegak = (Cs_iso^2 / (GM / dist^3))
			// H = h * dist
			if (parameters::Adiabatic || parameters::Polytropic) {
			/// Convert sound speed to isothermal sound speed cs,iso = cs /
			/// sqrt(gamma)
			const double gamma1 = pvte::get_gamma1(data, nr, n_az);
			const double h = cs * std::sqrt(dist / (constants::G * m_cm * gamma1));

			if(parameters::heating_star_enabled || parameters::self_gravity){
			data[t_data::ASPECTRATIO](nr, n_az) = h;
			}
			const double H = dist * h;
			data[t_data::SCALE_HEIGHT](nr, n_az) = H;

			} else { // locally isothermal
			const double h = cs * std::sqrt(dist / (constants::G * m_cm));
			if(parameters::heating_star_enabled || parameters::self_gravity){
			data[t_data::ASPECTRATIO](nr, n_az) = h;
			}
			const double H = dist * h;
			data[t_data::SCALE_HEIGHT](nr, n_az) = H;
			}
		}
	}
	break;
    }
    default: {
	#pragma omp parallel for collapse(2)
	for (unsigned int nr = 0; nr < Nr; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {
		data[t_data::ASPECTRATIO](nr, naz) =
			data[t_data::SCALE_HEIGHT](nr, naz) / Rmed[nr];
	    }
	}
    }
	}
}

void calculate_viscous_torque(t_data &data, unsigned int timestep,
			      bool force_update)
{
    (void)timestep;
    (void)force_update;

    if (!parameters::write_at_every_timestep) {
	const double denom = (double)parameters::NINTERM;
	data[t_data::VISCOUS_TORQUE] /= denom;
    }
}

void calculate_gravitational_torque(t_data &data, unsigned int timestep,
				    bool force_update)
{
    (void)timestep;
    (void)force_update;

    if (!parameters::write_at_every_timestep) {
	const double denom = (double)parameters::NINTERM;
	data[t_data::GRAVITATIONAL_TORQUE_NOT_INTEGRATED] /= denom;
    }
}

void calculate_advection_torque(t_data &data, unsigned int timestep,
				bool force_update)
{
    (void)timestep;
    (void)force_update;

    if (!parameters::write_at_every_timestep) {
	const double denom = (double)parameters::NINTERM;
	data[t_data::ADVECTION_TORQUE] /= denom;
    }
}



void CalculateMonitorQuantitiesAfterHydroStep(t_data &data,
						     int nTimeStep, double dt)
{
    if (data[t_data::ADVECTION_TORQUE].get_write()) {
	gas_torques::calculate_advection_torque(data, dt / parameters::DT);
    }
    if (data[t_data::VISCOUS_TORQUE].get_write()) {
	gas_torques::calculate_viscous_torque(data, dt / parameters::DT);
    }
    if (data[t_data::GRAVITATIONAL_TORQUE_NOT_INTEGRATED].get_write()) {
	gas_torques::calculate_gravitational_torque(data, dt / parameters::DT);
    }

    if (data[t_data::ALPHA_GRAV_MEAN].get_write()) {
	calculate_alpha_grav_mean_sumup(data, nTimeStep, dt / parameters::DT);
    }
    if (data[t_data::ALPHA_REYNOLDS_MEAN].get_write()) {
	calculate_alpha_reynolds_mean_sumup(data, nTimeStep,
							dt / parameters::DT);
    }
}


} // namespace quantities
