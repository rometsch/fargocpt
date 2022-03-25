#include "quantities.h"

#include "Theo.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "parameters.h"
#include "stress.h"
#include "util.h"
#include <math.h>
#include <mpi.h>
#include <vector>

extern boolean Corotating;
extern double M0;

namespace quantities
{


/**
	Calculates total gas mass.
*/
double gas_total_mass(t_data &data, const double quantitiy_radius)
{
	double local_mass = 0.0;
	double global_mass = 0.0;

	// calculate mass of this process' cells
	for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::DENSITY].get_size_azimuthal();
		 ++n_azimuthal) {
		if (Rmed[n_radial] <= quantitiy_radius) {
		local_mass += Surf[n_radial] *
				  data[t_data::DENSITY](n_radial, n_azimuthal);
		}
	}
	}

	MPI_Allreduce(&local_mass, &global_mass, 1, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);

	return global_mass;
}

double gas_aspect_ratio(t_data &data, const double quantitiy_radius)
{
    double local_mass = 0.0;
	const double gas_total_mass = quantities::gas_total_mass(data, quantitiy_radius);
    double aspect_ratio = 0.0;
    double local_aspect_ratio = 0.0;

    // Loop thru all cells excluding GHOSTCELLS & CPUOVERLAP cells (otherwise
    // they would be included twice!)
	for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::DENSITY].get_size_azimuthal();
	     ++n_azimuthal) {
	    // eccentricity and semi major axis weighted with cellmass
		if (Rmed[n_radial] <= quantitiy_radius) {
		local_mass = data[t_data::DENSITY](n_radial, n_azimuthal) *
			     Surf[n_radial];
		local_aspect_ratio +=
			data[t_data::ASPECTRATIO](n_radial, n_azimuthal) *
		    local_mass;
	    }
	}
    }

    // synchronize threads
    MPI_Allreduce(&local_aspect_ratio, &aspect_ratio, 1, MPI_DOUBLE, MPI_SUM,
		  MPI_COMM_WORLD);

    aspect_ratio /= gas_total_mass;

    return aspect_ratio;
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
    const unsigned int send_size = local_array_end - local_array_start;

    std::vector<double> local_mass(send_size);

    for (unsigned int n_radial = local_array_start; n_radial < local_array_end;
	 ++n_radial) {
	local_mass[n_radial - local_array_start] = 0.0;
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::DENSITY].get_size_azimuthal();
	     ++n_azimuthal) {
	    local_mass[n_radial - local_array_start] +=
		Surf[n_radial] * data[t_data::DENSITY](n_radial, n_azimuthal);
	}
    }
    double radius = 0.0;
    double current_mass = 0.0;

    MPI_Gatherv(&local_mass[0], send_size, MPI_DOUBLE, GLOBAL_bufarray,
		RootNradialLocalSizes, RootNradialDisplacements, MPI_DOUBLE, 0,
		MPI_COMM_WORLD);

    if (CPU_Master) {
	int j = -1;
	for (int rank = 0; rank < CPU_Number; ++rank) {
	    int id = RootRanksOrdered[rank];
	    for (int i = RootIMIN[id]; i <= RootIMAX[id]; ++i) {
		++j;
		current_mass += GLOBAL_bufarray[i];
		if (current_mass > 0.99 * total_mass) {
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

	for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::DENSITY].get_size_azimuthal();
	     ++n_azimuthal) {
		if (Rmed[n_radial] <= quantitiy_radius) {
		local_angular_momentum +=
		    Surf[n_radial] * 0.5 *
		    (data[t_data::DENSITY](n_radial, n_azimuthal) +
		     data[t_data::DENSITY](
			 n_radial,
			 n_azimuthal == 0
			     ? data[t_data::DENSITY].get_max_azimuthal()
			     : n_azimuthal - 1)) *
		    Rmed[n_radial] *
		    (data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) +
		     OmegaFrame * Rmed[n_radial]);
	    }
	    // local_angular_momentum +=
	    // Surf[n_radial]*data[t_data::DENSITY](n_radial,n_azimuthal)*Rmed[n_radial]*(0.5*(data[t_data::V_AZIMUTHAL](n_radial,n_azimuthal)+data[t_data::V_AZIMUTHAL](n_radial,n_azimuthal
	    // == data[t_data::V_AZIMUTHAL].get_max_azimuthal() ? 0 :
	    // n_azimuthal+1))+OmegaFrame*Rmed[n_radial]);
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

	for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size;
	 ++n_radial) {
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

	for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size;
	 ++n_radial) {
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

	for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size;
	 ++n_radial) {
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

    double v_radial_center, v_azimuthal_center;

	for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::DENSITY].get_size_azimuthal();
	     ++n_azimuthal) {
	    // centered-in-cell radial velocity
		if (Rmed[n_radial] <= quantitiy_radius) {
		v_radial_center =
		    (Rmed[n_radial] - Rinf[n_radial]) *
			data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) +
		    (Rsup[n_radial] - Rmed[n_radial]) *
			data[t_data::V_RADIAL](n_radial, n_azimuthal);
		v_radial_center /= (Rsup[n_radial] - Rinf[n_radial]);

		// centered-in-cell azimuthal velocity
		v_azimuthal_center =
		    0.5 *
			(data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) +
			 data[t_data::V_AZIMUTHAL](
			     n_radial, n_azimuthal == data[t_data::V_AZIMUTHAL]
							  .get_max_azimuthal()
					   ? 0
					   : n_azimuthal + 1)) +
		    Rmed[n_radial] * OmegaFrame;

		local_kinematic_energy +=
		    0.5 * Surf[n_radial] *
		    data[t_data::DENSITY](n_radial, n_azimuthal) *
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

    double v_radial_center;

	for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::DENSITY].get_size_azimuthal();
	     ++n_azimuthal) {
		if (Rmed[n_radial] <= quantitiy_radius) {
		// centered-in-cell radial velocity
		v_radial_center =
		    (Rmed[n_radial] - Rinf[n_radial]) *
			data[t_data::V_RADIAL](n_radial + 1, n_azimuthal) +
		    (Rsup[n_radial] - Rmed[n_radial]) *
			data[t_data::V_RADIAL](n_radial, n_azimuthal);
		v_radial_center /= (Rsup[n_radial] - Rinf[n_radial]);

		local_kinematic_energy +=
		    0.5 * Surf[n_radial] *
		    data[t_data::DENSITY](n_radial, n_azimuthal) *
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
double gas_azimuthal_kinematic_energy(t_data &data, const double quantitiy_radius)
{
    double local_kinematic_energy = 0.0;
    double global_kinematic_energy = 0.0;

    double v_azimuthal_center;

	for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::DENSITY].get_size_azimuthal();
	     ++n_azimuthal) {
		if (Rmed[n_radial] <= quantitiy_radius) {
		// centered-in-cell azimuthal velocity
		v_azimuthal_center =
		    0.5 *
			(data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) +
			 data[t_data::V_AZIMUTHAL](
			     n_radial, n_azimuthal == data[t_data::V_AZIMUTHAL]
							  .get_max_azimuthal()
					   ? 0
					   : n_azimuthal + 1)) +
		    Rmed[n_radial] * OmegaFrame;

		local_kinematic_energy +=
		    0.5 * Surf[n_radial] *
		    data[t_data::DENSITY](n_radial, n_azimuthal) *
		    std::pow(v_azimuthal_center, 2);
	    }
	}
    }

    MPI_Reduce(&local_kinematic_energy, &global_kinematic_energy, 1, MPI_DOUBLE,
	       MPI_SUM, 0, MPI_COMM_WORLD);

    return global_kinematic_energy;
}

/**
	Calculates gas gravitational energy
*/
double gas_gravitational_energy(t_data &data, const double quantitiy_radius)
{
    double local_gravitational_energy = 0.0;
    double global_gravitational_energy = 0.0;

	for (unsigned int n_radial = radial_first_active;
	 n_radial < radial_active_size;
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::DENSITY].get_size_azimuthal();
	     ++n_azimuthal) {
		if (Rmed[n_radial] <= quantitiy_radius) {
		local_gravitational_energy +=
		    -Surf[n_radial] *
		    data[t_data::DENSITY](n_radial, n_azimuthal) *
		    data[t_data::POTENTIAL](n_radial, n_azimuthal);
	    }
	}
    }

    MPI_Reduce(&local_gravitational_energy, &global_gravitational_energy, 1,
	       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    return global_gravitational_energy;
}

/**
	Calculates disk (grid) eccentricity
*/
void calculate_disk_quantities(t_data &data, unsigned int timestep,
			       bool force_update)
{
    static int last_timestep_calculated = -1;
    double angle, r_x, r_y, j;
    double v_xmed, v_ymed;
    double e_x, e_y;
    double total_mass = 0.0;

    if (!force_update) {
	if (last_timestep_calculated == (int)timestep) {
	    return;
	} else {
	    last_timestep_calculated = timestep;
	}
    }
    // calculations outside the loop for speedup
    double sinFrameAngle = std::sin(FrameAngle);
    double cosFrameAngle = std::cos(FrameAngle);
    for (unsigned int n_radial = 0;
	 n_radial < data[t_data::DENSITY].get_size_radial(); ++n_radial) {
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::DENSITY].get_size_azimuthal();
	     ++n_azimuthal) {
	    total_mass =
		hydro_center_mass +
		data[t_data::DENSITY](n_radial, n_azimuthal) * Surf[n_radial];

	    // location of the cell
	    angle = (double)n_azimuthal /
		    (double)data[t_data::V_RADIAL].get_size_azimuthal() * 2.0 *
		    M_PI;
	    r_x = Rmed[n_radial] * std::cos(angle);
	    r_y = Rmed[n_radial] * std::sin(angle);

	    // averaged velocities
	    v_xmed =
		std::cos(angle) * 0.5 *
		    (data[t_data::V_RADIAL](n_radial, n_azimuthal) +
		     data[t_data::V_RADIAL](n_radial + 1, n_azimuthal)) -
		std::sin(angle) *
		    (0.5 *
			 (data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) +
			  data[t_data::V_AZIMUTHAL](
			      n_radial, n_azimuthal == data[t_data::V_AZIMUTHAL]
							   .get_max_azimuthal()
					    ? 0
					    : n_azimuthal + 1)) +
		     OmegaFrame * Rmed[n_radial]);
	    v_ymed =
		std::sin(angle) * 0.5 *
		    (data[t_data::V_RADIAL](n_radial, n_azimuthal) +
		     data[t_data::V_RADIAL](n_radial + 1, n_azimuthal)) +
		std::cos(angle) *
		    (0.5 *
			 (data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) +
			  data[t_data::V_AZIMUTHAL](
			      n_radial, n_azimuthal == data[t_data::V_AZIMUTHAL]
							   .get_max_azimuthal()
					    ? 0
					    : n_azimuthal + 1)) +
		     OmegaFrame * Rmed[n_radial]);

	    // specific angular momentum for each cell j = j*e_z
	    j = r_x * v_ymed - r_y * v_xmed;
	    // Runge-Lenz vector Ax = x*vy*vy-y*vx*vy-G*m*x/d;
	    e_x =
		j * v_ymed / (constants::G * total_mass) - r_x / Rmed[n_radial];
	    e_y = -1.0 * j * v_xmed / (constants::G * total_mass) -
		  r_y / Rmed[n_radial];

	    data[t_data::ECCENTRICITY](n_radial, n_azimuthal) =
		std::sqrt(std::pow(e_x, 2) + std::pow(e_y, 2));

	    if (FrameAngle != 0.0) {
		// periastron grid is rotated to non-rotating coordinate system
		// to prevent phase jumps of atan2 in later transformations like
		// you would have had if you back-transform the output
		// periastron values
		data[t_data::PERIASTRON](n_radial, n_azimuthal) =
		    std::atan2(e_y * cosFrameAngle + e_x * sinFrameAngle,
			       e_x * cosFrameAngle - e_y * sinFrameAngle);
	    } else {
		data[t_data::PERIASTRON](n_radial, n_azimuthal) =
		    std::atan2(e_y, e_x);
	    }
	}
    }
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

    for (unsigned int n_radial = 0;
	 n_radial < data[t_data::ALPHA_GRAV].get_size_radial(); ++n_radial) {
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::ALPHA_GRAV].get_size_azimuthal();
	     ++n_azimuthal) {
	    /*data[t_data::ALPHA_GRAV](n_radial, n_azimuthal) = -2.0/3.0 *
	     * (data[t_data::T_GRAVITATIONAL](n_radial,n_azimuthal)+data[t_data::T_REYNOLDS](n_radial,
	     * n_azimuthal))/(data[t_data::DENSITY](n_radial,n_azimuthal)*pow2(data[t_data::SOUNDSPEED](n_radial,
	     * n_azimuthal)));*/
	    data[t_data::ALPHA_GRAV](n_radial, n_azimuthal) =
		2.0 / 3.0 *
		data[t_data::T_GRAVITATIONAL](n_radial, n_azimuthal) /
		(data[t_data::DENSITY](n_radial, n_azimuthal) *
		 std::pow(data[t_data::SOUNDSPEED](n_radial, n_azimuthal), 2));
	}
    }
}

void calculate_alpha_grav_mean_sumup(t_data &data, unsigned int timestep,
				     double dt)
{
    calculate_alpha_grav(data, timestep, true);

    for (unsigned int n_radial = 0;
	 n_radial < data[t_data::ALPHA_GRAV_MEAN].get_size_radial();
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::ALPHA_GRAV_MEAN].get_size_azimuthal();
	     ++n_azimuthal) {
	    data[t_data::ALPHA_GRAV_MEAN](n_radial, n_azimuthal) +=
		data[t_data::ALPHA_GRAV](n_radial, n_azimuthal) * dt;
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

    for (unsigned int n_radial = 0;
	 n_radial < data[t_data::ALPHA_REYNOLDS].get_size_radial();
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::ALPHA_REYNOLDS].get_size_azimuthal();
	     ++n_azimuthal) {
	    data[t_data::ALPHA_REYNOLDS](n_radial, n_azimuthal) =
		2.0 / 3.0 * data[t_data::T_REYNOLDS](n_radial, n_azimuthal) /
		(data[t_data::DENSITY](n_radial, n_azimuthal) *
		 std::pow(data[t_data::SOUNDSPEED](n_radial, n_azimuthal), 2));
	}
    }
}

void calculate_alpha_reynolds_mean_sumup(t_data &data, unsigned int timestep,
					 double dt)
{
    calculate_alpha_reynolds(data, timestep, true);

    for (unsigned int n_radial = 0;
	 n_radial < data[t_data::ALPHA_REYNOLDS_MEAN].get_size_radial();
	 ++n_radial) {
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal <
		 data[t_data::ALPHA_REYNOLDS_MEAN].get_size_azimuthal();
	     ++n_azimuthal) {
	    data[t_data::ALPHA_REYNOLDS_MEAN](n_radial, n_azimuthal) +=
		data[t_data::ALPHA_REYNOLDS](n_radial, n_azimuthal) * dt;
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

    for (unsigned int n_radial = 1;
	 n_radial < data[t_data::TOOMRE].get_size_radial(); ++n_radial) {
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::TOOMRE].get_size_azimuthal();
	     ++n_azimuthal) {
	    // kappa^2 = 1/r^3 d((r^2 Omega)^2)/dr = 1/r^3 d((r*v_phi)^2)/dr
	    kappa = std::sqrt(std::fabs(
		std::pow(InvRmed[n_radial], 3) *
		(std::pow(data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal) *
			      Rmed[n_radial],
			  2) -
		 std::pow(data[t_data::V_AZIMUTHAL](n_radial - 1, n_azimuthal) *
			      Rmed[n_radial - 1],
			  2)) *
		InvDiffRmed[n_radial]));

	    // Q = (c_s kappa) / (Pi G Sigma)
	    // data[t_data::TOOMRE](n_radial, n_azimuthal) =
	    // data[t_data::SOUNDSPEED](n_radial,
	    // n_azimuthal)*calculate_omega_kepler(Rmed[n_radial])/(PI*G*data[t_data::DENSITY](n_radial,
	    // n_azimuthal));
	    data[t_data::TOOMRE](n_radial, n_azimuthal) =
		data[t_data::SOUNDSPEED](n_radial, n_azimuthal) * kappa /
		(M_PI * constants::G *
		 data[t_data::DENSITY](n_radial, n_azimuthal));
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

    for (unsigned int n_radial = 0;
	 n_radial < data[t_data::LUMINOSITY_1D].get_size_radial(); ++n_radial) {
	// L = int( int(sigma T^4 r ,phi) ,r);
	data[t_data::LUMINOSITY_1D](n_radial) = 0.0;

	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::QMINUS].get_size_azimuthal();
	     ++n_azimuthal) {
	    double dr = (Rsup[n_radial] - Rinf[n_radial]);
	    data[t_data::LUMINOSITY_1D](n_radial) +=
		data[t_data::QMINUS](n_radial, n_azimuthal) * Rmed[n_radial] *
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

    for (unsigned int n_radial = 0;
	 n_radial < data[t_data::DISSIPATION_1D].get_size_radial();
	 ++n_radial) {
	data[t_data::DISSIPATION_1D](n_radial) = 0.0;

	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::QPLUS].get_size_azimuthal();
	     ++n_azimuthal) {
	    double dr = (Rsup[n_radial] - Rinf[n_radial]);

	    data[t_data::DISSIPATION_1D](n_radial) +=
		data[t_data::QPLUS](n_radial, n_azimuthal) * Rmed[n_radial] *
		dr * dphi;
	}
    }
}

void calculate_massflow(t_data &data, unsigned int timestep, bool force_update)
{
    (void)timestep;
    (void)force_update;

    double denom;
    denom = NINTERM * DT;

    // divide the data in massflow by the large timestep DT before writing out
    // to obtain the massflow from the mass difference
    for (unsigned int nRadial = 0;
	 nRadial < data[t_data::MASSFLOW].get_size_radial(); ++nRadial) {
	for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::MASSFLOW].get_size_azimuthal();
	     ++n_azimuthal) {
	    data[t_data::MASSFLOW](nRadial, n_azimuthal) *= 1. / denom;
	}
    }
}

void compute_aspectratio(t_data &data, unsigned int timestep, bool force_update)
{
	static int last_timestep_calculated = -1;

	if ((!force_update) && (last_timestep_calculated == (int)timestep)) {
	return;
	}

	switch(ASPECTRATIO_MODE){
		case 0:
		{
			for (unsigned int nRad = 0;
			 nRad < data[t_data::ASPECTRATIO].get_size_radial(); ++nRad) {
			for (unsigned int nAz = 0;
				 nAz < data[t_data::ASPECTRATIO].get_size_azimuthal();
				 ++nAz) {
			data[t_data::ASPECTRATIO](nRad, nAz) = data[t_data::SCALE_HEIGHT](nRad, nAz) / Rmed[nRad];
			}
			}

			break;
		}
		case 1:
		{
			static const unsigned int N_planets =
			data.get_planetary_system().get_number_of_planets();
			static std::vector<double> xpl(N_planets);
			static std::vector<double> ypl(N_planets);
			static std::vector<double> mpl(N_planets);
			static std::vector<double> rpl(N_planets);

			// setup planet data
			for (unsigned int k = 0; k < N_planets; k++) {
			t_planet &planet = data.get_planetary_system().get_planet(k);
			mpl[k] = planet.get_rampup_mass();
			xpl[k] = planet.get_x();
			ypl[k] = planet.get_y();
			rpl[k] = planet.get_planet_radial_extend();
			}

			const Pair r_cm = data.get_planetary_system().get_center_of_mass();
			const double m_cm = data.get_planetary_system().get_mass();

			for (unsigned int nRad = 0;
			 nRad < data[t_data::ASPECTRATIO].get_size_radial(); ++nRad) {
			for (unsigned int nAz = 0;
				 nAz < data[t_data::ASPECTRATIO].get_size_azimuthal();
				 ++nAz) {

				const int cell = get_cell_id(nRad, nAz);
				const double x = CellCenterX->Field[cell];
				const double y = CellCenterY->Field[cell];

				// cell_r is the distance to the closest body used for computing the sound speed
				// the cell belongs to a body, if it is inside its roche radius.
				// if no close body is found, the center of mass is used instead
				double cell_r = 0.0;
				double roche_radius;

				for (unsigned int k = 0; k < N_planets; k++) {

					// primary object uses next object to compute the roche radius
					// while all other objects use the primary object.
					if(k == 0){
						const double partner_dist = std::sqrt(std::pow(xpl[k] - xpl[1], 2) + std::pow(ypl[k] - ypl[1], 2));
						const double mass_q = mpl[k]/m_cm / (1.0 - mpl[k]/m_cm);
						roche_radius = eggleton_1983(mass_q , partner_dist);
					} else {
						const double partner_dist = std::sqrt(std::pow(xpl[k] - xpl[0], 2) + std::pow(ypl[k] - ypl[0], 2));
						const double mass_q = mpl[k]/m_cm / (1.0 - mpl[k]/m_cm);
						roche_radius = eggleton_1983(mass_q , partner_dist);
					}

				/// since the mass is distributed homogeniously distributed
				/// inside the cell, we assume that the planet is always at
				/// least cell_size / 2 plus planet radius away from the gas
				/// this is an rough estimate without explanation
				/// alternatively you can think about it yourself
				const double min_dist =
					0.5 * std::max(Rsup[nRad] - Rinf[nRad],
						   Rmed[nRad] * dphi) +
					rpl[k];

				const double dx = x - xpl[k];
				const double dy = y - ypl[k];

				const double dist = std::max(
					std::sqrt(std::pow(dx, 2) + std::pow(dy, 2)), min_dist);

				if(dist < roche_radius){
					cell_r = dist;
				}
				}

				if(cell_r == 0.0){
					cell_r = std::sqrt(std::pow(x - r_cm.x, 2) + std::pow(y - r_cm.y, 2));
				}

				data[t_data::ASPECTRATIO](nRad, nAz) = data[t_data::SCALE_HEIGHT](nRad, nAz) / cell_r;
			}
			}
			break;
		}
		case 2:
		{
			const Pair r_cm = data.get_planetary_system().get_center_of_mass();

			for (unsigned int nRad = 0;
			 nRad < data[t_data::ASPECTRATIO].get_size_radial(); ++nRad) {
			for (unsigned int nAz = 0;
				 nAz < data[t_data::ASPECTRATIO].get_size_azimuthal();
				 ++nAz) {

				const int cell = get_cell_id(nRad, nAz);
				const double x = CellCenterX->Field[cell];
				const double y = CellCenterY->Field[cell];


				/// since the mass is distributed homogeniously distributed
				/// inside the cell, we assume that the planet is always at
				/// least cell_size / 2 plus planet radius away from the gas
				/// this is an rough estimate without explanation
				/// alternatively you can think about it yourself
				//const double min_dist =
				//	0.5 * std::max(Rsup[nRad] - Rinf[nRad],
				//		   Rmed[nRad] * dphi);

				const double dx = x - r_cm.x;
				const double dy = y - r_cm.y;

				//const double dist = std::max(
				//	std::sqrt(std::pow(dx, 2) + std::pow(dy, 2)), min_dist);

				const double dist = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));
				data[t_data::ASPECTRATIO](nRad, nAz) = data[t_data::SCALE_HEIGHT](nRad, nAz) / dist;
			}
			}
			break;
		}
		default:
		{
				for (unsigned int nRad = 0;
				 nRad < data[t_data::ASPECTRATIO].get_size_radial(); ++nRad) {
				for (unsigned int nAz = 0;
					 nAz < data[t_data::ASPECTRATIO].get_size_azimuthal();
					 ++nAz) {
				data[t_data::ASPECTRATIO](nRad, nAz) = data[t_data::SCALE_HEIGHT](nRad, nAz) / Rmed[nRad];
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

    double denom;

    if (!parameters::write_at_every_timestep) {
	denom = (double)NINTERM;
	// divide the data in massflow by the large timestep DT before writing
	// out to obtain the massflow from the mass difference
	for (unsigned int nRadial = 0;
	     nRadial < data[t_data::VISCOUS_TORQUE].get_size_radial();
	     ++nRadial) {
	    for (unsigned int nAzimuthal = 0;
		 nAzimuthal < data[t_data::VISCOUS_TORQUE].get_size_azimuthal();
		 ++nAzimuthal) {
		data[t_data::VISCOUS_TORQUE](nRadial, nAzimuthal) *= 1. / denom;
	    }
	}
    }
}

void calculate_gravitational_torque(t_data &data, unsigned int timestep,
				    bool force_update)
{
    (void)timestep;
    (void)force_update;

    double denom;

    if (!parameters::write_at_every_timestep) {
	denom = (double)NINTERM;
	// divide the data in massflow by the large timestep DT before writing
	// out to obtain the massflow from the mass difference
	for (unsigned int nRadial = 0;
	     nRadial < data[t_data::GRAVITATIONAL_TORQUE_NOT_INTEGRATED]
			   .get_size_radial();
	     ++nRadial) {
	    for (unsigned int nAzimuthal = 0;
		 nAzimuthal < data[t_data::GRAVITATIONAL_TORQUE_NOT_INTEGRATED]
				  .get_size_azimuthal();
		 ++nAzimuthal) {
		data[t_data::GRAVITATIONAL_TORQUE_NOT_INTEGRATED](
		    nRadial, nAzimuthal) *= 1. / denom;
	    }
	}
    }
}

void calculate_advection_torque(t_data &data, unsigned int timestep,
				bool force_update)
{
    (void)timestep;
    (void)force_update;

    double denom;

    if (!parameters::write_at_every_timestep) {
	denom = (double)NINTERM;
	// divide the data in massflow by the large timestep DT before writing
	// out to obtain the massflow from the mass difference
	for (unsigned int nRadial = 0;
	     nRadial < data[t_data::ADVECTION_TORQUE].get_size_radial();
	     ++nRadial) {
	    for (unsigned int nAzimuthal = 0;
		 nAzimuthal <
		 data[t_data::ADVECTION_TORQUE].get_size_azimuthal();
		 ++nAzimuthal) {
		data[t_data::ADVECTION_TORQUE](nRadial, nAzimuthal) *=
		    1. / denom;
	    }
	}
    }
}

} // namespace quantities
