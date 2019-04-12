#include "quantities.h"

#include <mpi.h>
#include <math.h>
#include "global.h"
#include "constants.h"
#include "util.h"
#include "stress.h"
#include "parameters.h"
#include "logging.h"
#include "Theo.h"

extern boolean Corotating;
extern double M0;

namespace quantities {

/**
	Calculates total gas mass.
*/
double gas_total_mass(t_data &data)
{
	double local_mass = 0.0;
	double global_mass = 0.0;

	// calculate mass of this process' cells
	for (unsigned int n_radial = (CPU_Rank == 0) ? GHOSTCELLS_B : CPUOVERLAP; n_radial <= data[t_data::DENSITY].get_max_radial() - ( CPU_Rank == CPU_Highest ? GHOSTCELLS_B : CPUOVERLAP ); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
			local_mass += Surf[n_radial]*data[t_data::DENSITY](n_radial, n_azimuthal);
		}
	}

	MPI_Allreduce(&local_mass, &global_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return global_mass;
}

/**
	Calculates angular gas momentum.
*/
double gas_angular_momentum(t_data &data)
{
	double local_angular_momentum = 0.0;
	double global_angular_momentum = 0.0;

	for (unsigned int n_radial = (CPU_Rank == 0) ? (GHOSTCELLS_B) : CPUOVERLAP; n_radial <= data[t_data::DENSITY].get_max_radial() - ( CPU_Rank == CPU_Highest ? (GHOSTCELLS_B) : CPUOVERLAP ); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
			local_angular_momentum += Surf[n_radial]*0.5*(data[t_data::DENSITY](n_radial,n_azimuthal)+data[t_data::DENSITY](n_radial,n_azimuthal == 0 ? data[t_data::DENSITY].get_max_azimuthal() : n_azimuthal-1))*Rmed[n_radial]*(data[t_data::V_AZIMUTHAL](n_radial,n_azimuthal)+OmegaFrame*Rmed[n_radial]);
			//local_angular_momentum += Surf[n_radial]*data[t_data::DENSITY](n_radial,n_azimuthal)*Rmed[n_radial]*(0.5*(data[t_data::V_AZIMUTHAL](n_radial,n_azimuthal)+data[t_data::V_AZIMUTHAL](n_radial,n_azimuthal == data[t_data::V_AZIMUTHAL].get_max_azimuthal() ? 0 : n_azimuthal+1))+OmegaFrame*Rmed[n_radial]);
		}
	}

	MPI_Allreduce(&local_angular_momentum, &global_angular_momentum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return global_angular_momentum;
}

/**
	Calculates gas internal energy
*/
double gas_internal_energy(t_data &data)
{
	double local_internal_energy = 0.0;
	double global_internal_energy = 0.0;

	for (unsigned int n_radial = (CPU_Rank == 0) ? GHOSTCELLS_B : CPUOVERLAP; n_radial < data[t_data::ENERGY].get_max_radial() - ( CPU_Rank == CPU_Highest ? GHOSTCELLS_B : CPUOVERLAP ); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ENERGY].get_max_azimuthal(); ++n_azimuthal) {
			local_internal_energy += Surf[n_radial]*data[t_data::ENERGY](n_radial,n_azimuthal);
		}
	}

	MPI_Allreduce(&local_internal_energy, &global_internal_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return global_internal_energy;
}

/**
	Calculates gas kinematic energy
*/
double gas_kinematic_energy(t_data &data)
{
	double local_kinematic_energy = 0.0;
	double global_kinematic_energy = 0.0;

	double v_radial_center, v_azimuthal_center;

	for (unsigned int n_radial = (CPU_Rank == 0) ? GHOSTCELLS_B : CPUOVERLAP; n_radial <= data[t_data::DENSITY].get_max_radial() - ( CPU_Rank == CPU_Highest ? GHOSTCELLS_B : CPUOVERLAP ); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
			// centered-in-cell radial velocity
			v_radial_center = (Rmed[n_radial]-Rinf[n_radial])*data[t_data::V_RADIAL](n_radial+1,n_azimuthal) + (Rsup[n_radial]-Rmed[n_radial])*data[t_data::V_RADIAL](n_radial, n_azimuthal);
			v_radial_center /= (Rsup[n_radial]-Rinf[n_radial]);

			// centered-in-cell azimuthal velocity
			v_azimuthal_center = 0.5*(data[t_data::V_AZIMUTHAL](n_radial,n_azimuthal)+data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal == data[t_data::V_AZIMUTHAL].get_max_azimuthal() ? 0 : n_azimuthal+1)) + Rmed[n_radial]*OmegaFrame;

			local_kinematic_energy += 0.5*Surf[n_radial]*data[t_data::DENSITY](n_radial,n_azimuthal)*(pow2(v_radial_center) + pow2(v_azimuthal_center));
		}
	}

	MPI_Allreduce(&local_kinematic_energy, &global_kinematic_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return global_kinematic_energy;
}

/**
	Calculates gas kinematic energy
*/
double gas_radial_kinematic_energy(t_data &data)
{
	double local_kinematic_energy = 0.0;
	double global_kinematic_energy = 0.0;

	double v_radial_center;

	for (unsigned int n_radial = (CPU_Rank == 0) ? GHOSTCELLS_B : CPUOVERLAP; n_radial <= data[t_data::DENSITY].get_max_radial() - ( CPU_Rank == CPU_Highest ? GHOSTCELLS_B : CPUOVERLAP ); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
			// centered-in-cell radial velocity
			v_radial_center = (Rmed[n_radial]-Rinf[n_radial])*data[t_data::V_RADIAL](n_radial+1,n_azimuthal) + (Rsup[n_radial]-Rmed[n_radial])*data[t_data::V_RADIAL](n_radial, n_azimuthal);
			v_radial_center /= (Rsup[n_radial]-Rinf[n_radial]);

			local_kinematic_energy += 0.5*Surf[n_radial]*data[t_data::DENSITY](n_radial,n_azimuthal)*pow2(v_radial_center);
		}
	}

	MPI_Allreduce(&local_kinematic_energy, &global_kinematic_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return global_kinematic_energy;
}

/**
	Calculates gas kinematic energy
*/
double gas_azimuthal_kinematic_energy(t_data &data)
{
	double local_kinematic_energy = 0.0;
	double global_kinematic_energy = 0.0;

	double v_azimuthal_center;

	for (unsigned int n_radial = (CPU_Rank == 0) ? GHOSTCELLS_B : CPUOVERLAP; n_radial <= data[t_data::DENSITY].get_max_radial() - ( CPU_Rank == CPU_Highest ? GHOSTCELLS_B : CPUOVERLAP ); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
			// centered-in-cell azimuthal velocity
			v_azimuthal_center = 0.5*(data[t_data::V_AZIMUTHAL](n_radial,n_azimuthal)+data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal == data[t_data::V_AZIMUTHAL].get_max_azimuthal() ? 0 : n_azimuthal+1)) + Rmed[n_radial]*OmegaFrame;

			local_kinematic_energy += 0.5*Surf[n_radial]*data[t_data::DENSITY](n_radial,n_azimuthal)*pow2(v_azimuthal_center);
		}
	}

	MPI_Allreduce(&local_kinematic_energy, &global_kinematic_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return global_kinematic_energy;
}

/**
	Calculates gas gravitational energy
*/
double gas_gravitational_energy(t_data &data)
{
	double local_gravitational_energy = 0.0;
	double global_gravitational_energy = 0.0;

	for (unsigned int n_radial = (CPU_Rank == 0) ? GHOSTCELLS_B : CPUOVERLAP; n_radial <= data[t_data::DENSITY].get_max_radial() - ( CPU_Rank == CPU_Highest ? GHOSTCELLS_B : CPUOVERLAP ); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
			local_gravitational_energy += - Surf[n_radial]*data[t_data::DENSITY](n_radial,n_azimuthal)*data[t_data::POTENTIAL](n_radial,n_azimuthal);
		}
	}

	MPI_Allreduce(&local_gravitational_energy, &global_gravitational_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	return global_gravitational_energy;
}

/**
	Calculates total gas energy.
*/
double gas_total_energy(t_data &data)
{
	return gas_kinematic_energy(data)+gas_internal_energy(data)+gas_gravitational_energy(data);
}

/**
	Calculates disk (grid) eccentricity
*/
void calculate_disk_quantities(t_data &data, unsigned int timestep, bool force_update)
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
	double sinFrameAngle = sin(FrameAngle);
	double cosFrameAngle = cos(FrameAngle);
	for (unsigned int n_radial = 0; n_radial <= data[t_data::DENSITY].get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::DENSITY].get_max_azimuthal(); ++n_azimuthal) {
			total_mass = M + data[t_data::DENSITY](n_radial, n_azimuthal)*Surf[n_radial];

			// location of the cell
			angle = (double)n_azimuthal/(double)data[t_data::V_RADIAL].get_size_azimuthal()*2.0*PI;
			r_x = Rmed[n_radial]*cos(angle);
			r_y = Rmed[n_radial]*sin(angle);

			// averaged velocities
			v_xmed = cos(angle) * 0.5 * (data[t_data::V_RADIAL](n_radial, n_azimuthal)+data[t_data::V_RADIAL](n_radial+1, n_azimuthal))
				- sin(angle) *( 0.5 * (data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal)+data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal == data[t_data::V_AZIMUTHAL].get_max_azimuthal() ? 0 : n_azimuthal+1))+OmegaFrame*Rmed[n_radial]);
			v_ymed = sin(angle) * 0.5 * (data[t_data::V_RADIAL](n_radial, n_azimuthal)+data[t_data::V_RADIAL](n_radial+1, n_azimuthal))
				+ cos(angle) *( 0.5 * (data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal)+data[t_data::V_AZIMUTHAL](n_radial, n_azimuthal == data[t_data::V_AZIMUTHAL].get_max_azimuthal() ? 0 : n_azimuthal+1))+OmegaFrame*Rmed[n_radial]);

			// specific angular momentum for each cell j = j*e_z
			j = r_x * v_ymed - r_y * v_xmed;
			// Runge-Lenz vector Ax = x*vy*vy-y*vx*vy-G*m*x/d;
			e_x =        j * v_ymed / (constants::G * total_mass) - r_x/Rmed[n_radial];
			e_y = -1.0 * j * v_xmed / (constants::G * total_mass) - r_y/Rmed[n_radial];

			data[t_data::ECCENTRICITY](n_radial, n_azimuthal) = sqrt( pow2(e_x) + pow2(e_y) );

			if (FrameAngle != 0.0) {
			//periastron grid is rotated to non-rotating coordinate system to prevent phase jumps of atan2 in later transformations like you would have had if you back-transform the output periastron values
				data[t_data::PERIASTRON](n_radial, n_azimuthal) = atan2(e_y * cosFrameAngle + e_x * sinFrameAngle,e_x * cosFrameAngle - e_y * sinFrameAngle);
			} else {
				data[t_data::PERIASTRON](n_radial, n_azimuthal) = atan2(e_y,e_x);
			}
		}
	}
}

/**
	compute alpha gravitational

	alpha(R) = |d ln Omega/d ln R|^-1 (Tgrav)/(Sigma cs^2)
*/
void calculate_alpha_grav(t_data &data, unsigned int timestep, bool force_update)
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

	for (unsigned int n_radial = 0; n_radial <= data[t_data::ALPHA_GRAV].get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ALPHA_GRAV].get_max_azimuthal(); ++n_azimuthal) {
			/*data[t_data::ALPHA_GRAV](n_radial, n_azimuthal) = -2.0/3.0 * (data[t_data::T_GRAVITATIONAL](n_radial,n_azimuthal)+data[t_data::T_REYNOLDS](n_radial, n_azimuthal))/(data[t_data::DENSITY](n_radial,n_azimuthal)*pow2(data[t_data::SOUNDSPEED](n_radial, n_azimuthal)));*/			data[t_data::ALPHA_GRAV](n_radial, n_azimuthal) = 2.0/3.0 * data[t_data::T_GRAVITATIONAL](n_radial,n_azimuthal)/(data[t_data::DENSITY](n_radial,n_azimuthal)*pow2(data[t_data::SOUNDSPEED](n_radial, n_azimuthal)));
		}
	}
}

/**
	compute alpha gravitational

	alpha(R) = |d ln Omega/d ln R|^-1 (Tgrav)/(Sigma cs^2)
*/
void calculate_radial_alpha_grav(t_data &data, unsigned int timestep, bool force_update)
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

	// average gravitational stress
	data[t_data::T_GRAVITATIONAL_1D] = data[t_data::T_GRAVITATIONAL];
	// average density
	data[t_data::DENSITY_1D] = data[t_data::DENSITY];
	// average soundspeed
	data[t_data::SOUNDSPEED_1D] = data[t_data::SOUNDSPEED];

	for (unsigned int n_radial = 0; n_radial <= data[t_data::ALPHA_GRAV_1D].get_max_radial(); ++n_radial) {
 			data[t_data::ALPHA_GRAV_1D](n_radial) = 2.0/3.0 * data[t_data::T_GRAVITATIONAL_1D](n_radial)/(data[t_data::DENSITY_1D](n_radial)*pow2(data[t_data::SOUNDSPEED_1D](n_radial)));
	}
}

void calculate_alpha_grav_mean_sumup(t_data &data, unsigned int timestep, double dt)
{
	calculate_alpha_grav(data, timestep, true);

	for (unsigned int n_radial = 0; n_radial <= data[t_data::ALPHA_GRAV_MEAN].get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ALPHA_GRAV_MEAN].get_max_azimuthal(); ++n_azimuthal) {
			data[t_data::ALPHA_GRAV_MEAN](n_radial, n_azimuthal) += data[t_data::ALPHA_GRAV](n_radial, n_azimuthal) * dt;
		}
	}
}

void calculate_radial_alpha_grav_mean_sumup(t_data &data, unsigned int timestep, double dt)
{
	calculate_radial_alpha_grav(data, timestep, true);

	for (unsigned int n_radial = 0; n_radial <= data[t_data::ALPHA_GRAV_MEAN_1D].get_max_radial(); ++n_radial) {
		data[t_data::ALPHA_GRAV_MEAN_1D](n_radial) += data[t_data::ALPHA_GRAV_1D](n_radial) * dt;
	}
}

void calculate_alpha_grav_mean_finalize(t_data &data, double dt)
{
	for (unsigned int n_radial = 0; n_radial <= data[t_data::ALPHA_GRAV_MEAN].get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ALPHA_GRAV_MEAN].get_max_azimuthal(); ++n_azimuthal) {
			data[t_data::ALPHA_GRAV_MEAN](n_radial, n_azimuthal) /= dt;
		}
	}
}

void calculate_radial_alpha_grav_mean_finalize(t_data &data, double dt)
{
	for (unsigned int n_radial = 0; n_radial <= data[t_data::ALPHA_GRAV_MEAN_1D].get_max_radial(); ++n_radial) {
		data[t_data::ALPHA_GRAV_MEAN_1D](n_radial) /= dt;
	}
}

void calculate_alpha_grav_mean_reset(t_data &data)
{
	for (unsigned int n_radial = 0; n_radial <= data[t_data::ALPHA_GRAV_MEAN].get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ALPHA_GRAV_MEAN].get_max_azimuthal(); ++n_azimuthal) {
			data[t_data::ALPHA_GRAV_MEAN](n_radial, n_azimuthal) = 0;
		}
	}
}

void calculate_radial_alpha_grav_mean_reset(t_data &data)
{
	for (unsigned int n_radial = 0; n_radial <= data[t_data::ALPHA_GRAV_MEAN_1D].get_max_radial(); ++n_radial) {
		data[t_data::ALPHA_GRAV_MEAN_1D](n_radial) = 0;
	}
}

/**
	compute alpha Reynolds

	alpha(R) = |d ln Omega/d ln R|^-1 (Trey)/(Sigma cs^2)
*/
void calculate_alpha_reynolds(t_data &data, unsigned int timestep, bool force_update)
{
	static int last_timestep_calculated = -1;

	if (!force_update) {
		if (last_timestep_calculated == (int)timestep) {
			return;
		} else {
			last_timestep_calculated = timestep;
		}
	}

	stress::calculate_Reynolds_stress(data);

	// average gravitational stress
	data[t_data::T_REYNOLDS_1D] = data[t_data::T_REYNOLDS];
	// average density
	data[t_data::DENSITY_1D] = data[t_data::DENSITY];
	// average soundspeed
	data[t_data::SOUNDSPEED_1D] = data[t_data::SOUNDSPEED];

	for (unsigned int n_radial = 0; n_radial <= data[t_data::ALPHA_REYNOLDS_1D].get_max_radial(); ++n_radial) {
		data[t_data::ALPHA_REYNOLDS_1D](n_radial) = 2.0/3.0 * data[t_data::T_REYNOLDS_1D](n_radial)/(data[t_data::DENSITY_1D](n_radial)*pow2(data[t_data::SOUNDSPEED_1D](n_radial)));
	}
}

/**
	compute alpha Reynolds

	alpha(R) = |d ln Omega/d ln R|^-1 (Trey)/(Sigma cs^2)
*/
void calculate_radial_alpha_reynolds(t_data &data, unsigned int timestep, bool force_update)
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

	for (unsigned int n_radial = 0; n_radial <= data[t_data::ALPHA_REYNOLDS].get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ALPHA_REYNOLDS].get_max_azimuthal(); ++n_azimuthal) {
			data[t_data::ALPHA_REYNOLDS](n_radial, n_azimuthal) = 2.0/3.0 * data[t_data::T_REYNOLDS](n_radial,n_azimuthal)/(data[t_data::DENSITY](n_radial,n_azimuthal)*pow2(data[t_data::SOUNDSPEED](n_radial, n_azimuthal)));
		}
	}
}

void calculate_alpha_reynolds_mean_sumup(t_data &data, unsigned int timestep, double dt)
{
	calculate_alpha_reynolds(data, timestep, true);

	for (unsigned int n_radial = 0; n_radial <= data[t_data::ALPHA_REYNOLDS_MEAN].get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ALPHA_REYNOLDS_MEAN].get_max_azimuthal(); ++n_azimuthal) {
			data[t_data::ALPHA_REYNOLDS_MEAN](n_radial, n_azimuthal) += data[t_data::ALPHA_REYNOLDS](n_radial, n_azimuthal) * dt;
		}
	}
}

void calculate_radial_alpha_reynolds_mean_sumup(t_data &data, unsigned int timestep, double dt)
{
	calculate_radial_alpha_reynolds(data, timestep, true);

	for (unsigned int n_radial = 0; n_radial <= data[t_data::ALPHA_REYNOLDS_MEAN_1D].get_max_radial(); ++n_radial) {
		data[t_data::ALPHA_REYNOLDS_MEAN_1D](n_radial) += data[t_data::ALPHA_REYNOLDS_1D](n_radial) * dt;
	}
}

void calculate_alpha_reynolds_mean_finalize(t_data &data, double dt)
{
	for (unsigned int n_radial = 0; n_radial <= data[t_data::ALPHA_REYNOLDS_MEAN].get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ALPHA_REYNOLDS_MEAN].get_max_azimuthal(); ++n_azimuthal) {
			data[t_data::ALPHA_REYNOLDS_MEAN](n_radial, n_azimuthal) /= dt;
		}
	}
}

void calculate_radial_alpha_reynolds_mean_finalize(t_data &data, double dt)
{
	for (unsigned int n_radial = 0; n_radial <= data[t_data::ALPHA_REYNOLDS_MEAN_1D].get_max_radial(); ++n_radial) {
		data[t_data::ALPHA_REYNOLDS_MEAN_1D](n_radial) /= dt;
	}
}

void calculate_alpha_reynolds_mean_reset(t_data &data)
{
	for (unsigned int n_radial = 0; n_radial <= data[t_data::ALPHA_REYNOLDS_MEAN].get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::ALPHA_REYNOLDS_MEAN].get_max_azimuthal(); ++n_azimuthal) {
			data[t_data::ALPHA_REYNOLDS_MEAN](n_radial, n_azimuthal) = 0;
		}
	}
}

void calculate_radial_alpha_reynolds_mean_reset(t_data &data)
{
	for (unsigned int n_radial = 0; n_radial <= data[t_data::ALPHA_REYNOLDS_MEAN_1D].get_max_radial(); ++n_radial) {
		data[t_data::ALPHA_REYNOLDS_MEAN_1D](n_radial) = 0;
	}
}

/**
	Calculates Toomre Q parameter
*/
void calculate_toomre(t_data &data, unsigned int /* timestep */, bool /* force_update */)
{
	double kappa;

	for (unsigned int n_radial = 1; n_radial <= data[t_data::TOOMRE].get_max_radial(); ++n_radial) {
		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::TOOMRE].get_max_azimuthal(); ++n_azimuthal) {
			// kappa^2 = 1/r^3 d((r^2 Omega)^2)/dr = 1/r^3 d((r*v_phi)^2)/dr
			kappa = sqrt(fabs(pow3(InvRmed[n_radial])*(pow2(data[t_data::V_AZIMUTHAL](n_radial,n_azimuthal)*Rmed[n_radial])-pow2(data[t_data::V_AZIMUTHAL](n_radial-1,n_azimuthal)*Rmed[n_radial-1]))*InvDiffRmed[n_radial]));

			// Q = (c_s kappa) / (Pi G Sigma)
			//data[t_data::TOOMRE](n_radial, n_azimuthal) = data[t_data::SOUNDSPEED](n_radial, n_azimuthal)*omega_kepler(Rmed[n_radial])/(PI*G*data[t_data::DENSITY](n_radial, n_azimuthal));
			data[t_data::TOOMRE](n_radial, n_azimuthal) = data[t_data::SOUNDSPEED](n_radial, n_azimuthal)*kappa/(PI*constants::G*data[t_data::DENSITY](n_radial, n_azimuthal));
		}
	}
}

/**
	Calculates Toomre Q parameter
*/
void calculate_radial_toomre(t_data &data, unsigned int /* timestep */, bool /* force_update */)
{
	// calculate averaged density
	data[t_data::DENSITY_1D] = data[t_data::DENSITY];
	// calculate averaged soundspeed
	data[t_data::SOUNDSPEED_1D] = data[t_data::SOUNDSPEED];
	// calculate averaged v_azimuthal
	data[t_data::V_AZIMUTHAL_1D] = data[t_data::V_AZIMUTHAL];

	double kappa;

	for (unsigned int n_radial = 1; n_radial <= data[t_data::TOOMRE_1D].get_max_radial(); ++n_radial) {
		// kappa^2 = 1/r^3 d((r^2 Omega)^2)/dr = 1/r^3 d((r*v_phi)^2)/dr
		kappa = sqrt(fabs(pow3(InvRmed[n_radial])*(pow2(data[t_data::V_AZIMUTHAL_1D](n_radial)*Rmed[n_radial])-pow2(data[t_data::V_AZIMUTHAL_1D](n_radial-1)*Rmed[n_radial-1]))*InvDiffRmed[n_radial]));

		// Q = (c_s kappa) / (Pi G Sigma)
		data[t_data::TOOMRE_1D](n_radial) = data[t_data::SOUNDSPEED_1D](n_radial)*kappa/(PI*constants::G*data[t_data::DENSITY_1D](n_radial));
	}
}

void calculate_radial_luminosity(t_data &data, unsigned int timestep, bool force_update)
{
	static int last_timestep_calculated = -1;

	if ((!force_update) && (last_timestep_calculated == (int)timestep)) {
		return;
	}

	last_timestep_calculated = timestep;

	for (unsigned int n_radial = 0; n_radial <= data[t_data::LUMINOSITY_1D].get_max_radial(); ++n_radial) {
		// L = int( int(sigma T^4 r ,phi) ,r);
		data[t_data::LUMINOSITY_1D](n_radial) = 0.0;

		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::QMINUS].get_max_azimuthal(); ++n_azimuthal) {
			double dr = (Rsup[n_radial]-Rinf[n_radial]);
			double dphi = 2.0*PI/(double)data[t_data::QMINUS].get_size_azimuthal();
			data[t_data::LUMINOSITY_1D](n_radial) += data[t_data::QMINUS](n_radial, n_azimuthal)*Rmed[n_radial]*dr*dphi;
		}
	}
}

void calculate_radial_dissipation(t_data &data, unsigned int timestep, bool force_update)
{
	static int last_timestep_calculated = -1;

	if ((!force_update) && (last_timestep_calculated == (int)timestep)) {
		return;
	}

	last_timestep_calculated = timestep;

	for (unsigned int n_radial = 0; n_radial <= data[t_data::DISSIPATION_1D].get_max_radial(); ++n_radial) {
		data[t_data::DISSIPATION_1D](n_radial) = 0.0;

		for (unsigned int n_azimuthal = 0; n_azimuthal <= data[t_data::QPLUS].get_max_azimuthal(); ++n_azimuthal) {
			double dr = (Rsup[n_radial]-Rinf[n_radial]);
			double dphi = 2.0*PI/(double)data[t_data::QPLUS].get_size_azimuthal();
			data[t_data::DISSIPATION_1D](n_radial) += data[t_data::QPLUS](n_radial, n_azimuthal)*Rmed[n_radial]*dr*dphi;
		}
	}
}

void calculate_massflow(t_data &data, unsigned int timestep, bool force_update)
{
    (void) timestep;
    (void) force_update;

	// divide the data in massflow by the large timestep DT before writing out
	// to obtain the massflow from the mass difference
    for (unsigned int nRadial=0; nRadial<data[t_data::MASSFLOW_1D].get_size_radial(); ++nRadial) {
        data[t_data::MASSFLOW_1D](nRadial) *= 1./DT;
    }
}

}
