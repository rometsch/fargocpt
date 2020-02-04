#include "find_cell_id.h"
#include "global.h"
#include "parameters.h"
#include "constants.h"
#include <cmath>
#include "LowTasks.h"
#include "logging.h"

static double growth_factor;
static double log_growth_factor;
static double optimization_const;
static double cell_size;
static double inv_phi_cell_size;
static double phi_cell_size;


unsigned int clamp_r_id_to_rmed_grid(int cell_id)
{

	const int Nr_max = (int)NRadial-2;
	if(cell_id < 0)
	{
		cell_id = 0;
	}
	else if(cell_id > Nr_max)
	{
		cell_id = Nr_max;
	}

	return (unsigned int)cell_id;
}

unsigned int clamp_r_id_to_radii_grid(int cell_id)
{

	const int Nr_max = (int)NRadial-1;
	if(cell_id < 0)
	{
		cell_id = 0;
	}
	else if(cell_id > Nr_max)
	{
		cell_id = Nr_max;
	}

	return (unsigned int)cell_id;
}


unsigned int clamp_phi_id_to_grid(int cell_id)
{
	if(cell_id < 0)
	{
		cell_id += NAzimuthal;
	}
	if(cell_id >= (int)NAzimuthal)
	{
		cell_id -= NAzimuthal;
	}

	return (unsigned int)cell_id;
}

static void init_cell_finder_log(const double cell_growth_factor, const double first_cell_size){
	phi_cell_size = 2.0 * PI / (double)NAzimuthal;
	inv_phi_cell_size = (double)NAzimuthal / 2.0 / PI;

	growth_factor = cell_growth_factor;
	optimization_const = 3.0/2.0 / RMIN * (1-std::pow(growth_factor, 2.0)) / (1-std::pow(growth_factor, 3.0));
	log_growth_factor = std::log(growth_factor);

	// unused
	cell_size = 0.0;
	(void)first_cell_size;
	return;
}

static void init_cell_finder_arith(const double cell_growth_factor, const double first_cell_size){
	phi_cell_size = 2.0 * PI / (double)NAzimuthal;
	inv_phi_cell_size = (double)NAzimuthal / 2.0 / PI;
	growth_factor = cell_growth_factor;

	// unused
	(void)first_cell_size;
	log_growth_factor = 0.0;
	cell_size = 0.0;
}

static void init_cell_finder_exp(const double cell_growth_factor, const double first_cell_size){
	phi_cell_size = 2.0 * PI / (double)NAzimuthal;
	inv_phi_cell_size = (double)NAzimuthal / 2.0 / PI;
	cell_size = first_cell_size;
	growth_factor = cell_growth_factor;
	log_growth_factor = std::log(growth_factor);
	optimization_const = (growth_factor - 1.0)/cell_size;
}

static void init_cell_finder_custom(const double cell_growth_factor, const double first_cell_size){
	phi_cell_size = 2.0 * PI / (double)NAzimuthal;
	inv_phi_cell_size = (double)NAzimuthal / 2.0 / PI;
	// unused
	(void)cell_growth_factor;
	(void)first_cell_size;
	cell_size = 0.0;
	growth_factor = 0.0;
	log_growth_factor = 0.0;
	optimization_const = 0.0;
}



static int get_rmed_id_log(const double r){

		double did = log(r * optimization_const) / log_growth_factor;
		int id = (int)std::floor(did) - IMIN + 1;

#ifndef NDEBUG
		{
			bool lower_than_ceil = r < Rmed[id+1];
			bool higher_than_floor = r > Rmed[id];

			if(!(lower_than_ceil && higher_than_floor))
			{
				logging::print("Error: get_rmed_id has failed for %s grid!\n", parameters::radial_grid_names[parameters::radial_grid_type]);
				logging::print("Found id=%d for r=%.5e but check Rmed[%d] = %.5e < %.5e	< %.5e = Rmed[%d] is not true!\n",
							   id, r, id, Rmed[id], r, Rmed[id+1], id+1);
				die("END!\n");
			}
		}
#endif // DEBUG
		return id;

}

static int get_rmed_id_arith(const double r){
		double did = (r-RMIN) * growth_factor;
		int id = (int)std::floor(did) - IMIN + 1;
		if(Rmed[id] > r)
		{
			id--;
		}

#ifndef NDEBUG
		{
			bool lower_than_ceil = r < Rmed[id+1];
			bool higher_than_floor = r > Rmed[id];

			if(!(lower_than_ceil && higher_than_floor))
			{
				logging::print("Error: get_rmed_id has failed for %s grid!\n", parameters::radial_grid_names[parameters::radial_grid_type]);
				logging::print("Found id=%d for r=%.5e but check Rmed[%d] = %.5e < %.5e	< %.5e = Rmed[%d] is not true!\n",
							   id, r, id, Rmed[id], r, Rmed[id+1], id+1);
				die("END!\n");
			}
		}
#endif // DEBUG

		return id;
}

static int get_rmed_id_exp(const double r){
		double tmp = (r - RMIN) * optimization_const + 1.0;
		double did = std::log(tmp)/log_growth_factor;
		int id = std::floor(did) - IMIN + 1;

		if(Rmed[id] > r)
		{
			id--;
		}

#ifndef NDEBUG
		{
			bool lower_than_ceil = r < Rmed[id+1];
			bool higher_than_floor = r > Rmed[id];

			if(!(lower_than_ceil && higher_than_floor))
			{
				logging::print("Error: get_rmed_id has failed for %s grid!\n", parameters::radial_grid_names[parameters::radial_grid_type]);
				logging::print("Found id=%d for r=%.5e but check Rmed[%d] = %.5e < %.5e	< %.5e = Rmed[%d] is not true!\n",
							   id, r, id, Rmed[id], r, Rmed[id+1], id+1);
				die("END!\n");
			}
		}
#endif // DEBUG

		return id;
}

static int get_rmed_id_custom(const double r){
			int id = 0;
			while (Rmed[id] < r)
			{
				id++;
			}
			id--;

#ifndef NDEBUG
			{
				bool lower_than_ceil = r < Rmed[id+1];
				bool higher_than_floor = r > Rmed[id];

				if(!(lower_than_ceil && higher_than_floor))
				{
					logging::print("Error: get_rmed_id has failed for %s grid!\n", parameters::radial_grid_names[parameters::radial_grid_type]);
					logging::print("Found id=%d for r=%.5e but check Rmed[%d] = %.5e < %.5e	< %.5e = Rmed[%d] is not true!\n",
								   id, r, id, Rmed[id], r, Rmed[id+1], id+1);
					die("END!\n");
				}
			}
#endif // DEBUG

			return id;
}


static int get_rinf_id_log(const double r){
		double did = log(r / RMIN) / log_growth_factor;
		int id = (int)std::floor(did) - IMIN + 1;

#ifndef NDEBUG
		{
			bool lower_than_ceil = r < Rsup[id];
			bool higher_than_floor = r > Rinf[id];

			if(!(lower_than_ceil && higher_than_floor))
			{
				logging::print("Error: get_rinf_id_log has failed for %s grid!\n", parameters::radial_grid_names[parameters::radial_grid_type]);
				logging::print("Found id=%d for r=%.5e but check Rinf[%d] = %.5e < %.5e	< %.5e = Rsup[%d] is not true	%d	%d!\n",
							   id, r, id, Rinf[id], r, Rsup[id], id, higher_than_floor, lower_than_ceil);
				die("END!\n");
			}
		}
#endif // DEBUG
		return id;
}

static int get_rinf_id_arith(const double r){
		double did = (r-RMIN) * growth_factor;
		int id = (int)std::floor(did) - IMIN + 1;

#ifndef NDEBUG
		{
			bool lower_than_ceil = r < Rsup[id];
			bool higher_than_floor = r > Rinf[id];

			if(!(lower_than_ceil && higher_than_floor))
			{
				logging::print("Error: get_rinf_id_arith has failed for %s grid!\n", parameters::radial_grid_names[parameters::radial_grid_type]);
				logging::print("Found id=%d for r=%.5e but check Rinf[%d] = %.5e < %.5e	< %.5e = Rsup[%d] is not true!\n",
							   id, r, id, Rinf[id], r, Rsup[id], id);
				die("END!\n");
			}
		}
#endif // DEBUG
		return id;
}

static int get_rinf_id_exp(const double r){
		double tmp = (r - RMIN) * optimization_const + 1.0;
		double did = std::log(tmp)/log_growth_factor;
		int id = std::floor(did) - IMIN + 1;
#ifndef NDEBUG
		{
			bool lower_than_ceil = r < Rsup[id];
			bool higher_than_floor = r > Rinf[id];

			if(!(lower_than_ceil && higher_than_floor))
			{
				logging::print("Error: get_rinf_id_exp has failed for %s grid!\n", parameters::radial_grid_names[parameters::radial_grid_type]);
				logging::print("Found id=%d for r=%.5e but check Rinf[%d] = %.5e < %.5e	< %.5e = Rsup[%d] is not true!\n",
							   id, r, id, Rinf[id], r, Rsup[id], id);
				die("END!\n");
			}
		}
#endif // DEBUG
		return id;
}

static int get_rinf_id_custom(const double r){
	int id = 0;
	while (Rinf[id] < r)
	{
		id++;
	}
	id--;

#ifndef NDEBUG
	{
		bool lower_than_ceil = r < Rsup[id];
		bool higher_than_floor = r > Rinf[id];

		if(!(lower_than_ceil && higher_than_floor))
		{
			logging::print("Error: get_rinf_id_custom has failed for %s grid!\n", parameters::radial_grid_names[parameters::radial_grid_type]);
			logging::print("Found id=%d for r=%.5e but check Rinf[%d] = %.5e < %.5e	< %.5e = Rsup[%d] is not true!\n",
						   id, r, id, Rinf[id], r, Rsup[id], id);
			die("END!\n");
		}
	}
#endif // DEBUG

	return id;
}

unsigned int get_next_azimuthal_id(const unsigned int id_low)
{
	unsigned int id_high = id_low + 1;
	if(id_high == NAzimuthal)
	{
		id_high = 0;
	}

	return id_high;
}

int get_inf_azimuthal_id(const double phi){

	double did = phi * inv_phi_cell_size;
	int id_low = (int)std::floor(did);


#ifndef NDEBUG
	{
		bool lower_than_ceil = phi < phi_cell_size*(double)(id_low+1);
		bool higher_than_floor = phi > phi_cell_size*(double)(id_low);

		if(!(lower_than_ceil && higher_than_floor))
		{
			logging::print("Error: get_inf_azimuthal_id has failed for %s grid!\n", parameters::radial_grid_names[parameters::radial_grid_type]);
			logging::print("Found id=%d for phi=%.5e but check Phi[%d] = %.5e < %.5e	< %.5e = Phi[%d] is not true!\n",
						   id_low+1, phi, id_low, phi_cell_size*(double)id_low, phi, phi_cell_size*(double)(id_low+1), id_low+1);
			die("END!\n");
		}
	}
#endif // DEBUG

	return id_low;

}

int get_med_azimuthal_id(const double phi){

	double did = ((phi - 0.5*phi_cell_size) * inv_phi_cell_size);
	int id_low = (int)std::floor(did);


#ifndef NDEBUG
	{
		bool lower_than_ceil = phi - 0.5*phi_cell_size < phi_cell_size*(double)(id_low+1);

		bool higher_than_floor = phi - 0.5*phi_cell_size > phi_cell_size*(double)(id_low);

		if(!(lower_than_ceil && higher_than_floor))
		{
			logging::print("Error: get_med_azimuthal_id has failed for %s grid!\n", parameters::radial_grid_names[parameters::radial_grid_type]);
			logging::print("Found id=%d for phi=%.5e but check Phi[%d] = %.5e < %.5e	< %.5e = Phi[%d] is not true!	%.5e\n",
						   id_low, phi, id_low, phi_cell_size*(double)id_low + 0.5*phi_cell_size, phi, phi_cell_size*(double)(id_low + 1) + 0.5*phi_cell_size, id_low+1, did);
			die("END!\n");
		}
	}
#endif // DEBUG

	return id_low;

}



void init_cell_finder(parameters::t_radial_grid GRID, const double cell_growth_factor, const double first_cell_size)
{
	switch(GRID)
	{
		case parameters::t_radial_grid::logarithmic_spacing:
		{
			init_cell_finder_log(cell_growth_factor, first_cell_size);
			break;
		}

		case parameters::t_radial_grid::exponential_spacing:
		{
			init_cell_finder_exp(cell_growth_factor, first_cell_size);
			break;
		}

		case parameters::t_radial_grid::arithmetic_spacing:
		{
			init_cell_finder_arith(cell_growth_factor, first_cell_size);
			break;
		}

		case parameters::t_radial_grid::custom_spacing:
		{
			init_cell_finder_custom(cell_growth_factor, first_cell_size);
			break;
		}
		default:
		{
			die("Error in init_cell_finder: Invalid GRID option!\n");
			break;
		}
	}
}


int get_rmed_id(parameters::t_radial_grid GRID, const double r)
{

	int id;
	switch(GRID)
	{
		case parameters::t_radial_grid::logarithmic_spacing:
		{
			id = get_rmed_id_log(r);
			break;
		}

		case parameters::t_radial_grid::exponential_spacing:
		{
			id = get_rmed_id_exp(r);
			break;
		}

		case parameters::t_radial_grid::arithmetic_spacing:
		{
			id = get_rmed_id_arith(r);
			break;
		}

		case parameters::t_radial_grid::custom_spacing:
		{
			id = get_rmed_id_custom(r);
			break;
		}
		default:
		{
			die("Error in get_rmed_id: Invalid GRID option!\n");
			id = -1;
			break;
		}
	}

	return id;
}

int get_rinf_id(parameters::t_radial_grid GRID, const double r)
{
	int id;
	switch(GRID)
	{
		case parameters::t_radial_grid::logarithmic_spacing:
		{
			id = get_rinf_id_log(r);
			break;
		}

		case parameters::t_radial_grid::exponential_spacing:
		{
			id = get_rinf_id_exp(r);
			break;
		}

		case parameters::t_radial_grid::arithmetic_spacing:
		{
			id = get_rinf_id_arith(r);
			break;
		}

		case parameters::t_radial_grid::custom_spacing:
		{
			id = get_rinf_id_custom(r);
			break;
		}
		default:
		{
			die("Error in get_ring_id: Invalid GRID option!\n");
			id = -1;
			break;
		}
	}

	return id;
}