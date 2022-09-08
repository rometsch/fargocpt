#include "types.h"
#include "mpi.h"
#include "circumplanetary_mass.h"
#include "global.h"
#include "data.h"
#include "util.h"



/**
	Calculates the gas mass inside the planets Roche lobe
*/
void ComputeCircumPlanetaryMasses(t_data &data)
{
	const unsigned int Npl = data.get_planetary_system().get_number_of_planets();
    for (unsigned int k = 1; k < Npl; ++k) {

	// TODO: non global
	const double *cell_center_x = CellCenterX->Field;
	const double *cell_center_y = CellCenterY->Field;

	auto &planet = data.get_planetary_system().get_planet(k);
	const double planet_to_prim_dist = planet.get_distance_to_primary();
	const double roche_radius =
	    planet_to_prim_dist * planet.get_dimensionless_roche_radius();

	const double xpl = planet.get_x();
	const double ypl = planet.get_y();

	double mdcplocal = 0.0;

	for (unsigned int n_radial = radial_first_active;
	     n_radial < radial_active_size; ++n_radial) {
	    for (unsigned int n_azimuthal = 0;
		 n_azimuthal < data[t_data::SIGMA].get_size_azimuthal();
		 ++n_azimuthal) {
		unsigned int cell = get_cell_id(n_radial, n_azimuthal);
		const double dist = std::sqrt(
		    (cell_center_x[cell] - xpl) * (cell_center_x[cell] - xpl) +
		    (cell_center_y[cell] - ypl) * (cell_center_y[cell] - ypl));
		if (dist < roche_radius) {
		    mdcplocal += Surf[n_radial] *
				 data[t_data::SIGMA](n_radial, n_azimuthal);
		}
	    }
	}

	double mdcptotal = 0.0;
	MPI_Allreduce(&mdcplocal, &mdcptotal, 1, MPI_DOUBLE, MPI_SUM,
		      MPI_COMM_WORLD);

	planet.set_circumplanetary_mass(mdcptotal);
    }
}