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
    for (unsigned int k = 1;
	 k < data.get_planetary_system().get_number_of_planets(); ++k) {

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
	const unsigned int Nphi = data[t_data::SIGMA].get_size_azimuthal();

	for (unsigned int nr = radial_first_active; nr < radial_active_size; ++nr) {
		for (unsigned int naz = 0; naz < Nphi; ++naz) {
		unsigned int cell = get_cell_id(nr, naz);
		const double dist = std::sqrt(
		    (cell_center_x[cell] - xpl) * (cell_center_x[cell] - xpl) +
		    (cell_center_y[cell] - ypl) * (cell_center_y[cell] - ypl));
		if (dist < roche_radius) {
			mdcplocal += Surf[nr] *
				 data[t_data::SIGMA](nr, naz);
		}
	    }
	}

	double mdcptotal = 0.0;
	MPI_Allreduce(&mdcplocal, &mdcptotal, 1, MPI_DOUBLE, MPI_SUM,
		      MPI_COMM_WORLD);

	planet.set_circumplanetary_mass(mdcptotal);
    }
}