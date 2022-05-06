/**
	\file accretion.cpp
	\author Thomas Rometsch <thomas.rometsch@uni-tuebingen.de>

	Handling accretion onto planets and stars.
*/

#include <cmath>
#include <tuple>

#include "Theo.h"
#include "accretion.h"
#include "constants.h"
#include "find_cell_id.h"
#include "global.h"
#include "parameters.h"
#include "planet.h"
#include "viscosity.h"

namespace accretion
{
static std::tuple<unsigned int, unsigned int>
hill_radial_index(const double Rplanet, const double RHill)
{

    const bool is_vector = false;
    /* Calculate the indeces in radial direction where
       the Hill sphere starts and stops */
    unsigned i_min =
	clamp_r_id_to_radii_grid(get_rinf_id(Rplanet - RHill), is_vector);
    unsigned i_max =
	clamp_r_id_to_radii_grid(get_rinf_id(Rplanet + RHill) + 1, is_vector);
    std::tuple<unsigned int, unsigned int> ids(i_min, i_max);
    return ids;
}

static std::tuple<int, int> hill_azimuthal_index(const double angle,
						 const double Rplanet,
						 const double RHill)
{
    /* Calculate the index in azimuthal direction
       where the Hill sphere starts and stops */
    const int i_min = get_med_azimuthal_id(angle - 2.0 * RHill / Rplanet);
    const int i_max = get_med_azimuthal_id(angle + 2.0 * RHill / Rplanet) + 1;
    std::tuple<int, int> ids(i_min, i_max);
    return ids;
}

static void update_planet(t_planet &planet, const double dMplanet,
			  const double dPxPlanet, const double dPyPlanet)
{
    /* Update a planets mass and velocities. */
    const double VXplanet = planet.get_vx();
    const double VYplanet = planet.get_vy();

    double Mplanet = planet.get_mass();
    double PxPlanet = Mplanet * VXplanet;
    double PyPlanet = Mplanet * VYplanet;

    Mplanet += dMplanet;
    PxPlanet += dPxPlanet;
    PyPlanet += dPyPlanet;

    planet.set_vx(PxPlanet / Mplanet);
    planet.set_vy(PyPlanet / Mplanet);
    planet.set_mass(Mplanet);
}

static bool AccreteOntoSinglePlanet(t_data &data, t_planet &planet, double dt)
{
    bool mass_changed = false;
    const int ns = data[t_data::DENSITY].Nsec;
    double *dens = data[t_data::DENSITY].Field;
    double *energy = data[t_data::ENERGY].Field;
    const double *cell_center_x = CellCenterX->Field;
    const double *cell_center_y = CellCenterY->Field;
    const double *vrad = data[t_data::V_RADIAL].Field;
    const double *vtheta = data[t_data::V_AZIMUTHAL].Field;
    const double density_floor = parameters::sigma_floor * parameters::sigma0;

    const double Xplanet = planet.get_x();
    const double Yplanet = planet.get_y();
    const double Rplanet = planet.get_r();

    // Hereafter : initialization of W. Kley's parameters
    // remove a ratio of facc = planet.get_acc() of the mass inside the
    // Hill sphere every planet orbit
    // const double facc = dt * planet.get_acc() / std::pow(Rplanet, 2) /
    // planet.get_period() * std::log(2);
    const double facc =
	dt * planet.get_acc() / planet.get_orbital_period() * std::log(2);

    const double facc1 = 1.0 / 3.0 * facc;
    const double facc2 = 2.0 / 3.0 * facc;
	const double frac1 = parameters::accretion_radius_fraction;
	const double frac2 = 0.5 * parameters::accretion_radius_fraction;
    // W. Kley's parameters initialization finished

	const double RHill = planet.get_dimensionless_roche_radius() *
			 planet.get_distance_to_primary();
    // search radius is bigger fraction + 2 dphi cell sizes to capture all cells
    const double search_radius = RHill * frac1 + 2.0 * Rplanet / ns;

    // calculate range of indeces to iterate over
    const auto [i_min, i_max] = hill_radial_index(Rplanet, search_radius);
    const double angle = planet.get_phi();
    const auto [j_min, j_max] =
	hill_azimuthal_index(angle, Rplanet, search_radius);

    double dMplanet = 0.0;
    double dPxPlanet = 0.0;
    double dPyPlanet = 0.0;

    for (unsigned int i = i_min; i <= i_max; i++) {
	for (int j = j_min; j <= j_max; j++) {
	    // map azimuthal index to [0, ns]
	    int jf = clamp_phi_id_to_grid(j);
	    ;
	    // calculate cell 1d index
	    int l = jf + i * ns;
	    int lip = l + ns;
	    int ljp = l + 1;
	    if (jf == ns - 1) {
		ljp = i * ns;
	    }

	    const double xc = cell_center_x[l];
	    const double yc = cell_center_y[l];
	    const double dx = Xplanet - xc;
	    const double dy = Yplanet - yc;
	    const double distance = sqrt(dx * dx + dy * dy);

	    // interpolate velocities to cell centers
	    const double vtcell =
		0.5 * (vtheta[l] + vtheta[ljp]) + Rmed[i] * OmegaFrame;
	    const double vrcell = 0.5 * (vrad[l] + vrad[lip]);
	    // calculate cartesian velocities
	    const double vxcell = (vrcell * xc - vtcell * yc) / Rmed[i];
	    const double vycell = (vrcell * yc + vtcell * xc) / Rmed[i];

	    // only allow removal of mass down to density floor
	    const double facc_max = 1 - density_floor / dens[l];
	    // handle accretion zone 1
	    if (distance < frac1 * RHill) {
		const double facc_ceil = std::min(facc1, facc_max);
		const double deltaM = facc_ceil * dens[l] * Surf[i];
		dens[l] *= 1.0 - facc_ceil;
		if (parameters::Adiabatic) {
		    energy[l] *= 1.0 - facc_ceil;
		}
		if (radial_first_active < i &&
		    i < radial_active_size) { // Only add active cells to
					      // planet
		    dPxPlanet += deltaM * vxcell;
		    dPyPlanet += deltaM * vycell;
		    dMplanet += deltaM;
		}
	    }
	    // handle accretion zone 2
	    if (distance < frac2 * RHill) {
		const double facc_ceil = std::min(facc2, facc_max);
		const double deltaM = facc_ceil * dens[l] * Surf[i];
		dens[l] *= 1.0 - facc_ceil;
		if (parameters::Adiabatic) {
		    energy[l] *= 1.0 - facc2;
		}
		if (radial_first_active < i &&
		    i < radial_active_size) { // Only add active cells to
					      // planet
		    dPxPlanet += deltaM * vxcell;
		    dPyPlanet += deltaM * vycell;
		    dMplanet += deltaM;
		}
	    }
	}
    }

    // MPI reduce
    double temp;
    MPI_Allreduce(&dMplanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dMplanet = temp;

    // monitoring purpose only
    planet.add_accreted_mass(dMplanet);

    if (parameters::disk_feedback) { // only update planets if they feel the
				     // disk
	MPI_Allreduce(&dPxPlanet, &temp, 1, MPI_DOUBLE, MPI_SUM,
		      MPI_COMM_WORLD);
	dPxPlanet = temp;
	MPI_Allreduce(&dPyPlanet, &temp, 1, MPI_DOUBLE, MPI_SUM,
		      MPI_COMM_WORLD);
	dPyPlanet = temp;

	// update planet momentum
	update_planet(planet, dMplanet, dPxPlanet, dPyPlanet);
	mass_changed = dMplanet > 0;
    }
    return mass_changed;
}

static bool AccreteOntoSinglePlanetViscous(t_data &data, t_planet &planet,
					   double dt)
{
    bool mass_changed = false;
    const int ns = data[t_data::DENSITY].Nsec;
    double *dens = data[t_data::DENSITY].Field;
    double *energy = data[t_data::ENERGY].Field;
    const double *cell_center_x = CellCenterX->Field;
    const double *cell_center_y = CellCenterY->Field;
    const double *vrad = data[t_data::V_RADIAL].Field;
    const double *vtheta = data[t_data::V_AZIMUTHAL].Field;
    const double density_floor = parameters::sigma_floor * parameters::sigma0;

    const double Xplanet = planet.get_x();
    const double Yplanet = planet.get_y();
    const double Rplanet = planet.get_r();

    // remove a ratio of facc = planet.get_acc() of the mass inside the
    // Hill sphere every planet orbit
    // we use M_dot = 3 pi nu Sigma / (1 - sqrt(R_in / R))
    // to derive the fraction we need to remove
    const double facc = dt * 3.0 * M_PI * parameters::viscous_outflow_speed;

	const double frac = parameters::accretion_radius_fraction;

    const double RHill = planet.get_dimensionless_roche_radius() *
			 planet.get_distance_to_primary();
    // search radius is bigger fraction + 2 dphi cell sizes to capture all cells
    const double search_radius = RHill * frac + 2.0 * Rplanet / ns;

    const double dist_max = RHill * frac;
    const double f_const = 3.0 / M_PI / std::pow(dist_max, 2);

    // calculate range of indeces to iterate over
    const auto [i_min, i_max] = hill_radial_index(Rplanet, search_radius);
    const double angle = planet.get_phi();
    const auto [j_min, j_max] =
	hill_azimuthal_index(angle, Rplanet, search_radius);

    double dMplanet = 0.0;
    double dPxPlanet = 0.0;
    double dPyPlanet = 0.0;

    for (unsigned int i = i_min; i <= i_max; i++) {
	for (int j = j_min; j <= j_max; j++) {
	    // map azimuthal index to [0, ns]
	    int jf = clamp_phi_id_to_grid(j);

	    // calculate cell 1d index
	    int l = jf + i * ns;
	    int lip = l + ns;
	    int ljp = l + 1;
	    if (jf == ns - 1) {
		ljp = i * ns;
	    }

	    const double xc = cell_center_x[l];
	    const double yc = cell_center_y[l];
	    const double dx = Xplanet - xc;
	    const double dy = Yplanet - yc;
	    const double distance = sqrt(dx * dx + dy * dy);
	    const double nu = data[t_data::VISCOSITY].Field[l];

	    const double spread =
		f_const * (1.0 - distance / dist_max) * Surf[i];

	    // debug variables
	    // static double acc_max = 0.0;
	    // const double acc = 3.0 * nu * M_PI * spread *
	    // parameters::viscous_outflow_speed * planet.get_period() /
	    // std::log(2); if(acc > acc_max) 	acc_max = acc;

	    // interpolate velocities to cell centers
	    const double vtcell =
		0.5 * (vtheta[l] + vtheta[ljp]) + Rmed[i] * OmegaFrame;
	    const double vrcell = 0.5 * (vrad[l] + vrad[lip]);
	    // calculate cartesian velocities
	    const double vxcell = (vrcell * xc - vtcell * yc) / Rmed[i];
	    const double vycell = (vrcell * yc + vtcell * xc) / Rmed[i];

	    // only allow removal of mass down to density floor
	    const double facc_max_dens = 1 - density_floor / dens[l];
	    // handle accretion zone 1
	    if (distance < frac * RHill) {
		// printf("Cell1 (%d %d) acc = %.3e	(%.3e)	1/acc = %.3e
		// %.3e	i(%d %d) j(%d %d)\n", i, j, acc,
		// planet.get_acc(), 1.0/acc, 1.0/acc_max, i_min, i_max, j_min,
		// j_max);

		const double facc_tmp = facc * nu * spread;
		const double facc_ceil = std::min(facc_tmp, facc_max_dens);
		const double deltaM = facc_ceil * dens[l] * Surf[i];
		dens[l] *= 1.0 - facc_ceil;
		if (parameters::Adiabatic) {
		    // We reduce energy by the same fraction as density,
		    // so temperature should not change and this is safe to do.
		    energy[l] *= 1.0 - facc_ceil;
		}
		if (radial_first_active < i &&
		    i < radial_active_size) { // Only add active cells to
					      // planet
		    dPxPlanet += deltaM * vxcell;
		    dPyPlanet += deltaM * vycell;
		    dMplanet += deltaM;
		}
	    }
	}
    }

    // MPI reduce
    double temp;
    MPI_Allreduce(&dMplanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    dMplanet = temp;

    // monitoring purpose only
    planet.add_accreted_mass(dMplanet);

    if (parameters::disk_feedback) { // only update planets if they feel the
				     // disk
	MPI_Allreduce(&dPxPlanet, &temp, 1, MPI_DOUBLE, MPI_SUM,
		      MPI_COMM_WORLD);
	dPxPlanet = temp;
	MPI_Allreduce(&dPyPlanet, &temp, 1, MPI_DOUBLE, MPI_SUM,
		      MPI_COMM_WORLD);
	dPyPlanet = temp;

	// update planet momentum
	update_planet(planet, dMplanet, dPxPlanet, dPyPlanet);
	mass_changed = dMplanet > 0;
    }
    return mass_changed;
}

void AccreteOntoPlanets(t_data &data, const double dt)
{
    bool masses_changed = false;

    auto &planetary_system = data.get_planetary_system();
    for (unsigned int k = 0; k < planetary_system.get_number_of_planets();
	 k++) {

	auto &planet = planetary_system.get_planet(k);
	if (planet.get_acc() > 1e-10) {
	    const bool changed = AccreteOntoSinglePlanet(data, planet, dt);
	    masses_changed = masses_changed || changed;
	}

	// Negative accretion enables viscous accretion
	if (planet.get_acc() < 0.0) {
	    const bool changed =
		AccreteOntoSinglePlanetViscous(data, planet, dt);
	    masses_changed = masses_changed || changed;
	}
    }

    // update hydro center mass
    if (masses_changed) {
	planetary_system.update_global_hydro_frame_center_mass();
	planetary_system.update_roche_radii();
    }
}

} // namespace accretion
