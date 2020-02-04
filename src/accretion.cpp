/**
	\file accretion.cpp
	\author Thomas Rometsch <thomas.rometsch@uni-tuebingen.de>

	Handling accretion onto planets and stars.
*/

#include <tuple>
#include <cmath>

#include "constants.h"
#include "global.h"
#include "parameters.h"
#include "planet.h"

namespace accretion
{
	std::tuple<unsigned int,unsigned int> hill_radial_index(const double Rplanet,
										  const double RHill,
										  const unsigned int nr) {
		/* Calculate the indeces in radial direction where
		   the Hill sphere starts and stops */
		unsigned int i_min = 0;
		while ((Rsup[i_min] < Rplanet - RHill) && (i_min < nr)) {
			i_min++;
		}

		unsigned int i_max = nr > 0 ? nr -1 : 0;
		while ((Rinf[i_max] > Rplanet + RHill) && (i_max > 0)) {
			i_max--;
		}
		std::tuple<unsigned int, unsigned int> ids(i_min, i_max);
		return ids;
	}


	std::tuple<int, int> hill_azimuthal_index(const double angle,
											  const double Rplanet,
											  const double RHill,
											  const unsigned int ns) {
		/* Calculate the index in azimuthal direction
		   where the Hill sphere starts and stops */
		const int i_min = (int)((double)ns / 2.0 / PI * (angle - 2.0 * RHill / Rplanet));
		const int i_max = (int)((double)ns / 2.0 / PI * (angle + 2.0 * RHill / Rplanet));
		std::tuple<int, int> ids(i_min, i_max);
		return ids;
	}

	void update_planet(t_planet &planet,
					   const double dMplanet,
					   const double dPxPlanet,
					   const double dPyPlanet) {
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


	bool AccreteOntoSinglePlanet(t_data &data, t_planet &planet, double dt) {
		bool mass_changed = false;
		const int nr = data[t_data::DENSITY].Nrad;
		const int ns = data[t_data::DENSITY].Nsec;
		double* dens = data[t_data::DENSITY].Field;
		const double* cell_center_x = CellCenterX->Field;
		const double* cell_center_y = CellCenterY->Field;
		const double* vrad = data[t_data::V_RADIAL].Field;
		const double* vtheta = data[t_data::V_AZIMUTHAL].Field;

		if (planet.get_acc() > 1e-10) {
			// Hereafter : initialization of W. Kley's parameters
			// remove a ratio of facc = planet.get_acc() of the mass inside the
			// Hill sphere every free fall time at the Hill radius
			double facc = dt * planet.get_acc()
				* planet.get_omega() * sqrt(12.0)
				/2.0 / PI;
			const double facc1 = 1.0 / 3.0 * facc;
			const double facc2 = 2.0 / 3.0 * facc;
			const double frac1 = 0.75;
			const double frac2 = 0.45;


			// W. Kley's parameters initialization finished
			const double Xplanet = planet.get_x();
			const double Yplanet = planet.get_y();

			const double Rplanet = planet.get_r();
			const double RHill = planet.get_rhill();

			// calculate range of indeces to iterate over
			const auto [i_min, i_max] = hill_radial_index(Rplanet, RHill, nr);
			const double angle = planet.get_phi();
			const auto [j_min, j_max] = hill_azimuthal_index(angle, Rplanet, RHill, ns);


			double dMplanet = 0.0;
			double dPxPlanet = 0.0;
			double dPyPlanet = 0.0;

			for (int i = i_min; i <= i_max; i++) {
				for (int j = j_min; j <= j_max; j++) {
					// map azimuthal index to [0, ns]
					int jf = j;
					while (jf < 0)
						jf += ns;
					while (jf >= ns)
						jf -= ns;
					// calculate cell 1d index
					int l = jf + i * ns;
					int lip = l + ns;
					int ljp = l + 1;
					if (jf == ns - 1)
						ljp = i * ns;

					const double xc = cell_center_x[l];
					const double yc = cell_center_y[l];
					const double dx = Xplanet - xc;
					const double dy = Yplanet - yc;
					const double distance = sqrt(dx * dx + dy * dy);

					// interpolate velocities to cell centers
					const double vtcell = 0.5 * (vtheta[l] + vtheta[ljp])
						+ Rmed[i] * OmegaFrame;
					const double vrcell = 0.5 * (vrad[l] + vrad[lip]);
					// calculate cartesian velocities
					const double vxcell = (vrcell * xc - vtcell * yc) / Rmed[i];
					const double vycell = (vrcell * yc + vtcell * xc) / Rmed[i];

					double deltaM = 0.0;
					// handle accretion zone 1
					if (distance < frac1 * RHill) {
						deltaM = facc1 * dens[l] * Surf[i];
						if (i < One_or_active) {
							deltaM = 0.0;
						} else if (i >= Max_or_active) {
							deltaM = 0.0;
						} else {
							dens[l] *= (1.0 - facc1);
							dPxPlanet += deltaM * vxcell;
							dPyPlanet += deltaM * vycell;
							dMplanet += deltaM;
						}
					}
					// handle accretion zone 2
					if (distance < frac2 * RHill) {
						deltaM = facc2 * dens[l] * Surf[i];
						if (i < One_or_active) {
							deltaM = 0.0;
						} else if (i >= Max_or_active) {
							deltaM = 0.0;
						} else {
							dens[l] *= (1.0 - facc2);
							dPxPlanet += deltaM * vxcell;
							dPyPlanet += deltaM * vycell;
							dMplanet += deltaM;
						}
					}
				}
			}

			// MPI reduce
			double temp;
			MPI_Allreduce(&dMplanet, &temp, 1, MPI_DOUBLE, MPI_SUM,
						  MPI_COMM_WORLD);
			dMplanet = temp;
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



	void AccreteOntoPlanets(t_data &data, double dt)
	{
		bool masses_changed = false;

		auto &planetary_system = data.get_planetary_system();
		for (unsigned int k = 0; k < planetary_system.get_number_of_planets();
			 k++) {
			auto &planet = planetary_system.get_planet(k);
			const bool changed = AccreteOntoSinglePlanet(data, planet, dt);
			masses_changed = masses_changed || changed;
		}

		// update hydro center mass
		if (masses_changed) {
			planetary_system.update_global_hydro_frame_center_mass();
		}
	}

} // namespace accretion
