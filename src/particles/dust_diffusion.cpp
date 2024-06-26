#include "dust_diffusion.h"
#include "../find_cell_id.h"
#include "../global.h"
#include "../random/random_wrapper.h"
#include "../Theo.h"
#include "../parameters.h"
#include <filesystem>

namespace dust_diffusion
{

void init(t_data &data)
{
    compute_gas_diffusion_coefficient(data);
    compute_gas_density_radial_derivative(data);
}

/* Apply dust diffusion by modelling it as a Brownian motion.

Main reference is Charnoz et al. 2011 (DOI:10.1088/0004-637X/737/1/33).
See also Youdin & Lithwick 2007 (DOI:10.1016/j.icarus.2007.07.012).

Each timestep each particle gets a kick according to a probability and the
size of the kick is drawn from a normal distribution.

Kicks are only applied in the radial direction.
In the azimuthal direction, Keplerian shear dominates.
*/
void diffuse_dust(t_data &data, std::vector<t_particle> &particles,
		  const double dt, const unsigned int N_particles)
{
    if (parameters::calculate_disk) {
        compute_gas_diffusion_coefficient(data);
        compute_gas_density_radial_derivative(data);
    }

    #pragma omp parallel for
    for (unsigned int i = 0; i < N_particles; i++) {
	auto &particle = particles[i];
    const double deltar = kick_length(particle, data, dt);
    const double rold = particle.r;
    const double rnew = rold + deltar;
	particle.r = rnew;
    // correct azimuthal motion to avoid oscillations
    // for a particle on a circular orbit, v_phi^new = vK(r^new) instead of vK(r^old)
    particle.phi_dot *= std::pow(rold/rnew, 1.5); 
    }
}

/* Compute the Schmidt number from the Stokes number.

Use Youdin & Lithwick 2007 Eq. 37
Sc = (1 + St^2)^2/(1 + 4 St^2)
*/
static double compute_Schmidt_number(const double St)
{
    const double St2 = std::pow(St, 2);
    return std::pow(1 + St2, 2) / (1 + 4 * St2);
}

/* Kick a single dust particle according to Brownian motion according to
Charnoz et al. 2011 Eq. 21 (right now still 17 which ignores the varying dust
diffusion coeff!)

Use values from the grid cell the particle is currently in without
interpolation.
 */
double kick_length(t_particle &particle, t_data &data, const double dt)
{
    const double r = particle.get_distance_to_star();
    const double phi = particle.get_angle();
    const double St = particle.stokes;
    const double Sc = compute_Schmidt_number(St);

    const unsigned int n_rad = get_rinf_id(r);
    const unsigned int n_az = get_inf_azimuthal_id(phi);

    const double Dg = data[t_data::GAS_DIFFUSION_COEFFICIENT](n_rad, n_az);
    const double Dd = Dg / Sc;

    const double rho = data[t_data::RHO](n_rad, n_az);
    const double drho_dr = data[t_data::DRHO_DR](n_rad, n_az);

    const double OmegaK = calculate_omega_kepler(r);
    // multiply by an additional dt * OmegaK to match 1D radial dust diffusion eq.
    const double mean = Dd / rho * drho_dr * dt * dt * OmegaK;
    const double sigma = std::sqrt(2 * Dd * dt);

    const double snv = fargo_random::get_std_normal();
    // Correct for missing diffusion in direction tangential to r
    const double corr_2d = r * ( std::sqrt(1 + std::pow(sigma*snv/r,2)) - 1 );
    const double deltar = mean + snv * sigma + corr_2d;



    return deltar;
}

/* Gas diffusion coefficient Dg = alpha * cs * H
Charnoz et al. 2011
*/
void compute_gas_diffusion_coefficient(t_data &data)
{
    auto &Dg = data[t_data::GAS_DIFFUSION_COEFFICIENT];

    const unsigned int N_rad_max = Dg.get_max_radial();
    const unsigned int N_az_max = Dg.get_max_azimuthal();

    #pragma omp parallel for collapse (2)
    for (unsigned int n_rad = 0; n_rad <= N_rad_max; ++n_rad) {
	for (unsigned int n_az = 0; n_az <= N_az_max; ++n_az) {
	    const double alpha = parameters::viscous_alpha;
	    const double cs = data[t_data::SOUNDSPEED](n_rad, n_az);
	    const double h = data[t_data::ASPECTRATIO](n_rad, n_az);
	    const double r = Rmed(n_rad);
	    Dg(n_rad, n_az) = alpha * cs * h * r;
	}
    }
}

/* Compute the central derivate assuming a constantly spaced grid */
static double radial_central_derivative(t_polargrid &val, unsigned int n_rad,
				 unsigned int n_az)
{
    const double val_plus = val(n_rad + 1, n_az);
    const double val_minus = val(n_rad - 1, n_az);
    const double r_plus = Rmed[n_rad + 1];
    const double r_minus = Rmed[n_rad - 1];
    return (val_plus - val_minus) / (r_plus - r_minus);
}

/* Compute the radial derivative of the mass volume density.
Use central difference and exclude the ghost cells.
*/
void compute_gas_density_radial_derivative(t_data &data)
{
    auto &rho = data[t_data::RHO];
    auto &deriv = data[t_data::DRHO_DR];

    const unsigned int N_rad_max = deriv.get_max_radial();
    const unsigned int N_az_max = deriv.get_max_azimuthal();

    #pragma omp parallel for collapse (2)
    for (unsigned int n_rad = 1; n_rad <= N_rad_max - 1; ++n_rad) {
	for (unsigned int n_az = 0; n_az <= N_az_max; ++n_az) {
	    deriv(n_rad, n_az) = radial_central_derivative(rho, n_rad, n_az);
	}
    }
}

} // namespace dust_diffusion