#ifndef PARTICLES_H
#define PARTICLES_H

#include "data.h"
#include "particle.h"

namespace particles
{

extern t_particle *particles;
extern unsigned int particles_size;
extern unsigned int local_number_of_particles;

void init(t_data &data);
void restart(unsigned int timestep);
// void calculate_accelerations_from_star_and_planets(t_data& data);
void calculate_accelerations_from_star_and_planets(
    double &ar, double &aphi, const double r, const double r_dot,
    const double phi, const double phi_dot, t_data &data);
void calculate_accelerations_from_star_and_planets_cart(double &ax, double &ay,
							const double x,
							const double y,
							t_data &data);
void calculate_derivitives_from_star_and_planets(double &grav_r_ddot,
						 double &grav_l_dot,
						 const double r,
						 const double phi,
						 t_data &data);

void update_velocities_from_gas_drag(t_data &data, double dt);
void update_velocities_from_gas_drag_cart(t_data &data, double dt);
void update_velocity_from_disk_gravity_cart_old(t_data &data, double dt);

void update_velocity_from_disk_gravity(const int n_radial_a_minus,
				       const int n_radial_a_plus,
				       const int n_azimuthal_b_minus,
				       const int n_azimuthal_b_plus,
				       const double r, const double phi,
				       const int particle_id, const double dt);
void update_velocity_from_disk_gravity_cart(
    const int n_radial_a_minus, const int n_radial_a_plus,
    const int n_azimuthal_b_minus, const int n_azimuthal_b_plus, const double r,
    const double phi, const int particle_id, const double dt);
void check_tstop(t_data &data);
void integrate(t_data &data, double dt);
void integrate_implicit(t_data &data, const double dt);
void integrate_explicit(t_data &data, const double dt);
void integrate_explicit_adaptive(t_data &data, const double dt);
void integrate_exponential_midpoint(t_data &data, const double dt);
void integrate_semiimplicit(t_data &data, const double dt);
void write(unsigned int timestep);
void move(void);
void rotate(double Omega, double dt);

} // namespace particles

#endif // PARTICLES_H
