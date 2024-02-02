#pragma once

#include "data.h"

namespace parameters
{

extern double aspectratio_ref;
extern int aspectratio_mode;
extern double constant_viscosity;
extern double viscous_alpha;
// Is viscous accretion turned on at at either boundary or Nbodies?
extern bool VISCOUS_ACCRETION;
extern double sigma_slope;
extern double IMPOSEDDISKDRIFT;
extern double flaring_index;
extern double ADIABATICINDEX;
extern double POLYTROPIC_CONSTANT;

extern bool CartesianParticles;
// Compute the gravity acceleration acting on the particles in Cartesian coordinates.
extern bool ParticleGravityCalcInCartesian;

extern double monitor_timestep;
extern unsigned int Nmonitor;
extern unsigned int Nsnap;
extern double quantities_radius_limit;
extern double disk_radius_mass_fraction;

// cells per scaleheight
extern double cps;

extern int ShockTube;
extern bool SpreadingRing;

// energy euations
extern bool Adiabatic;
extern bool Polytropic;
extern bool Locally_Isothermal;

extern bool variableGamma;
extern double hydrogenMassFraction;

/// Type of radial Grid
enum t_radial_grid {
    logarithmic_spacing,
    arithmetic_spacing,
    exponential_spacing,
    custom_spacing
};
extern t_radial_grid radial_grid_type;
extern const char *radial_grid_names[];
extern double exponential_cell_size_factor;

// energy equation
/// mean molecular mass
extern double MU;
/// minimum temperature allowed
extern double minimum_temperature;
/// maximum temperature allowed
extern double maximum_temperature;
/// enable viscous heating
extern bool heating_viscous_enabled;
/// viscous heating factor
extern double heating_viscous_factor;

extern bool heating_star_enabled;

/// local radiative cooling enabled
extern bool cooling_surface_enabled;
/// local radiative cooling factor
extern double surface_cooling_factor;
/// beta cooling enabled
extern bool cooling_beta_enabled;
/// beta cooling ramp up time
extern double cooling_beta_ramp_up;
/// beta cooling constant
extern double cooling_beta;
/// beta cooling to aspect ratio profile
extern bool cooling_beta_reference;
/// beta cooling to initial profile
extern bool cooling_beta_model;
/// beta cooling to minimum temperature
extern bool cooling_beta_floor;
/// local Scurve cooling enabled
extern bool cooling_scurve_enabled;
extern bool cooling_scurve_type;






// initialisation
enum t_initialize_condition {
    initialize_condition_profile,
    initialize_condition_profile_Nbody_centered,
    initialize_condition_read1D,
    initialize_condition_read2D
};

/// initialize condition for sigma
extern t_initialize_condition sigma_initialize_condition;
/// filename to read sigma profile from
extern std::string sigma_filename;
/// random seed
extern int random_seed;
/// randomize sigma?
extern bool sigma_randomize;
/// sigma random factor
extern double sigma_random_factor;
/// open simplex feature size
extern double sigma_feature_size;
/// sigma floor (in multiples of Sigma0)
extern double sigma_floor;
/// adjust sigma0 to have total diskmass = sigma_diskmass
extern bool sigma_adjust;
/// total diskmass
extern double sigma_diskmass;
/// Sigma0
extern double sigma0;

/// initiliaze condition for energy
extern t_initialize_condition energy_initialize_condition;
/// filename to read energy profile from
extern std::string energy_filename;

extern bool keep_mass_constant;


/// circumbinary ring
extern bool cbd_ring;
extern double cbd_ring_position;
extern double cbd_ring_width;
extern double cbd_decay_width;
extern double cbd_decay_exponent;
extern double cbd_ring_enhancement_factor;

extern double center_mass_density_correction_factor;

/// enable profile cutoff
extern bool profile_cutoff_outer;
/// profile cutoff point
extern double profile_cutoff_point_outer;
/// profile cutoff width
extern double profile_cutoff_width_outer;

/// enable profile cutoff
extern bool profile_cutoff_inner;
/// profile cutoff point
extern double profile_cutoff_point_inner;
/// profile cutoff width
extern double profile_cutoff_width_inner;

// type of artifical viscosity
enum t_artificial_viscosity {
    artificial_viscosity_none, // no artificial viscosity
    artificial_viscosity_TW,   // artificial viscosity based on Tscharnuter &
			       // Winkler, 1979
    artificial_viscosity_SN    // artificial viscosity based on Stone & Norman,
			       // 1991 (ZEUS)
};

/// type of artificial viscosity
extern t_artificial_viscosity artificial_viscosity;
/// factor/constant for artificial for viscosity. (CVNR)
extern double artificial_viscosity_factor;
/// enable/disable dissipation thru artificial viscosity
extern bool artificial_viscosity_dissipation;

extern double thickness_smoothing;
extern double thickness_smoothing_sg;
// evalutate the smoothing length at the planet location for compatibility with literature results
extern bool compatibility_smoothing_planetloc;
// do not smooth the gravity from the central star for compatibility with literature results
extern bool compatibility_no_star_smoothing;
// substract the azimuthal average of the density in the calculation of torques on the planets to be consistent with no self-gravity of the disk
extern bool correct_disk_selfgravity;


/// calculate disk (hydro)
extern bool calculate_disk;

// control centering of frame
extern unsigned int n_bodies_for_hydroframe_center;
extern unsigned int corotation_reference_body;
extern bool corotating;

extern int AlphaMode;
extern double alphaCold;
extern double alphaHot;

/// disk applies force onto star
/// and correction terms due to off barycenter system are applied
/// replaces feels_disk
extern bool disk_feedback;
extern bool accretion_feedback;
extern bool planet_orbit_disk_test;

extern bool fast_transport;
extern int hydro_integrator;

extern int indirect_term_mode;

/// factor for conversation from surface density to density
extern double density_factor;

/// factor for tau calculation
extern double tau_factor;
extern double tau_min;

extern bool v_azimuthal_with_quadropole_support;

/// factor for kappa calculation
extern double kappa_factor;

extern bool do_init_secondary_disk;

// self gravity
extern bool self_gravity;
// type of opacity
enum t_sg {
    sg_B,      // fourier trafo based on Baruteau PhD thesis with first order smoothing length
    sg_S,     // symmetric smoothing length based on masterthesis by Tobias Moldenhauer
    sg_BK       // exact solution for the kernel using bessel functions by Steven Rendon Restrepo
};
extern t_sg self_gravity_mode;
extern unsigned int self_gravity_steps_between_kernel_update;
extern double self_gravity_aspectratio_change_threshold;

extern bool body_force_from_potential;

// output
extern bool write_torques;

/// write disk quantities
extern bool bitwise_exact_restarting;
extern bool write_disk_quantities;
extern bool write_at_every_timestep;
extern bool write_lightcurves;
extern std::vector<double> lightcurves_radii;
extern bool write_massflow;

// runtime output
extern unsigned int log_after_steps;
extern double log_after_real_seconds;

// type of opacity
enum t_opacity {
    opacity_lin,      // opacity based on Lin & Papaloizou, 1985
    opacity_bell,     // opacity based on Bell & Lin, 1994
    opacity_const_op, // constant opacity
    opacity_simple    // eq. 30 from 'Thermohydrodynamics of Circumstellar Disks
		      // with High-Mass Planets
    // Gennaro D'Angelo1, Thomas Henning, and Willy Kley, 2003'
};

extern t_opacity opacity;

// For use of constant opacity
extern double kappa_const;

/// initialize pure keplerian
extern bool initialize_pure_keplerian;
extern bool initialize_vradial_zero;

/// boundary layer parameters
extern double radial_viscosity_factor;

extern double accretion_radius_fraction;
extern double klahr_smoothing_radius;

extern double visc_accret_massflow_test; // gas massflow test specific parameter

/// CFL Factor
extern double CFL;
extern double CFL_max_var;
extern double HEATING_COOLING_CFL_LIMIT;

/// (total) number of particles
extern unsigned int number_of_particles;
/// enable particle integration
extern bool integrate_particles;
/// particle radius
extern double particle_radius;
/// number of particle species
extern unsigned int particle_species_number;
/// factor of increase from particle species to next larger
extern double particle_radius_increase_factor;
/// particle eccentricity
extern double particle_eccentricity;
/// particle density
extern double particle_density;
/// particle slope
extern double particle_slope;
/// particle minimum radius
extern double particle_minimum_radius;
/// particle maximum radius
extern double particle_maximum_radius;
/// particle escape radius
extern double particle_minimum_escape_radius_sq;
extern double particle_maximum_escape_radius_sq;
/// particle gas drag
extern bool particle_gas_drag_enabled;
/// particle disk gravity
extern bool particle_disk_gravity_enabled;
/// particle dust diffusion
extern bool particle_dust_diffusion;
/// particle integrator
enum t_particle_integrator {
    integrator_adaptive,     // adaptive Cash-Karp integrator
    integrator_exponential_midpoint // exponential midpoint integrator
				     // presented in Mignone et al. 2019: A
				     // PARTICLE MODULE FOR THE PLUTO CODE: III
				     // - DUST
};
extern t_particle_integrator particle_integrator;

void read(const std::string &filename, t_data &data);
void summarize_parameters();
void write_grid_data_to_file();
void exitOnDeprecatedSetting(std::string setting_name, std::string reason,
			     std::string instruction);
} // namespace parameters
