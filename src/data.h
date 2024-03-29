#pragma once

#include <stddef.h>
#include <stdlib.h>

#include "nbody/planetary_system.h"
#include "polargrid.h"
#include "radialgrid.h"
#include "massflow_tracker.h"

class t_data
{
  public:
    t_data();
    ~t_data();

    enum t_polargrid_type {
	SIGMA,		      // surface density
	V_RADIAL,	      // radial velocity
	V_AZIMUTHAL,	      // azimuthal velocity
	ENERGY,		      // energy
	TEMPERATURE,	      // temperature
	PRESSURE,	      // pressure
	SOUNDSPEED,	      // soundspeed c_s (introduced by Fargo-ADSG)
	GAMMAEFF,	      // effective adiabatic indices
	GAMMA1,		      // first adiabatic indices
	ALPHA,			//local alpha
	MU,		      // mean molecular weights
	ACCEL_RADIAL,	      // gas accerleration in r direction
	ACCEL_AZIMUTHAL,      // gas acceleration in phi direction
	SG_ACCEL_RAD,	      // gas accerleration in r direction due to self-gravity, cell-centered
	SG_ACCEL_AZI,      // gas acceleration in phi direction due to self-gravity, cell-centered
	TOOMRE,		      // Toomre parameter Q
	ECCENTRICITY_X,	      // disk eccentricity vector component
	ECCENTRICITY_Y,	      // disk eccentricity vector component
	ALPHA_GRAV,	      // alpha grav
	ALPHA_GRAV_MEAN,      // alpha grav time average
	ALPHA_REYNOLDS,	      // alpha reynolds
	ALPHA_REYNOLDS_MEAN,  // alpha reynolds time average
	V_RADIAL0,	      // radial velocity at beginning
	V_AZIMUTHAL0,	      // azimuthal velocity at beginning
	SIGMA0,		      // surface density at beginning
	ENERGY0,	      // azimuthal velocity at beginning
	KAPPA,		      // opacity
	TAU_COOL,	      // cooling time
	QPLUS,		      // heating terms
	QMINUS,		      // cooling terms
	P_DIVV,		      // pdV work
	VISCOSITY,	      // (kinematic) viscosity
	SIGMA_ART_VISC, // (artificial) viscosity
	VISCOSITY_CORRECTION_FACTOR_R,
	VISCOSITY_CORRECTION_FACTOR_PHI,
	ART_VISCOSITY_CORRECTION_FACTOR_R,
	ART_VISCOSITY_CORRECTION_FACTOR_PHI,
	VISCOSITY_SIGMA_RP,
	VISCOSITY_SIGMA,
	TAU_R_R,	 // viscous stress tensor r,r component
	TAU_R_PHI,	 // viscous stress tensor r,phi component
	TAU_PHI_PHI,	 // viscous stress tensor phi,phi component
	DIV_V,		 // divergence of velocity used for viscosity
	T_REYNOLDS,	 // Reynold Stress tensor
	T_GRAVITATIONAL, // graviational Stress tensor
	POTENTIAL, // the gravitational potential felt by any fluid elements
		   // (introduced by Fargo-ADSG)
		   // 1)
	DENSITY_INT,  //
	Q_R,	      //
	Q_PHI,	      //
	TAU,	      // optical depth
	TAU_EFF,      // effective optical depth (in vertical direction)
	SCALE_HEIGHT, // scale height H
	ASPECTRATIO,  // aspectratio h = H / r
	VISIBILITY,   //
	TORQUE,
	RHO,
	MASSFLOW,
	ADVECTION_TORQUE,
	VISCOUS_TORQUE,
	GRAVITATIONAL_TORQUE_NOT_INTEGRATED,
	WORKER_SCALAR_ARRAY,
	GAS_DIFFUSION_COEFFICIENT,
	DRHO_DR,

	// number of t_polargrid_types
	N_POLARGRID_TYPES
    };

    enum t_radialgrid_type {
	LUMINOSITY_1D, //
	DISSIPATION_1D,
	TORQUE_1D,
	SIGMA_1D,

	// number of t_radialgrid_types
	N_RADIALGRID_TYPES
    };

    inline t_polargrid &operator[](t_polargrid_type polargrid_type)
    {
	return m_polargrids[polargrid_type];
    }
    inline t_radialgrid &operator[](t_radialgrid_type radialgrid_type)
    {
	return m_radialgrids[radialgrid_type];
    }

    void set_size(unsigned int global_n_radial,
		  unsigned int global_n_azimiuthal, unsigned int n_radial,
		  unsigned int n_azimuthal);
    inline t_planetary_system &get_planetary_system()
    {
	return m_planetary_system;
    }
	inline t_massflow_tracker &get_massflow_tracker(){
		return m_massflow_tracker;
	}

    void print_memory_usage(unsigned int n_radial, unsigned int n_azimuthal);

    double pdivv_total;
	inline ptrdiff_t get_n_radial(void) { return m_n_radial; }
	inline ptrdiff_t get_n_azimuthal(void) { return m_n_azimuthal; }
	inline ptrdiff_t get_global_n_radial(void) { return m_global_n_radial; }
    inline ptrdiff_t get_global_n_azimuthal(void)
    {
	return m_global_n_azimuthal;
	}

  private:
    ptrdiff_t m_n_radial;
    ptrdiff_t m_n_azimuthal;
    ptrdiff_t m_global_n_radial;
    ptrdiff_t m_global_n_azimuthal;
    t_polargrid m_polargrids[N_POLARGRID_TYPES];
    t_radialgrid m_radialgrids[N_RADIALGRID_TYPES];
    t_planetary_system m_planetary_system;
	t_massflow_tracker m_massflow_tracker;

};