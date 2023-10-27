
#include "data.h"
#include "logging.h"
#include "output.h"
#include "quantities.h"
#include "units.h"
#include <mpi.h>

/// constructor
t_data::t_data()
{
    // -- polar grids --
    m_polargrids[SIGMA].set_scalar(true);
    m_polargrids[SIGMA].set_name("Sigma");
    m_polargrids[SIGMA].set_unit(units::surface_density);

    m_polargrids[V_RADIAL].set_vector(true);
    m_polargrids[V_RADIAL].set_name("vrad");
    m_polargrids[V_RADIAL].set_unit(units::velocity);

    m_polargrids[V_AZIMUTHAL].set_vector(false);
    m_polargrids[V_AZIMUTHAL].set_name("vtheta");
    m_polargrids[V_AZIMUTHAL].set_unit(units::velocity);

    m_polargrids[ENERGY].set_scalar(true);
    m_polargrids[ENERGY].set_name("energy");
    m_polargrids[ENERGY].set_unit(units::energy_density);

    m_polargrids[TEMPERATURE].set_scalar(true);
    m_polargrids[TEMPERATURE].set_name("Temperature");
    m_polargrids[TEMPERATURE].set_unit(units::temperature);

    m_polargrids[PRESSURE].set_scalar(true);
    m_polargrids[PRESSURE].set_name("pressure");
    m_polargrids[PRESSURE].set_unit(units::pressure);

    m_polargrids[GAMMAEFF].set_scalar(true);
    m_polargrids[GAMMAEFF].set_name("gammaeff");

    m_polargrids[GAMMA1].set_scalar(true);
    m_polargrids[GAMMA1].set_name("gamma1");

    m_polargrids[ALPHA].set_scalar(true);
    m_polargrids[ALPHA].set_name("alpha");
	m_polargrids[ALPHA].set_do_before_write(&quantities::fill_alpha_array);

    m_polargrids[MU].set_scalar(true);
    m_polargrids[MU].set_name("mu");

    m_polargrids[ACCEL_RADIAL].set_vector(true);
    m_polargrids[ACCEL_RADIAL].set_name("a_rad");
    m_polargrids[ACCEL_RADIAL].set_unit(units::acceleration);

    m_polargrids[ACCEL_AZIMUTHAL].set_vector(false);
    m_polargrids[ACCEL_AZIMUTHAL].set_name("a_azi");
    m_polargrids[ACCEL_AZIMUTHAL].set_unit(units::acceleration);

    // sg accel is cell-centered
    m_polargrids[SG_ACCEL_RAD].set_vector(false);
    m_polargrids[SG_ACCEL_RAD].set_name("a_sg_rad");
    m_polargrids[SG_ACCEL_RAD].set_unit(units::acceleration);

    // sg accel is cell-centered
    m_polargrids[SG_ACCEL_AZI].set_vector(false);
    m_polargrids[SG_ACCEL_AZI].set_name("a_sg_azi");
    m_polargrids[SG_ACCEL_AZI].set_unit(units::acceleration);

    m_polargrids[TOOMRE].set_scalar(true);
    m_polargrids[TOOMRE].set_name("Toomre");
    m_polargrids[TOOMRE].set_do_before_write(&quantities::calculate_toomre);


    m_polargrids[ECCENTRICITY_X].set_scalar(true);
    m_polargrids[ECCENTRICITY_X].set_name("Eccentricity_X");
    m_polargrids[ECCENTRICITY_X].set_do_before_write(
	&quantities::calculate_disk_ecc_vector);

    m_polargrids[ECCENTRICITY_Y].set_scalar(true);
    m_polargrids[ECCENTRICITY_Y].set_name("Eccentricity_Y");
    /// Always computed together with ECCENTRICITY_X

    m_polargrids[ALPHA_GRAV].set_scalar(true);
    m_polargrids[ALPHA_GRAV].set_name("alpha_grav");
    m_polargrids[ALPHA_GRAV].set_do_before_write(
	&quantities::calculate_alpha_grav);

    m_polargrids[ALPHA_REYNOLDS].set_scalar(true);
    m_polargrids[ALPHA_REYNOLDS].set_name("alpha_reynolds");
    m_polargrids[ALPHA_REYNOLDS].set_do_before_write(
	&quantities::calculate_alpha_reynolds);

    m_polargrids[ALPHA_GRAV_MEAN].set_scalar(true);
    m_polargrids[ALPHA_GRAV_MEAN].set_name("alpha_grav_mean");
    m_polargrids[ALPHA_GRAV_MEAN].set_clear_after_write(true);

    m_polargrids[ALPHA_REYNOLDS_MEAN].set_scalar(true);
    m_polargrids[ALPHA_REYNOLDS_MEAN].set_name("alpha_reynolds_mean");
    m_polargrids[ALPHA_REYNOLDS_MEAN].set_clear_after_write(true);

    m_polargrids[V_RADIAL0].set_vector(true);
    m_polargrids[V_RADIAL0].set_name("vrad0");
    m_polargrids[V_RADIAL0].set_unit(units::velocity);

    m_polargrids[V_AZIMUTHAL0].set_vector(
	m_polargrids[V_AZIMUTHAL].is_vector());
    m_polargrids[V_AZIMUTHAL0].set_name("vtheta0");
    m_polargrids[V_AZIMUTHAL0].set_unit(units::velocity);

    m_polargrids[SIGMA0].set_vector(false);
    m_polargrids[SIGMA0].set_name("dens0");
    m_polargrids[SIGMA0].set_unit(units::surface_density);

    m_polargrids[ENERGY0].set_vector(false);
    m_polargrids[ENERGY0].set_name("energy0");
    m_polargrids[ENERGY0].set_unit(units::energy_density);

    m_polargrids[SOUNDSPEED].set_scalar(true);
    m_polargrids[SOUNDSPEED].set_name("soundspeed");
    m_polargrids[SOUNDSPEED].set_unit(units::velocity);

    m_polargrids[KAPPA].set_scalar(true);
    m_polargrids[KAPPA].set_name("kappa");
    m_polargrids[KAPPA].set_unit(units::opacity);

    m_polargrids[TAU_COOL].set_scalar(true);
    m_polargrids[TAU_COOL].set_name("tau_cool");
    m_polargrids[TAU_COOL].set_unit(units::time);

    m_polargrids[QPLUS].set_scalar(true);
    m_polargrids[QPLUS].set_name("Qplus");
    m_polargrids[QPLUS].set_unit(units::energy_flux);

    m_polargrids[QMINUS].set_scalar(true);
    m_polargrids[QMINUS].set_name("Qminus");
    m_polargrids[QMINUS].set_unit(units::energy_flux);

    m_polargrids[P_DIVV].set_scalar(true);
    m_polargrids[P_DIVV].set_name("P_divV");
    m_polargrids[P_DIVV].set_unit(units::energy_flux);

    m_polargrids[VISCOSITY].set_scalar(true);
    m_polargrids[VISCOSITY].set_name("viscosity");
    m_polargrids[VISCOSITY].set_unit(units::kinematic_viscosity);

	m_polargrids[SIGMA_ART_VISC].set_scalar(true);
	m_polargrids[SIGMA_ART_VISC].set_name("artificial_viscosity");

    m_polargrids[VISCOSITY_CORRECTION_FACTOR_R].set_scalar(true);
    m_polargrids[VISCOSITY_CORRECTION_FACTOR_R].set_name(
	"viscosity_correction_factor_r");

    m_polargrids[VISCOSITY_CORRECTION_FACTOR_PHI].set_scalar(true);
    m_polargrids[VISCOSITY_CORRECTION_FACTOR_PHI].set_name(
	"viscosity_correction_factor_phi");

	m_polargrids[ART_VISCOSITY_CORRECTION_FACTOR_R].set_scalar(true);
	m_polargrids[ART_VISCOSITY_CORRECTION_FACTOR_R].set_name(
	"artificial_viscosity_correction_factor_r");

	m_polargrids[ART_VISCOSITY_CORRECTION_FACTOR_PHI].set_scalar(true);
	m_polargrids[ART_VISCOSITY_CORRECTION_FACTOR_PHI].set_name(
	"artificial_viscosity_correction_factor_phi");

    m_polargrids[VISCOSITY_SIGMA_RP].set_vector(true);
    m_polargrids[VISCOSITY_SIGMA_RP].set_name(
	"viscosity_correction_factor_phi_helper_rp");

    m_polargrids[VISCOSITY_SIGMA].set_scalar(true);
    m_polargrids[VISCOSITY_SIGMA].set_name(
	"viscosity_correction_factor_phi_helper_pp");

    m_polargrids[ADVECTION_TORQUE].set_scalar(true);
    m_polargrids[ADVECTION_TORQUE].set_name("ADVECTION_TORQUE");
    m_polargrids[ADVECTION_TORQUE].set_unit(units::torque);
    m_polargrids[ADVECTION_TORQUE].set_do_before_write(
	&quantities::calculate_advection_torque);
	m_polargrids[ADVECTION_TORQUE].set_clear_after_write(true);
    m_polargrids[ADVECTION_TORQUE].set_integrate_azimuthally_for_1D_write(true);

    m_polargrids[VISCOUS_TORQUE].set_scalar(true);
    m_polargrids[VISCOUS_TORQUE].set_name("VISCOUS_TORQUE");
    m_polargrids[VISCOUS_TORQUE].set_unit(units::torque);
    m_polargrids[VISCOUS_TORQUE].set_do_before_write(
	&quantities::calculate_viscous_torque);
	m_polargrids[VISCOUS_TORQUE].set_clear_after_write(true);
    m_polargrids[VISCOUS_TORQUE].set_integrate_azimuthally_for_1D_write(true);

    m_polargrids[GRAVITATIONAL_TORQUE_NOT_INTEGRATED].set_scalar(true);
    m_polargrids[GRAVITATIONAL_TORQUE_NOT_INTEGRATED].set_name(
	"GRAVITATIONAL_TORQUE_NOT_INTEGRATED");
    m_polargrids[GRAVITATIONAL_TORQUE_NOT_INTEGRATED].set_unit(units::torque);
    m_polargrids[GRAVITATIONAL_TORQUE_NOT_INTEGRATED].set_do_before_write(
	&quantities::calculate_gravitational_torque);
	m_polargrids[GRAVITATIONAL_TORQUE_NOT_INTEGRATED].set_clear_after_write(
	true);
    m_polargrids[GRAVITATIONAL_TORQUE_NOT_INTEGRATED]
	.set_integrate_azimuthally_for_1D_write(true);

    m_polargrids[WORKER_SCALAR_ARRAY].set_scalar(true);
    m_polargrids[WORKER_SCALAR_ARRAY].set_name("WORKER_SCALAR_ARRAY");

    m_polargrids[GAS_DIFFUSION_COEFFICIENT].set_scalar(true);
    m_polargrids[GAS_DIFFUSION_COEFFICIENT].set_name("GAS_DIFFUSION_COEFFICIENT");

    m_polargrids[DRHO_DR].set_scalar(true);
    m_polargrids[DRHO_DR].set_name("DRHO_DR");

    // tau_r_r is cell centered
    m_polargrids[TAU_R_R].set_scalar(true);
    m_polargrids[TAU_R_R].set_name("tau_r_r");

    // tau_phi_phi is cell centered
    m_polargrids[TAU_PHI_PHI].set_scalar(true);
    m_polargrids[TAU_PHI_PHI].set_name("tau_phi_phi");

    // tau_r_phi is on the edges
    m_polargrids[TAU_R_PHI].set_vector(true);
    m_polargrids[TAU_R_PHI].set_name("tau_r_phi");

    m_polargrids[DIV_V].set_scalar(true);
    m_polargrids[DIV_V].set_name("div_v");

    m_polargrids[T_GRAVITATIONAL].set_scalar(true);
    m_polargrids[T_GRAVITATIONAL].set_name("T_gravitational");
    m_polargrids[T_GRAVITATIONAL].set_unit(units::stress);

    m_polargrids[T_REYNOLDS].set_scalar(true);
    m_polargrids[T_REYNOLDS].set_name("T_Reynolds");
    m_polargrids[T_REYNOLDS].set_unit(units::stress);

    m_polargrids[POTENTIAL].set_scalar(true);
    m_polargrids[POTENTIAL].set_name("potential");
    m_polargrids[POTENTIAL].set_unit(units::potential);

    m_polargrids[Q_PHI].set_vector(false);
    m_polargrids[Q_PHI].set_name("q_phi");

    m_polargrids[Q_R].set_vector(false);
    m_polargrids[Q_R].set_name("q_r");

    m_polargrids[DENSITY_INT].set_scalar(true);
    m_polargrids[DENSITY_INT].set_name("dens_int");
    m_polargrids[DENSITY_INT].set_unit(units::surface_density);

    m_polargrids[TAU].set_scalar(true);
    m_polargrids[TAU].set_name("tau");

    m_polargrids[TAU_EFF].set_scalar(true);
    m_polargrids[TAU_EFF].set_name("tau_eff");

    m_polargrids[SCALE_HEIGHT].set_scalar(true);
    m_polargrids[SCALE_HEIGHT].set_name("scale_height");

    m_polargrids[ASPECTRATIO].set_scalar(true);
    m_polargrids[ASPECTRATIO].set_name("aspectratio");
    m_polargrids[ASPECTRATIO].set_do_before_write(
	&quantities::compute_aspectratio);

    m_polargrids[VISIBILITY].set_scalar(true);
    m_polargrids[VISIBILITY].set_name("visiblity");

    m_polargrids[TORQUE].set_scalar(true);
    m_polargrids[TORQUE].set_name("torque");
    m_polargrids[TORQUE].set_unit(units::torque);

    m_polargrids[RHO].set_scalar(true);
    m_polargrids[RHO].set_name("rho");
    m_polargrids[RHO].set_unit(units::density);

    m_polargrids[MASSFLOW].set_scalar(false);
    m_polargrids[MASSFLOW].set_name("MassFlow");
    m_polargrids[MASSFLOW].set_unit(units::mass_accretion_rate);
    m_polargrids[MASSFLOW].set_do_before_write(&quantities::calculate_massflow);
    m_polargrids[MASSFLOW].set_clear_after_write(true);
    m_polargrids[MASSFLOW].set_integrate_azimuthally_for_1D_write(true);

    // -- radial grids --
    m_radialgrids[LUMINOSITY_1D].set_scalar(true);
    m_radialgrids[LUMINOSITY_1D].set_name("1D_Luminosity");
    m_radialgrids[LUMINOSITY_1D].set_unit(units::power);
    m_radialgrids[LUMINOSITY_1D].set_do_before_write(
	&quantities::calculate_radial_luminosity);

    m_radialgrids[DISSIPATION_1D].set_scalar(true);
    m_radialgrids[DISSIPATION_1D].set_name("1D_Dissipation");
    m_radialgrids[DISSIPATION_1D].set_unit(units::power);
    m_radialgrids[DISSIPATION_1D].set_do_before_write(
	&quantities::calculate_radial_dissipation);

    m_radialgrids[TORQUE_1D].set_scalar(true);
    m_radialgrids[TORQUE_1D].set_name("1D_torque");
    m_radialgrids[TORQUE_1D].set_unit(units::torque);

    m_radialgrids[SIGMA_1D].set_scalar(true);
    m_radialgrids[SIGMA_1D].set_name("Sigma");
    m_radialgrids[SIGMA_1D].set_unit(units::surface_density);

    pdivv_total = 0.0;
}

/// destructor
t_data::~t_data() {}

void t_data::set_size(unsigned int global_n_radial,
		      unsigned int global_n_azimiuthal, unsigned int n_radial,
		      unsigned int n_azimuthal)
{
    m_global_n_radial = global_n_radial;
    m_global_n_azimuthal = global_n_azimiuthal;
    m_n_radial = n_radial;
    m_n_azimuthal = n_azimuthal;

    for (unsigned int i = 0; i < N_POLARGRID_TYPES; ++i) {
	m_polargrids[i].set_size(m_n_radial, m_n_azimuthal);
    }

    for (unsigned int i = 0; i < N_RADIALGRID_TYPES; ++i) {
	m_radialgrids[i].set_size(m_n_radial);
    }
}

void t_data::print_memory_usage(unsigned int n_radial, unsigned int n_azimuthal)
{
    double local_memory_usage = 0;
    double global_memory_usage = 0;

    for (unsigned int i = 0; i < N_POLARGRID_TYPES; ++i) {
	local_memory_usage +=
	    m_polargrids[i].get_memory_usage(n_radial, n_azimuthal);
    }

    for (unsigned int i = 0; i < N_RADIALGRID_TYPES; ++i) {
	local_memory_usage += m_radialgrids[i].get_memory_usage(n_radial);
    }

    MPI_Allreduce(&local_memory_usage, &global_memory_usage, 1, MPI_DOUBLE,
		  MPI_SUM, MPI_COMM_WORLD);

    logging::print(
	LOG_INFO
	"Need about %.0lf bytes = %.2lf KB = %.2lf MB = %.2lf GB of memory on this process.\n",
	local_memory_usage, local_memory_usage / 1024.0,
	local_memory_usage / 1048576.0, local_memory_usage / 1073741824.0);
    logging::print_master(
	LOG_INFO
	"Need about %.0lf bytes = %.2lf KB = %.2lf MB = %.2lf GB of memory in total.\n",
	global_memory_usage, global_memory_usage / 1024.0,
	global_memory_usage / 1048576.0, global_memory_usage / 1073741824.0);
}
