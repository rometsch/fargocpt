#pragma once

#include "radialarray.h"
#include "types.h"
#include <mpi.h>
#include <csignal>
#include <vector>
#include "hydro_dt_logger.h"
#include "polargrid.h"

extern int CPU_Rank;
extern int CPU_Number;
extern int Thread_Number;
extern int CPU_Master;
extern int CPU_Next;
extern int CPU_Prev;
extern int CPU_Highest;
extern int CPU_Friend, CPU_NoFriend;
extern double *dens_friend;
extern double *SGP_buffft_Accr_friend, *SGP_buffft_Acct_friend;
extern double *ffttohydro_transfer, *ffttohydro_transfer_friend;

extern ptrdiff_t local_Nx;
extern ptrdiff_t local_i_start;
extern ptrdiff_t total_local_size;

extern ptrdiff_t local_i_start_friend;
extern ptrdiff_t local_Nx_friend;
extern ptrdiff_t total_local_size_friend;

extern ptrdiff_t local_Ny_after_transpose;
extern ptrdiff_t local_j_start_after_transpose;

extern ptrdiff_t transfer_size;
extern ptrdiff_t transfer_size_friend;

extern ptrdiff_t ifront;

extern int Zero_or_active_friend;

extern int hydro_totalsize, active_hydro_totalsize,
    active_hydro_totalsize_friend;
extern unsigned int IMIN;
extern unsigned int IMAX;
extern unsigned int Zero_or_active;
extern unsigned int Zero_no_ghost;
extern unsigned int One_no_ghost_vr;
extern unsigned int Max_or_active;
extern unsigned int Max_no_ghost;
extern unsigned int MaxMo_no_ghost_vr;

extern unsigned int radial_first_active;
extern unsigned int radial_active_size;
extern unsigned int GlobalNRadial;

extern int *RootNradialLocalSizes;    // Needed for MPI_Gatherv
extern int *RootNradialDisplacements; // Needed for MPI_Gatherv
extern int *RootIMAX;
extern int *RootIMIN;
extern int *RootRanksOrdered;

extern t_radialarray Rmed;
extern t_radialarray &Rb;
extern t_radialarray InvRmed;
extern t_radialarray &InvRb;
extern t_radialarray Rinf;
extern t_radialarray &Ra;
extern t_radialarray InvRinf;
extern t_radialarray &InvRa;

extern t_radialarray Rsup;
extern t_radialarray Surf;
extern t_radialarray InvSurf;
extern t_radialarray InvDiffRmed;
extern t_radialarray InvDiffRsup;
extern t_radialarray InvDiffRsupRb;
extern t_radialarray TwoDiffRaSq;
extern t_radialarray TwoDiffRbSq;
extern t_radialarray FourThirdInvRbInvdphiSq;
extern t_radialarray Radii;
extern t_radialarray GlobalRmed;
extern t_radialarray SigmaInf;
extern t_radialarray GLOBAL_bufarray;
extern t_radialarray GLOBAL_AxiSGAccr;

extern t_radialarray EnergyMed;
extern t_radialarray SigmaMed;

extern double hydro_center_mass;

extern int debug;
extern int StabilizeViscosity;
extern int flux_limiter_type;
extern int CentrifugalBalance, SloppyCFL;
extern MPI_Status global_MPI_Status;
extern t_polargrid *CellCenterX, *CellCenterY;

extern double dphi;
extern double invdphi;
extern double dt_parabolic_local;
extern hydro_dt_logger dt_logger;

extern unsigned int NAzimuthal;
extern unsigned int NRadial;
extern double RMIN;
extern double RMAX;

extern BoundaryFlow MassDelta;

extern volatile sig_atomic_t SIGTERM_RECEIVED, PRINT_SIG_INFO;

extern int ECC_GROWTH_MONITOR;
extern double ecc_old, peri_old;
extern double delta_ecc_source, delta_peri_source;
extern double delta_ecc_art_visc, delta_peri_art_visc;
extern double delta_ecc_visc, delta_peri_visc;
extern double delta_ecc_transport, delta_peri_transport;
extern double delta_ecc_damp, delta_peri_damp;

extern double binary_quadropole_moment;

extern std::vector<double> g_xpl;
extern std::vector<double> g_ypl;
extern std::vector<double> g_mpl;
extern std::vector<double> g_rpl;
extern std::vector<double> g_cubic_smoothing_radius;
