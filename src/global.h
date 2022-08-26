#ifndef GLOBAL_H
#define GLOBAL_H

#include "radialarray.h"
#include "types.h"
#include <mpi.h>
#include <csignal>
#include <vector>
#include "hydro_dt_logger.h"

extern int CPU_Rank;
extern int CPU_Number;
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
extern unsigned int Zero_no_ghost;
extern unsigned int Zero_or_active;
extern unsigned int Max_or_active;
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

extern double OmegaFrame, PhysicalTime, PhysicalTimeInitial, FrameAngle;

extern double hydro_center_mass;

extern int debug;
extern int GotoNextOutput, ViscosityAlpha, CartesianParticles,
	ParticlesInCartesian, StabilizeViscosity, StabilizeArtViscosity;
extern int flux_limiter_type;
extern int CentrifugalBalance, SloppyCFL;
extern MPI_Status global_MPI_Status;
extern t_polargrid *CellCenterX, *CellCenterY;

extern char *OUTPUTDIR;
extern char *PLANETCONFIG;

extern char *PRESCRIBED_BOUNDARY_OUTER_FILE;
extern int PRESCRIBED_TIME_SEGMENT_NUMBER;

extern double dphi;
extern double invdphi;
extern double DT;
extern double last_dt;
extern double dt_parabolic_local;
extern double hydro_dt;

extern int N_output;
extern unsigned int N_outer_loop;
extern unsigned int N_hydro_iter;

extern hydro_dt_logger dt_logger;

extern unsigned int NINTERM;
extern unsigned int NTOT;

extern unsigned int NAzimuthal;
extern unsigned int NRadial;
extern double RMIN;
extern double RMAX;

extern double quantities_radius_limit;

extern double ASPECTRATIO_REF;
extern int ASPECTRATIO_MODE;
extern int EXPLICIT_VISCOSITY;
extern double STS_NU;
extern double VISCOSITY;
extern double ALPHAVISCOSITY;
extern int VISCOUS_ACCRETION;
extern double SIGMASLOPE;
extern double OMEGAFRAME;
extern double IMPOSEDDISKDRIFT;
extern double FLARINGINDEX;

extern double ADIABATICINDEX;
extern double POLYTROPIC_CONSTANT;

extern BoundaryFlow MassDelta;

extern double dtemp;
extern int debug_outputs;
extern volatile sig_atomic_t SIGTERM_RECEIVED, PRINT_SIG_INFO;


#endif // GLOBAL_H
