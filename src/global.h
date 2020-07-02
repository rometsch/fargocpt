#ifndef GLOBAL_H
#define GLOBAL_H

#include "radialarray.h"
#include "types.h"
#include <mpi.h>

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
extern unsigned int Zero_or_active;
extern unsigned int Max_or_active;
extern unsigned int One_or_active;
extern unsigned int MaxMO_or_active; /* MO: Minus One */
extern unsigned int GlobalNRadial;

void sum_without_ghost_cells(double &accumulator, const double &addend,
			     const unsigned int &n_radial);

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
extern t_radialarray Radii;
extern t_radialarray GlobalRmed;
extern t_radialarray SigmaInf;
extern t_radialarray GLOBAL_bufarray;

extern t_radialarray EnergyMed;
extern t_radialarray SigmaMed;

extern double OmegaFrame, PhysicalTime, PhysicalTimeInitial, FrameAngle;
extern int TimeStep;
extern double HillRadius, mdcp, mdcp0, exces_mdcp;
extern double hydro_center_mass;

extern int debug;
extern int GotoNextOutput, ViscosityAlpha, CartesianParticles, 
    ParticlesInCartesian;
extern int CentrifugalBalance, SloppyCFL;
extern MPI_Status global_MPI_Status;
extern t_polargrid *CellCenterX, *CellCenterY;

extern char *OUTPUTDIR;
extern char *PLANETCONFIG;

extern double DT;
extern unsigned int NINTERM;
extern unsigned int NTOT;
extern unsigned int N_iter;

extern unsigned int NAzimuthal;
extern unsigned int NRadial;
extern double RMIN;
extern double RMAX;

extern double ROCHESMOOTHING;
extern double ASPECTRATIO_REF;
extern double VISCOSITY;
extern double ALPHAVISCOSITY;
extern double SIGMASLOPE;
extern double OMEGAFRAME;
extern double IMPOSEDDISKDRIFT;
extern double FLARINGINDEX;
extern double ADIABATICINDEX;
extern double POLYTROPIC_CONSTANT;

extern BoundaryFlow MassDelta;

#endif // GLOBAL_H
