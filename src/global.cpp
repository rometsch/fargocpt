/**
	\file global.c

	Declares all global variables.
*/

#include <mpi.h>
#include "types.h"
#include "radialarray.h"

/** number of this process, not an unsigned integer because MPI excepts it to be signed */
int CPU_Rank;

/** total number of processes, not an unsigned integer because MPI excepts it to be signed */
int CPU_Number;

/** is this process the master process */
int CPU_Master;

/** number of upper next CPU */
int CPU_Next;

/** number of lower next CPU */
int CPU_Prev;

/** number of upper most CPU */
int CPU_Highest;

/* ------------------------------------- */
/* Variables specific to fftw mesh split */
/* ------------------------------------- */
int CPU_Friend, CPU_NoFriend;
double *dens_friend;
double *SGP_buffft_Accr_friend, *SGP_buffft_Acct_friend;
double *ffttohydro_transfer, *ffttohydro_transfer_friend;

ptrdiff_t local_Nx;
ptrdiff_t local_i_start;
ptrdiff_t total_local_size;

ptrdiff_t local_i_start_friend;
ptrdiff_t local_Nx_friend;
ptrdiff_t total_local_size_friend;

ptrdiff_t local_Ny_after_transpose;
ptrdiff_t local_j_start_after_transpose;

ptrdiff_t transfer_size;
ptrdiff_t transfer_size_friend;

ptrdiff_t ifront;

int Zero_or_active_friend;
int hydro_totalsize, active_hydro_totalsize, active_hydro_totalsize_friend;
/* ------------------------------------- */

/** radial index of the inner most cell in global radial mesh of this process including ghost cells (!) */
unsigned int IMIN;

/** radial index of the outer most cell in global radial mesh of this process including ghost cells (!) */
unsigned int IMAX;

/** radial index of the inner most cell in global radial mesh of this process excluding ghost cells (!) */
unsigned int Zero_or_active;

/** radial index of the inner most cell in global radial mesh of this process excluding ghost cells (!) */
unsigned int Max_or_active;

unsigned int One_or_active;
unsigned int MaxMO_or_active;		/* MO: Minus One */
unsigned int GlobalNRadial;

int *GlobalNradialLocalSizes;			// Needed for MPI_Gatherv
int *GlobalNradialDisplacements;	// Needed for MPI_Gatherv

/** Rmed is the location of be the center of mass of the cell.
    Its definition is in fact : 0.5 * [ (4/3) \pi Rsup[i]^3 - (4/3) \pi Rinf[i]^3 ] / [ \pi Rsup[i]^2 - \pi Rinf[i]^2]
    or: one half of the elementary volume divided by the elementary surface.

    Note that this represents the position of the center of mass only for d\theta << \pi . If d\theta becomes large, then the center of mass can be out of the cell, closer to the center (it reached the center for d\theta=2\pi), but Rmed stays between Rinf and Rsup, and doesn't depend on d\theta in the code. (Aurélien Crida) */
t_radialarray Rmed;
t_radialarray &Rb = Rmed;

/** inverse of Rmed */
t_radialarray InvRmed;
t_radialarray &InvRb = InvRmed;

/** inner radius of a cell */
t_radialarray Rinf;
t_radialarray &Ra = Rinf;

/** inverse of Rinf */
t_radialarray InvRinf;
t_radialarray &InvRa = InvRinf;

/** outer radius of a cell */
t_radialarray Rsup;

/** surface of a cell */
t_radialarray Surf;

/** inverse of Surf */
t_radialarray InvSurf;

t_radialarray EnergyMed;
t_radialarray SigmaMed;

/** inverse of (RMed[i] - RMed[i - 1]) */
t_radialarray InvDiffRmed;

t_radialarray InvDiffRsup;
t_radialarray Radii;
t_radialarray GlobalRmed;
t_radialarray SigmaInf;
t_radialarray GLOBAL_bufarray;

double OmegaFrame, FrameAngle, PhysicalTime, PhysicalTimeInitial;
int TimeStep;
double HillRadius, mdcp, mdcp0, exces_mdcp;
int debug, OnlyInit;
int GotoNextOutput, ViscosityAlpha, RocheSmoothing, ThicknessSmoothingAtCell, ThicknessSmoothingAtPlanet;
int CentrifugalBalance, ExcludeHill, SloppyCFL;
MPI_Status global_MPI_Status;
t_polargrid *CellAbscissa, *CellOrdinate;

int OverridesOutputdir;

char *OUTPUTDIR;
char *PLANETCONFIG;

double DT;
unsigned int NINTERM;
unsigned int NTOT;
unsigned int N_iter;

double MASSTAPER;

unsigned int NRadial;
unsigned int NAzimuthal;
double RMIN;
double RMAX;

double ROCHESMOOTHING;
double ASPECTRATIO_REF;
double VISCOSITY;
double ALPHAVISCOSITY;
double SIGMASLOPE;
double OMEGAFRAME;
double IMPOSEDDISKDRIFT;
double FLARINGINDEX;
double ADIABATICINDEX; // Also used for polytropic energy equation
double POLYTROPIC_CONSTANT;

BoundaryFlow MassDelta;
bool Stockholm;
