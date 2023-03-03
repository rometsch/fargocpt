#include "split.h"
#include "LowTasks.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "parameters.h"
#include "selfgravity.h"
#include "time.h"
#include <fftw3-mpi.h>
#include <vector>

#ifndef NDEBUG
#undef FFTW_MEASURE
#define FFTW_MEASURE FFTW_ESTIMATE
#endif

#include <algorithm>
#include <vector>

void SplitDomain()
{
    /* Variables specific to standard mesh split */
    int remainder = 0;
    unsigned int size_low = 0, size_high = 0;
    /* Variables specific to fftw mesh split */
    int IMAX_friend;
    int one_if_odd;
    ptrdiff_t sg_split[5], hydro_split[5];

    GlobalNRadial = NRadial;

    /* Standard domain decomposition */
    if (!parameters::self_gravity) {
	logging::print_master(
	    LOG_DEBUG "SplitDomain: Doing standard domain decomposition.\n");

	size_low = NRadial / CPU_Number;
	size_high = size_low + 1;
	remainder = NRadial % CPU_Number;

	if (size_low < 2 * CPUOVERLAP) {
	    logging::print_master(
		LOG_ERROR
		"The number of processes is too large or the mesh is radially too narrow.\n");
	    PersonalExit(1);
	}

	if (CPU_Rank < remainder) {
	    IMIN = size_high * CPU_Rank;
	    IMAX = IMIN + size_high - 1;
	} else {
	    IMIN = size_high * remainder + (CPU_Rank - remainder) * size_low;
	    IMAX = IMIN + size_low - 1;
	}

	if (CPU_Rank > 0)
	    IMIN -= CPUOVERLAP;

	if (CPU_Rank < CPU_Number - 1)
	    IMAX += CPUOVERLAP;

	/* set local NRadial */
	NRadial = IMAX - IMIN + 1;

	Zero_no_ghost = (CPU_Rank == 0) ? 1 : 0;
	One_no_ghost_vr = (CPU_Rank == 0) ? 2 : 1;

	Max_no_ghost = NRadial - ((CPU_Rank == CPU_Number - 1) ? 1 : 0);
	MaxMo_no_ghost_vr = NRadial+1 - ((CPU_Rank == CPU_Number - 1) ? 2 : 1);

	Zero_or_active = (CPU_Rank == 0) ? 0 : CPUOVERLAP;
	radial_first_active = (CPU_Rank == 0) ? GHOSTCELLS_B : CPUOVERLAP;
	Max_or_active =
	    NRadial - ((CPU_Rank == CPU_Number - 1) ? 0 : CPUOVERLAP);
	radial_active_size =
	    NRadial -
	    ((CPU_Rank == CPU_Number - 1) ? GHOSTCELLS_B : CPUOVERLAP);

	CPU_Highest = CPU_Number - 1;

	/* Set next and previous CPU */
	if (CPU_Rank > 0)
	    CPU_Prev = CPU_Rank - 1;

	if (CPU_Rank < CPU_Highest)
	    CPU_Next = CPU_Rank + 1;
    } else {
	// Domain decomposition imposed by fftw_mpi
	logging::print_master(
	    LOG_DEBUG
	    "SplitDomain: Doing fftw-mesh implied domain decomposition.\n");

	logging::print_master(
	    LOG_INFO
	    "SplitDomain: Doing some FFT test measures to optimize FFT calculations. This may take a few seconds.\n");

	// calculate fft mesh sizes
	total_local_size = fftw_mpi_local_size_2d_transposed(
	    2 * GlobalNRadial, NAzimuthal, MPI_COMM_WORLD, &local_Nx,
	    &local_i_start, &local_Ny_after_transpose,
	    &local_j_start_after_transpose);

	if (CPU_Number % 2 == 0) {
	    if ((CPU_Rank == CPU_Number / 2 - 1) &&
		(local_i_start + local_Nx - 1 < GlobalNRadial - 1)) {
		logging::print(
		    LOG_ERROR
		    "ERROR: Bad choice of number of processes for current value of GLOBALNRAD.\n");
		PersonalExit(1);
	    }
	}

	if (local_Nx < 2 * CPUOVERLAP) {
	    logging::print(
		LOG_ERROR
		"ERROR: At least one of the processes has a ring number (%u) less than the total number of ghosts (%u): I must exit.\n",
		local_Nx, 2 * CPUOVERLAP);
	    PersonalExit(1);
	}

	one_if_odd = (CPU_Number % 2 == 0 ? 0 : 1);
	ifront = local_Nx / 2 - 1 + (local_Nx % 2 == 0 ? 0 : 1);

	/* Each CPU communicates with its 'CPU_Friend' the info from the hydro
	 * mesh to the fftw mesh, and vice versa */
	CPU_Friend = (CPU_Rank + (CPU_Number + one_if_odd) / 2) %
		     (CPU_Number + one_if_odd);

	if (CPU_Number % 2 == 1) {
	    CPU_NoFriend = (CPU_Number - 1) / 2;
	    CPU_Highest = CPU_NoFriend;
	} else {
	    CPU_NoFriend =
		CPU_Number; // does not exist, on purpose. Do not erase!!
	    CPU_Highest = CPU_Number - 1;
	}

	/* Each CPU communicates boundaries with CPU_Prev and CPU_Next, except
	 * CPU 0 and CPU_Highest */
	if (CPU_Rank < (CPU_Number - one_if_odd) / 2) {
	    if (CPU_Rank == 0) {
		CPU_Prev = CPU_Rank;
		CPU_Next = CPU_Friend;
	    } else {
		CPU_Prev = CPU_Friend - 1;
		CPU_Next = CPU_Friend;
	    }
	} else {
	    if (CPU_Rank == CPU_Highest) {
		if (CPU_Number % 2 == 0) {
		    CPU_Prev = CPU_Friend;
		} else {
		    CPU_Prev = 2 * CPU_Rank;
		}
		CPU_Next = CPU_Rank;
	    } else {
		CPU_Prev = CPU_Friend;
		CPU_Next = CPU_Friend + 1;
	    }
	}

	/* The hydro mesh dd can now be defined */
	if (CPU_Rank < (CPU_Number + one_if_odd) / 2) {
	    IMIN = local_i_start;
	    IMAX = local_i_start +
		   (local_Nx + (local_Nx % 2 == 0 ? 0 : 1)) / 2 - 1;

	    if ((CPU_Rank == CPU_Highest) && (IMAX >= GlobalNRadial))
		IMAX = GlobalNRadial - 1;

	    if ((CPU_Number % 2 == 0) && (CPU_Rank == CPU_Number / 2 - 1)) {
		transfer_size =
		    2 *
		    (GlobalNRadial - local_i_start - ifront - 1 + CPUOVERLAP) *
		    2 * (NAzimuthal / 2 + 1);
	    } else {
		transfer_size = 2 * (local_Nx - ifront - 1 + CPUOVERLAP) * 2 *
				(NAzimuthal / 2 + 1);
	    }

	    if (CPU_Rank != CPU_NoFriend) {
		ffttohydro_transfer =
		    (double *)malloc(sizeof(double) * transfer_size);
		if (ffttohydro_transfer == NULL) {
		    logging::print(
			LOG_ERROR
			"No enougth memory for ffttohydro_transfer tab in split.c ...\n");
		    PersonalExit(1);
		}
	    }

	    sg_split[0] = IMAX;
	    sg_split[1] = local_i_start;
	    sg_split[2] = local_Nx;
	    sg_split[3] = total_local_size;
	    sg_split[4] = transfer_size;

	    if (CPU_Rank != CPU_NoFriend)
		MPI_Send(sg_split, 5, MPI_AINT, CPU_Friend, 10, MPI_COMM_WORLD);
	} else {
	    MPI_Recv(hydro_split, 5, MPI_AINT, CPU_Friend, 10, MPI_COMM_WORLD,
		     &global_MPI_Status);
	    IMAX_friend = hydro_split[0];
	    local_i_start_friend = hydro_split[1];
	    local_Nx_friend = hydro_split[2];
	    total_local_size_friend = hydro_split[3];
	    transfer_size_friend = hydro_split[4];

	    SGP_buffft_Accr_friend =
		(double *)malloc(sizeof(double) * total_local_size_friend);
	    SGP_buffft_Acct_friend =
		(double *)malloc(sizeof(double) * total_local_size_friend);
	    if ((SGP_buffft_Accr_friend == NULL) ||
		(SGP_buffft_Acct_friend == NULL)) {
		logging::print(
		    LOG_ERROR
		    "No enougth memory for SGP_buffft_Acc_friend in split.c ...\n");
		PersonalExit(1);
	    }

	    ffttohydro_transfer_friend =
		(double *)malloc(sizeof(double) * transfer_size_friend);
	    if (ffttohydro_transfer_friend == NULL) {
		logging::print(
		    LOG_ERROR
		    "No enougth memory for ffttohydro_transfer_friend tab in split.c ...\n");
		PersonalExit(1);
	    }

	    IMIN = IMAX_friend + 1;
	    IMAX = local_i_start_friend + local_Nx_friend - 1;

	    if ((CPU_Rank == CPU_Highest) && (IMAX >= GlobalNRadial))
		IMAX = GlobalNRadial - 1;
	}

	if (CPU_Rank > 0)
	    IMIN -= CPUOVERLAP;

	if (CPU_Rank != CPU_Highest)
	    IMAX += CPUOVERLAP;

	NRadial = IMAX - IMIN + 1;

	if (NRadial < 2 * CPUOVERLAP) {
	    logging::print(LOG_ERROR "CPU %d: local_Nx = %d\n", CPU_Rank,
			   NRadial);
	    logging::print(
		LOG_ERROR
		"ERROR: One of the processes has a ring number less than the total number of ghosts: I must exit.\n");
	    PersonalExit(1);
	}

	Zero_no_ghost = (CPU_Rank == 0) ? 1 : 0;
	One_no_ghost_vr = (CPU_Rank == 0) ? 2 : 1;

	Max_no_ghost = NRadial - ((CPU_Rank == CPU_Number - 1) ? 0 : 1);
	MaxMo_no_ghost_vr = NRadial+1 - ((CPU_Rank == CPU_Number - 1) ? 1 : 2);

	Zero_or_active = CPUOVERLAP * (CPU_Rank > 0 ? 1 : 0);
	radial_first_active = 1 + (CPUOVERLAP - 1) * (CPU_Rank > 0 ? 1 : 0);
	Max_or_active =
	    NRadial - CPUOVERLAP * (CPU_Rank != CPU_Highest ? 1 : 0);
	radial_active_size =
	    NRadial - 1 - (CPUOVERLAP - 1) * (CPU_Rank != CPU_Highest ? 1 : 0);
	hydro_totalsize = NRadial * NAzimuthal;
	active_hydro_totalsize =
	    (Max_or_active - Zero_or_active + 1) * NAzimuthal;

	if ((CPU_Number % 2 == 0) || (CPU_Rank != CPU_NoFriend)) {
	    if (CPU_Rank >= (CPU_Number + one_if_odd) / 2) {
		MPI_Ssend(&active_hydro_totalsize, 1, MPI_INT, CPU_Friend, 20,
			  MPI_COMM_WORLD);
		MPI_Ssend(&Zero_or_active, 1, MPI_INT, CPU_Friend, 22,
			  MPI_COMM_WORLD);
	    } else {
		MPI_Recv(&active_hydro_totalsize_friend, 1, MPI_INT, CPU_Friend,
			 20, MPI_COMM_WORLD, &global_MPI_Status);
		MPI_Recv(&Zero_or_active_friend, 1, MPI_INT, CPU_Friend, 22,
			 MPI_COMM_WORLD, &global_MPI_Status);
		dens_friend = (double *)malloc(sizeof(double) *
					       active_hydro_totalsize_friend);
		if (dens_friend == NULL) {
		    logging::print(
			LOG_ERROR
			"No enougth memory for dens_friend in split.c ...\n");
		    PersonalExit(1);
		}
	    }
	}
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    /// Added by Lucas
    // Needed for MPI_Gatherv
    // Initializes Radial sizes of non vector arrays and required displacements
    if (CPU_Master) {
	RootNradialLocalSizes = new int[CPU_Number];
	RootNradialDisplacements = new int[CPU_Number];
	RootIMAX = new int[CPU_Number];
	RootIMIN = new int[CPU_Number];
	RootRanksOrdered = new int[CPU_Number];
    }

    const int local_array_start = IMIN + ((CPU_Rank == 0) ? 0 : CPUOVERLAP);
    const int local_array_end =
	IMAX - (CPU_Rank == CPU_Highest ? 0 : CPUOVERLAP);

    const int send_size = local_array_end - local_array_start + 1;

    MPI_Gather(&IMAX, 1, MPI_INT, RootIMAX, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&IMIN, 1, MPI_INT, RootIMIN, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&send_size, 1, MPI_INT, RootNradialLocalSizes, 1, MPI_INT, 0,
	       MPI_COMM_WORLD);

    MPI_Gather(&CPU_Next, 1, MPI_INT, RootRanksOrdered, 1, MPI_INT, 0,
	       MPI_COMM_WORLD);

    if (CPU_Master) {

	std::vector<int> tmp(CPU_Number);
	tmp[0] = 0;
	for (int i = 1; i < CPU_Number; ++i) {
	    tmp[i] = RootRanksOrdered[tmp[i - 1]];
	}
	for (int i = 0; i < CPU_Number; ++i) {
	    RootRanksOrdered[i] = tmp[i];
	}

	// Displacement for the first chunk of data - 0
	for (int i = 0; i < CPU_Number; i++) {
	    RootNradialDisplacements[i] =
		(i > 0) ? (RootNradialDisplacements[i - 1] +
			   RootNradialLocalSizes[i - 1])
			: 0;
	    RootIMIN[i] =
		RootNradialDisplacements[i] + ((i == 0) ? GHOSTCELLS_B : 0);
	    RootIMAX[i] = RootNradialDisplacements[i] +
			  RootNradialLocalSizes[i] - 1 -
			  (i == CPU_Highest ? 1 : 0);
	}
    }
    ///////////////////////////////////////////////////////////////////////////////////////
    /// End Added by Lucas

    /* print debugging */
    logging::print(LOG_DEBUG "SplitDomain: DomainSplit Information:\n");
    if (!parameters::self_gravity)
	logging::print_master(LOG_DEBUG "SplitDomain: %d = %d * %d + %d\n",
			      GlobalNRadial, CPU_Number, size_low, remainder);
    logging::print(LOG_DEBUG "SplitDomain: IMIN: %d\n", IMIN);
    logging::print(LOG_DEBUG "SplitDomain: IMAX: %d\n", IMAX);
    logging::print(LOG_DEBUG "SplitDomain: NRAD: %d\n", NRadial);
    logging::print(LOG_DEBUG "SplitDomain: Zero_or_active: %d\n",
		   Zero_or_active);
    logging::print(LOG_DEBUG "SplitDomain: One_or_active: %d\n",
		   radial_first_active);
    logging::print(LOG_DEBUG "SplitDomain: Max_or_active: %d\n", Max_or_active);
    logging::print(LOG_DEBUG "SplitDomain: MaxMO_or_active: %d\n",
		   radial_active_size);
    logging::print(LOG_DEBUG "SplitDomain: GLOB: %d\n", GlobalNRadial);
    logging::print(LOG_DEBUG "SplitDomain: CPUOVERLAP: %d\n", CPUOVERLAP);
    if (parameters::self_gravity) {
	logging::print(LOG_DEBUG "SplitDomain: LocalNx: %ld\n", local_Nx);
	logging::print(LOG_DEBUG "SplitDomain: LocalIStart: %ld\n",
		       local_i_start);
	logging::print(LOG_DEBUG "SplitDomain: total_local_size: %ld\n",
		       total_local_size);

	logging::print(LOG_DEBUG "Friend: %ld\n", CPU_Friend);
	logging::print(LOG_DEBUG "SplitDomain: LocalNx_friend: %ld\n",
		       local_Nx_friend);
	logging::print(LOG_DEBUG "SplitDomain: LocalIStart_friend: %ld\n",
		       local_i_start_friend);
	logging::print(LOG_DEBUG "SplitDomain: total_local_size_friend: %ld\n",
		       total_local_size_friend);
    }
}

void FreeSplitDomain()
{

    if (dens_friend != nullptr) {
	free(dens_friend);
    }

    if (CPU_Master) {
	delete[] RootNradialLocalSizes;
	delete[] RootNradialDisplacements;
	delete[] RootIMAX;
	delete[] RootIMIN;
	delete[] RootRanksOrdered;
    }
}
