
#include "parallel.h"
#include <unistd.h>
#include "global.h"
#include "commbound.h"
#include "selfgravity.h"
#include "split.h"
#include "LowTasks.h"
#include "logging.h"
#include "fld.h"
#include "options.h"
#include <iostream>

// copy paste from https://github.com/jeffhammond/HPCInfo/blob/master/mpi/with-threads/mpi-openmp.c
#define MPI_THREAD_STRING(level)  \
		( level==MPI_THREAD_SERIALIZED ? "THREAD_SERIALIZED" : \
				( level==MPI_THREAD_MULTIPLE ? "THREAD_MULTIPLE" : \
						( level==MPI_THREAD_FUNNELED ? "THREAD_FUNNELED" : \
								( level==MPI_THREAD_SINGLE ? "THREAD_SINGLE" : "THIS_IS_IMPOSSIBLE" ) ) ) )

// TODO: NRad darf nicht größer sein als MAX1D
// copy operator for t_polargrid
// write all polargrids on error

void init_parallel(int argc, char *argv[]) {
	int CPU_NameLength;
    char CPU_Name[MPI_MAX_PROCESSOR_NAME + 1];

	int provided;
	int requested = MPI_THREAD_FUNNELED;
	MPI_Init_thread(&argc, &argv, requested, &provided);
	if(provided < requested) {
		die("MPI_Init_thread provided %s when %s was requested.  Exiting. \n",
					   MPI_THREAD_STRING(provided), MPI_THREAD_STRING(requested));
	}

    // initialize MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &CPU_Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &CPU_Number);
    selfgravity::mpi_init();

    // are we master CPU?
    CPU_Master = (CPU_Rank == 0 ? 1 : 0);

	pid_t pid = getpid();
    if (CPU_Master) {
		if (options::pidfile != "") {
			std::ofstream pidfile;
			pidfile.open(options::pidfile);
			pidfile << pid;
			pidfile.close();
		} else {
			printf("%ld\n",(long) pid);
		}
    }

	// make sure the pid of the master is printed before additional info
	MPI_Barrier(MPI_COMM_WORLD);

#ifdef _OPENMP
	#pragma omp parallel
	{
		// print information on which processor we're running
		MPI_Get_processor_name(CPU_Name, &CPU_NameLength);
		CPU_Name[CPU_NameLength] = '\0';
		int omp_id  = omp_get_thread_num();
		int omp_num = omp_get_num_threads();
		logging::print(LOG_INFO "MPI rank # %2d OpenMP thread # %2d of %2d on %s\n", CPU_Rank, omp_id, omp_num, CPU_Name);
		fflush(stdout);
		if (omp_id == 0) {
			Thread_Number = omp_num;
		}
	}
#else
	Thread_Number = 1;
	printf("MPI rank # %2d \n", CPU_Rank);
	fflush(stdout);
#endif


#ifndef _OPENMP
	// print information on which processor we're running
    MPI_Get_processor_name(CPU_Name, &CPU_NameLength);
    CPU_Name[CPU_NameLength] = '\0';
    logging::print(LOG_INFO "fargo: running on %s\n", CPU_Name);
#endif

	MPI_Barrier(MPI_COMM_WORLD);

}


void finalize_parallel() {
	DeallocateBoundaryCommunicationBuffers();
    selfgravity::mpi_finalize();
    FreeSplitDomain();
	MPI_Finalize();
	fld::finalize();
}