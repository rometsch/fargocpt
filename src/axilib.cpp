#include "axilib.h"
#include "LowTasks.h"
#include "global.h"
#include "logging.h"

/**
	Make an 1D profile of a 2D grid by calculating mean values over
   azimuthal cells

	\param grid 2D grid
	\param array 1D array
*/
void mpi_make1Dprofile(double *grid, double *array)
{
    MPI_Request req;
    unsigned int nRadial, nAzimuthal, cell;

    double *tempArray;

    /* allocate memory for temporary array */
    tempArray = (double *)malloc(sizeof(double) * NRadial);

    if (tempArray == NULL) {
	logging::print(LOG_ERROR,
		       "Not enough memory in axisglib.c ; suspect...");
	PersonalExit(1);
    }

    /* clear temporary array */
    for (nRadial = 0; nRadial < NRadial; nRadial++)
	tempArray[nRadial] = 0.;

    /* computate mean values */
    for (nRadial = Zero_or_active; nRadial < Max_or_active; nRadial++) {
	for (nAzimuthal = 0; nAzimuthal < NAzimuthal; nAzimuthal++) {
	    cell = nRadial * NAzimuthal + nAzimuthal;
	    tempArray[nRadial] += grid[cell];
	}
	tempArray[nRadial] /= (double)NAzimuthal;
    }

    /* if only 1 CPU is calculating everything is easy :) */
    if (CPU_Number == 1) {
	for (nRadial = 0; nRadial < GlobalNRadial; nRadial++)
	    array[nRadial] = tempArray[nRadial];

	/* if more than 1 cpu... */
    } else if (CPU_Number > 1) {
	if (CPU_Rank == 0) {
	    /* first cpu initializes array and put its data in */
	    for (nRadial = 0; nRadial < GlobalNRadial; nRadial++) {
		if (nRadial < Max_or_active)
		    array[nRadial] = tempArray[nRadial];
		else
		    array[nRadial] = 0.;
	    }
	    /* send it to next cpu */
	    MPI_Isend(array, GlobalNRadial, MPI_DOUBLE, CPU_Next, 0,
		      MPI_COMM_WORLD, &req);
	    MPI_Wait(&req, &global_MPI_Status);
	} else if (CPU_Rank != 0) {
	    /* wait for data from previous cpu */
	    MPI_Irecv(array, GlobalNRadial, MPI_DOUBLE, CPU_Prev, 0,
		      MPI_COMM_WORLD, &req);
	    MPI_Wait(&req, &global_MPI_Status);
	    /* add local data */
	    for (nRadial = Zero_or_active; nRadial < Max_or_active; nRadial++)
		array[nRadial + IMIN] = tempArray[nRadial];
	    /* send it to next cpu */
	    if (CPU_Rank != CPU_Highest) {
		MPI_Isend(array, GlobalNRadial, MPI_DOUBLE, CPU_Next, 0,
			  MPI_COMM_WORLD, &req);
		MPI_Wait(&req, &global_MPI_Status);
	    }
	}

	/* synchronize all threads */
	// MPI_Barrier(MPI_COMM_WORLD);

	/* highest cpu has correct array, so broadcast to everyone */
	MPI_Bcast(array, GlobalNRadial, MPI_DOUBLE, CPU_Highest,
		  MPI_COMM_WORLD);
    }

    /* free temporary array */
    free(tempArray);
}
