/**
	\file commbound.c

	Functions for communication of boundaries between CPUs
*/

#include <string.h>

#include "LowTasks.h"
#include "commbound.h"
#include "constants.h"
#include "global.h"
#include "logging.h"
#include "parameters.h"
#include "time.h"

/// pointer to inner boundary send puffer
static double *SendInnerBoundary;

/// pointer to outer boundary send buffer
static double *SendOuterBoundary;

/// pointer to inner boundary receive buffer
static double *RecvInnerBoundary;

/// pointer to outer boundary receive buffer
static double *RecvOuterBoundary;

/// size of each buffer for boundary communication
static size_t bufferSize;

/// flag if buffers for boundary communication already have been allocated
static int buffersAllocated = 0;

void DeallocateBoundaryCommunicationBuffers()
{
    /* deallocate the buffers */
    free(SendInnerBoundary);
    free(SendOuterBoundary);
    free(RecvInnerBoundary);
    free(RecvOuterBoundary);
}

/**
	Calculate size of buffers needed for boundary communication and allocate
   them.
*/
void AllocateBoundaryCommunicationBuffers()
{
    /* number of polar grids to communicate; standard is 3 (density, vrad,
     * vtheta) */
    bufferSize = 3;

    /* if we calculate adiabatic, energy is additional needed */
    if (parameters::Adiabatic) {
	bufferSize++;
    }

    /* each variable has to be stored NAzimuthal times overlap cells. */
    bufferSize *= NAzimuthal * CPUOVERLAP;

    /* allocate the buffers */
    SendInnerBoundary = (double *)malloc(bufferSize * sizeof(double));
    SendOuterBoundary = (double *)malloc(bufferSize * sizeof(double));
    RecvInnerBoundary = (double *)malloc(bufferSize * sizeof(double));
    RecvOuterBoundary = (double *)malloc(bufferSize * sizeof(double));

    if ((SendInnerBoundary == NULL) || (SendOuterBoundary == NULL) ||
	(RecvInnerBoundary == NULL) || (RecvOuterBoundary == NULL)) {
	logging::print(
	    LOG_ERROR
	    "CPU %d had not enough memory to allocate communicators.\n",
	    CPU_Rank);
	PersonalExit(0);
    }

    // remember that we already allocated buffers
    buffersAllocated = 1;
}

/**
	Communicate boundaries.

	\param Density
	\param Vrad
	\param Vtheta
	\param Energy
*/
void CommunicateBoundaries(t_polargrid *Density, t_polargrid *Vrad,
			   t_polargrid *Vtheta, t_polargrid *Energy)
{
    MPI_Request req1, req2, req3, req4;

    // check if buffers have already been allocted
    if (!buffersAllocated) {
	AllocateBoundaryCommunicationBuffers();
    }

    ptrdiff_t l = CPUOVERLAP * NAzimuthal;
    ptrdiff_t oo = (Density->Nrad - CPUOVERLAP) * NAzimuthal;
    ptrdiff_t o = (Density->Nrad - 2 * CPUOVERLAP) * NAzimuthal;

    // copy data into send buffers
    memcpy(SendInnerBoundary, Density->Field + l, l * sizeof(double));
    memcpy(SendInnerBoundary + l, Vrad->Field + l, l * sizeof(double));
    memcpy(SendInnerBoundary + 2 * l, Vtheta->Field + l, l * sizeof(double));
    memcpy(SendOuterBoundary, Density->Field + o, l * sizeof(double));
    memcpy(SendOuterBoundary + l, Vrad->Field + o, l * sizeof(double));
    memcpy(SendOuterBoundary + 2 * l, Vtheta->Field + o, l * sizeof(double));

    if (parameters::Adiabatic) {
	memcpy(SendInnerBoundary + 3 * l, Energy->Field + l,
	       l * sizeof(double));
	memcpy(SendOuterBoundary + 3 * l, Energy->Field + o,
	       l * sizeof(double));
    }

    /* Note that boundary exchange is independant from chosen domain
     * decomposition */
    /* send / receive data */
    if (CPU_Rank % 2 == 0) {
	if (CPU_Rank != 0) {
	    MPI_Isend(SendInnerBoundary, bufferSize, MPI_DOUBLE, CPU_Prev, 0,
		      MPI_COMM_WORLD, &req1);
	    MPI_Irecv(RecvInnerBoundary, bufferSize, MPI_DOUBLE, CPU_Prev, 0,
		      MPI_COMM_WORLD, &req2);
	}

	if (CPU_Rank != CPU_Highest) {
	    MPI_Isend(SendOuterBoundary, bufferSize, MPI_DOUBLE, CPU_Next, 0,
		      MPI_COMM_WORLD, &req3);
	    MPI_Irecv(RecvOuterBoundary, bufferSize, MPI_DOUBLE, CPU_Next, 0,
		      MPI_COMM_WORLD, &req4);
	}
    } else {
	if (CPU_Rank != CPU_Highest) {
	    MPI_Irecv(RecvOuterBoundary, bufferSize, MPI_DOUBLE, CPU_Next, 0,
		      MPI_COMM_WORLD, &req3);
	    MPI_Isend(SendOuterBoundary, bufferSize, MPI_DOUBLE, CPU_Next, 0,
		      MPI_COMM_WORLD, &req4);
	}

	if (CPU_Rank != 0) {
	    MPI_Irecv(RecvInnerBoundary, bufferSize, MPI_DOUBLE, CPU_Prev, 0,
		      MPI_COMM_WORLD, &req1);
	    MPI_Isend(SendInnerBoundary, bufferSize, MPI_DOUBLE, CPU_Prev, 0,
		      MPI_COMM_WORLD, &req2);
	}
    }

    if (CPU_Rank != 0) {
	MPI_Wait(&req1, &global_MPI_Status);
	MPI_Wait(&req2, &global_MPI_Status);
	memcpy(Density->Field, RecvInnerBoundary, l * sizeof(double));
	memcpy(Vrad->Field, RecvInnerBoundary + l, l * sizeof(double));
	memcpy(Vtheta->Field, RecvInnerBoundary + 2 * l, l * sizeof(double));
	if (parameters::Adiabatic)
	    memcpy(Energy->Field, RecvInnerBoundary + 3 * l,
		   l * sizeof(double));
    }

    if (CPU_Rank != CPU_Highest) {
	MPI_Wait(&req3, &global_MPI_Status);
	MPI_Wait(&req4, &global_MPI_Status);
	memcpy(Density->Field + oo, RecvOuterBoundary, l * sizeof(double));
	memcpy(Vrad->Field + oo, RecvOuterBoundary + l, l * sizeof(double));
	memcpy(Vtheta->Field + oo, RecvOuterBoundary + 2 * l,
	       l * sizeof(double));
	if (parameters::Adiabatic)
	    memcpy(Energy->Field + oo, RecvOuterBoundary + 3 * l,
		   l * sizeof(double));
    }
}
