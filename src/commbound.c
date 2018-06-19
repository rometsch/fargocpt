/** \file commbound.c

Contains the functions used to synchronize buffer zones on all
processes.  In addition to the main function that allows the
synchronization (note that even processes first send their inner zones
to the previous process, then receive their inner buffer from this
process, while odd processes first receive their outer buffer from the
next process, then send their outer zones to the next process. This
file also contains the function that allocates the memory for the
communication (once for all the run).
*/

#include "fargo.h"

static real *SendInnerBoundary;
static real *SendOuterBoundary;
static real *RecvInnerBoundary;
static real *RecvOuterBoundary;

static int allocated_com = 0;
static int size_com;

void AllocateComm () {
  size_com = 3;
  if (AdvecteLabel == YES) size_com = 4;
  size_com *= NSEC * CPUOVERLAP;
  SendInnerBoundary = malloc (size_com * sizeof(real));
  SendOuterBoundary = malloc (size_com * sizeof(real));
  RecvInnerBoundary = malloc (size_com * sizeof(real));
  RecvOuterBoundary = malloc (size_com * sizeof(real));
  if ((SendInnerBoundary == NULL) ||\
      (SendOuterBoundary == NULL) ||\
      (RecvInnerBoundary == NULL) ||\
      (RecvOuterBoundary == NULL)) {
    fprintf (stderr, "CPU %d had not enough memory to allocate communicators.\n", CPU_Rank);
    prs_exit(0);
  }
  allocated_com = 1;
}

void CommunicateBoundaries (Density, Vrad, Vtheta, Label)
     PolarGrid *Density, *Vrad, *Vtheta, *Label;
{
  MPI_Request req1, req2, req3, req4;
  int l, prev, next, oo, o, nr;
  if (!allocated_com) AllocateComm ();
  prev = CPU_Rank-1;
  next = CPU_Rank+1;
  l = CPUOVERLAP*NSEC;
  nr = Density->Nrad;
  oo = (nr-CPUOVERLAP)*NSEC;
  o = (nr-2*CPUOVERLAP)*NSEC;
  memcpy (SendInnerBoundary, Density->Field+l, l*sizeof(real));
  memcpy (SendInnerBoundary+l, Vrad->Field+l, l*sizeof(real));
  memcpy (SendInnerBoundary+2*l, Vtheta->Field+l, l*sizeof(real));
  memcpy (SendOuterBoundary, Density->Field+o, l*sizeof(real));
  memcpy (SendOuterBoundary+l, Vrad->Field+o, l*sizeof(real));
  memcpy (SendOuterBoundary+2*l, Vtheta->Field+o, l*sizeof(real));
  if (AdvecteLabel == YES) {
    memcpy (SendInnerBoundary+3*l, Label->Field+l, l*sizeof(real));
    memcpy (SendOuterBoundary+3*l, Label->Field+o, l*sizeof(real));
  }
  if (CPU_Rank % 2 == 0) {
    if (CPU_Rank > 0) {
      MPI_Isend (SendInnerBoundary, size_com, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, &req1);
      MPI_Irecv (RecvInnerBoundary, size_com, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, &req2);
    }
    if (CPU_Rank < CPU_Number-1) {
      MPI_Isend (SendOuterBoundary, size_com, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &req3);
      MPI_Irecv (RecvOuterBoundary, size_com, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &req4);
    }
  } else {
    if (CPU_Rank < CPU_Number-1) {
      MPI_Irecv (RecvOuterBoundary, size_com, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &req3);
      MPI_Isend (SendOuterBoundary, size_com, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &req4);
    }
    if (CPU_Rank > 0) {
      MPI_Irecv (RecvInnerBoundary, size_com, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, &req1);
      MPI_Isend (SendInnerBoundary, size_com, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, &req2);
    }
  }
  if (CPU_Rank > 0) {
    MPI_Wait (&req1, &fargostat);
    MPI_Wait (&req2, &fargostat);
    memcpy (Density->Field, RecvInnerBoundary, l*sizeof(real));
    memcpy (Vrad->Field, RecvInnerBoundary+l, l*sizeof(real));
    memcpy (Vtheta->Field, RecvInnerBoundary+2*l, l*sizeof(real));
    if (AdvecteLabel == YES) {
      memcpy (Label->Field, RecvInnerBoundary+3*l, l*sizeof(real));
    }
  }
  if (CPU_Rank < CPU_Number-1) {
    MPI_Wait (&req3, &fargostat);
    MPI_Wait (&req4, &fargostat);
    memcpy (Density->Field+oo, RecvOuterBoundary, l*sizeof(real));
    memcpy (Vrad->Field+oo, RecvOuterBoundary+l, l*sizeof(real));
    memcpy (Vtheta->Field+oo, RecvOuterBoundary+2*l, l*sizeof(real));
    if (AdvecteLabel == YES) {
      memcpy (Label->Field+oo, RecvOuterBoundary+3*l, l*sizeof(real));
    }
  }
}
