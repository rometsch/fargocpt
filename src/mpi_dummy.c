/** \file mpi_dummy.c

Fake MPI functions library for sequential built.
It is used instead of the true MPI library in the
case of a sequential built (see makefile).

*/


#include <stdio.h>
#include "mpi_dummy.h"

void MPI_Comm_rank (a, b)
     int a;
     int *b;
{
  *b = 0;			/* Only one process, with rank zero... */
}

void MPI_Comm_size (a, b)
     int a;
     int *b;
{
  *b = 1;			/* Only one process in the world communicator... */
}

void MPI_Init (argc, argv)
     int *argc;
     char **argv[];
{
  fprintf (stderr, "\n       !!!! WARNING !!!!\n\n");
  fprintf (stderr, "This is a sequential built of the %s code\n", *argv[0]);
  fprintf (stderr, "If you planned to run the MPI-parallel version,\n");
  fprintf (stderr, "then you MUST rebuild the executable. Go into the\n");
  fprintf (stderr, "source directory (normally src/), then issue:\n");
  fprintf (stderr, "\ngmake BUILD=parallel\n");
  fprintf (stderr, "\nAny further invocation of gmake will refer to a parallel built.\n");
}

void MPI_Finalize ()
{
}

void MPI_Bcast ()
{
}

void MPI_Isend ()
{
}

void MPI_Irecv ()
{
}

void MPI_Send ()
{
}

void MPI_Recv ()
{
}

void MPI_Barrier ()
{
}

void MPI_Wait ()
{
}

void MPI_Allreduce (void *ptr, void *ptr2, int count, int type, int foo3, int foo4)
{
  int i;
  for (i = 0; i < count; i++) {
    switch (type) {
    case MPI_DOUBLE:
      *(((double *)ptr2)+i) = (double)(*(((double *)ptr)+i));
      break;
    case MPI_INT:
      *(((int *)ptr2)+i) = (int)(*(((int *)ptr)+i));
      break;
    }
  }
}
