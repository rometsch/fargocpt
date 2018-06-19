/** \file mpi_dummy.h

Declaration of fake MPI functions for sequential built.
There are used instead of the true MPI functions in the
case of a sequential built (see makefile).

*/


#define MPI_COMM_WORLD 0
#define MPI_DOUBLE 2
#define MPI_CHAR 1
#define MPI_LONG 3
#define MPI_INT 0
#define MPI_MIN 0
#define MPI_MAX 0
#define MPI_SUM 0

typedef int MPI_Request;
typedef int MPI_Status;

void MPI_Comm_rank ();
void MPI_Barrier ();
void MPI_Comm_size ();
void MPI_Init ();
void MPI_Finalize ();
void MPI_Bcast ();
void MPI_Isend ();
void MPI_Irecv ();
void MPI_Allreduce ();
void MPI_Send ();
void MPI_Recv ();
void MPI_Wait ();

