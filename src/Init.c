/** \file Init.c

Contains the functions needed to initialize the hydrodynamics arrays.
These can be initialized by reading a given output (in the case of a
restart) or by calling a function, InitEuler (), which contains
analytic prescription for the different hydrodynamics fields. Note
that this function InitEuler() is located in SourceEuler.c, which
itself calls InitGas(), in the file Pframeforce.c.
Also, note that the present file contains InitLabel(), which sets
the initial value of a passive scalar.
*/

#include "fargo.h"

extern boolean Restart;
extern int     NbRestart;

void ReadfromFile (array, fileprefix, filenumber)
PolarGrid *array;
char *fileprefix;
int filenumber;
{
  int nr,ns,c, foo=0;
  real *field;
  char name[256];
  FILE *input;
/* Simultaneous read access to the same file have been observed to give wrong results. */
/* A sequential reading is imposed below. */
				/* If current CPU has a predecessor, wait for a message from him */
  if (CPU_Rank > 0) MPI_Recv (&foo, 1, MPI_INT, CPU_Rank-1, 10, MPI_COMM_WORLD, &fargostat);
  sprintf (name, "%s%s%d.dat", OUTPUTDIR, fileprefix, filenumber);
  input = fopen (name, "r");
  if (input == NULL) {
    fprintf (stderr, "WARNING ! Can't read %s. Restarting with t=0 settings.\n", name); 
    if (CPU_Rank < CPU_Number-1) MPI_Send (&foo, 1, MPI_INT, CPU_Rank+1, 10, MPI_COMM_WORLD);
    return;
  }
  field = array->Field;
  nr = array->Nrad;
  ns = array->Nsec;
  for (c = 0; c < IMIN; c++) {
    fread (field, sizeof(real), ns, input); /* Can't read at once in order not to overflow 'field' */
  }
  fread (field, sizeof(real), nr*ns, input);
  fclose (input);
  /* Next CPU is waiting. Tell it to start now by sending the message that it expects */
  if (CPU_Rank < CPU_Number-1) MPI_Send (&foo, 1, MPI_INT, CPU_Rank+1, 10, MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);	/* previous CPUs do not touch anything meanwhile */
}

void InitLabel (array)
PolarGrid *array;
{
  int nr,ns,i,j,l;
  real *field;
  field = array->Field;
  nr = array->Nrad;
  ns = array->Nsec;
  for (i = 0; i <= nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      field[l] = (Rmed[i]-RMIN)/(RMAX-RMIN);
    }
  }
}

void Initialization (gas_density, gas_v_rad, gas_v_theta, gas_label)
PolarGrid *gas_density, *gas_v_rad, *gas_v_theta, *gas_label;
{
  ReadPrevDim ();
  InitEuler (gas_density, gas_v_rad, gas_v_theta);
  InitLabel (gas_label);
  if (Restart == YES) {
    CheckRebin (NbRestart);
    MPI_Barrier (MPI_COMM_WORLD); /* Don't start reading before master has finished rebining... */
				  /* It shouldn't be a problem though since a sequential read is */
                                  /* imposed in the ReadfromFile function below */
    mastererr ("Reading restart files...");
    fflush (stderr);
    ReadfromFile (gas_density, "gasdens", NbRestart);
    ReadfromFile (gas_v_rad, "gasvrad", NbRestart);
    ReadfromFile (gas_v_theta, "gasvtheta", NbRestart);
    ReadfromFile (gas_label, "gaslabel", NbRestart);
    if (StoreSigma) RefillSigma (gas_density);
    fprintf (stderr, "done\n");
    fflush (stderr);
  }
  WriteDim (); 
}
