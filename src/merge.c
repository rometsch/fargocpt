/** \file merge.c

Contains the function that merges the output of different processors.
The resulting merged file is undistinguishable from the file that
would have been produced by a sequential run.

*/

#include "fargo.h"

void merge (nb)
     int nb;
{
  int i,j;
  char radix[512];
  char command[1024];
  if (!CPU_Master) return;
  message ("Merging output files...");
  for (j = 0; j < 3+(AdvecteLabel == YES); j++) {
    switch (j) {
    case 0: strcpy (radix, "dens");
      break;
    case 1: strcpy (radix, "vrad");
      break;
    case 2: strcpy (radix, "vtheta");
      break;
    case 3: strcpy (radix, "label");
      break;
    }
    for (i = 1; i < CPU_Number; i++) {
      sprintf (command, "cd %s; cat gas%s%d.dat.%05d >> gas%s%d.dat",\
	       OUTPUTDIR, radix, nb, i, radix, nb);
      system (command);
    }
    sprintf (command, "cd %s; rm -f gas%s%d.dat.0*",\
	     OUTPUTDIR, radix, nb);
    system (command);
  }
  message ("done\n");
  fflush (stdout);
}
