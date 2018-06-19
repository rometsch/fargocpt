/** \file rebin.c

Resample the hydrodynamical fields at a restart with a different resolution.
Note that at the restart, even if NRAD and NSEC coincide with previous values,
the data is resampled if the radii do not coincide (for instance, if we switch
from ARITHMETIC to LOGARITHMIC spacing).

*/


#include "fargo.h"

static real OldRadii[MAX1D], OldRmed[MAX1D], New_r[MAX1D];
static int OldNRAD, OldNSEC;

void ReadPrevDim () {
  FILE *DIM, *RAD;
  char name_dim[1024], name_rad[1024];
  int i, foo;
  float Foo;
  double value;
  if (!CPU_Master) return;
  OldNRAD = 0;
  sprintf (name_dim, "%sdims.dat", OUTPUTDIR);
  sprintf (name_rad, "%sused_rad.dat", OUTPUTDIR);
  DIM = fopen (name_dim, "r");
  if (DIM == NULL) return;
  RAD = fopen (name_rad, "r");
  if (RAD == NULL) return;
  fscanf (DIM,"%d %d %d %d %f %d %d %d\n",&foo,&foo,&foo,&foo,&Foo,&foo,&OldNRAD,&OldNSEC);
  for (i = 0; i <= OldNRAD; i++) {
    fscanf (RAD, "%lf", &value);
    OldRadii[i] = (real)value;
  }
  fclose (DIM);
  fclose (RAD);
  for (i = 0; i < OldNRAD; i++) {
    OldRmed[i] = 2.0/3.0*(OldRadii[i+1]*OldRadii[i+1]*OldRadii[i+1]-OldRadii[i]*OldRadii[i]*OldRadii[i]);
    OldRmed[i] = OldRmed[i] / (OldRadii[i+1]*OldRadii[i+1]-OldRadii[i]*OldRadii[i]);
  }
}

void CheckRebin (nb) 
int nb;
{
  boolean RebinNeeded = NO, found;
  char radix[1024], filename[1024];
  int i, iold, jold, l, j, lo, loip, lojp, loipjp, type;
  real ifrac, jfrac, jreal, r, angle, dangle;
  FILE *ARR;
  real *Old_r, *OldArray, *NewArray;
  if (!CPU_Master) return;
  if (NSEC != OldNSEC) RebinNeeded = YES;
  if (GLOBALNRAD != OldNRAD) RebinNeeded = YES;
  for (i = 0; i <= NRAD; i++) {
    if (fabs((Radii[i]-OldRadii[i])/Radii[i]) > 1e-9) RebinNeeded = YES;
  }
  if (!RebinNeeded) return;
  printf ("Restart/Old mesh mismatch. Rebin needed.\n");
  OldArray = (real *)malloc(OldNRAD*OldNSEC*sizeof(real));
  NewArray = (real *)malloc(GLOBALNRAD*NSEC*sizeof(real));
  if ((OldArray == NULL) || (NewArray == NULL)) {
    mastererr ("Not enough memory left for rebining.\n");
    mastererr ("Aborted.\n");
    prs_exit (1);
  }
  for (type = 0; type < 4; type++) {
    Old_r = OldRmed;
    memcpy (New_r, GlobalRmed, (GLOBALNRAD+1)*sizeof(double));
    dangle = 0.0;
    switch (type) {
    case 0: sprintf (radix, "dens"); break;
    case 1: sprintf (radix, "vrad"); 
      Old_r = OldRadii; 
      memcpy (New_r, Radii, (GLOBALNRAD+1)*sizeof(double));
      break;
    case 2: sprintf (radix, "vtheta"); dangle = 0.5; break;
    case 3: sprintf (radix, "label"); break;
    }
    for (i = 0; i < GLOBALNRAD; i++) {
      if (New_r[i] < Old_r[0]) New_r[i] = Old_r[0];
      if (New_r[i] > Old_r[OldNRAD-1]) New_r[i] = Old_r[OldNRAD-1];
    }
    sprintf (filename, "%sgas%s%d.dat", OUTPUTDIR, radix, nb);
    ARR = fopen (filename, "r");
    if (ARR != NULL) {
      fread (OldArray, sizeof(real), OldNRAD*OldNSEC, ARR);
      fclose (ARR);
      for (i = 0; i < GLOBALNRAD; i++) {
	r = New_r[i];
	iold = 0;
	found = NO;
	while ((iold < OldNRAD) && (!found)) {
	  if (Old_r[iold+1] <= r) iold++;
	  else found = YES;
	}
	if (r <= Old_r[0]) {
	  iold = 0;
	  ifrac = 0.0;
	} else if (r >= Old_r[OldNRAD-1]) {
	  iold = OldNRAD-2;
	  ifrac = 1.0;
	} else
	  ifrac = (r-Old_r[iold])/(Old_r[iold+1]-Old_r[iold]);
	printf ("%d\t%d\t%g\n", i, iold, ifrac);
	for (j = 0; j < NSEC; j++) {
	  l = j+i*NSEC;
	  angle = ((real)j-dangle)/(real)NSEC*2.0*M_PI;
	  jreal = angle/2.0/M_PI*(real)OldNSEC+dangle;
	  while (jreal < 0.0) jreal += (real)OldNSEC;
	  jold = (int)jreal;
	  jold = jold % OldNSEC;
	  jfrac = jreal-(real)jold;
	  lo = jold+iold*OldNSEC;
	  loip = lo+OldNSEC;
	  lojp = (jold == OldNSEC-1 ? lo-OldNSEC+1 : lo+1);
	  loipjp = lojp+OldNSEC;
	  NewArray[l] = OldArray[lo]*(1.0-ifrac)*(1.-jfrac)+\
	    OldArray[lojp]*(1.0-ifrac)*jfrac+\
	    OldArray[loip]*ifrac*(1.0-jfrac)+\
	    OldArray[loipjp]*ifrac*jfrac;
	}
      }
      ARR = fopenp (filename, "w");
      fwrite (NewArray, sizeof(real), NRAD*NSEC, ARR);
      fclose (ARR);
    }
    else {
      mastererr("Could not rebin %s. File not found\n", filename);
    }
  }
}
