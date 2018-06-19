/** \file Interpret.c

Contains the functions required to read the parameter file,
and functions that provide runtime information. The function
var() associates a string to a global variable. The function
ReadVariables() reads the content of a parameter file.
In addition, this file contains a function that prints the
command line usage to the standard output, a function that
provides verbose information about the setup (if the -v switch
is set on the command line), and functions that act as a
chronometer (if the -t switch is set on the command line). 
*/

#include "fargo.h"
#define MAXVARIABLES 500

extern int      begin_i;
extern boolean  OpenInner;
static Param    VariableSet[MAXVARIABLES];
static int      VariableIndex = 0;
static int	FirstStep = YES;
static clock_t  First, Preceeding, Current, FirstUser, CurrentUser, PreceedingUser;
static long	Ticks;
boolean         FastTransport = YES, GuidingCenter = NO;
boolean         IsDisk = YES, NonReflecting = NO, Corotating = NO, OuterSourceMass = NO;
boolean         Write_Density = YES, Write_Velocity = YES, Indirect_Term = YES;

void
var(name, ptr, type, necessary, deflt)
char           *name;
char           *ptr;
int             type;
int             necessary;
char           *deflt;
{
  real            valuer;
  int             valuei;
  double	  temp;
  sscanf (deflt, "%lf", &temp);
  valuer = (real) (temp);
  valuei = (int) valuer;
  strcpy(VariableSet[VariableIndex].name, name);
  VariableSet[VariableIndex].variable = ptr;
  VariableSet[VariableIndex].type = type;
  VariableSet[VariableIndex].necessary = necessary;
  VariableSet[VariableIndex].read = NO;
  if (necessary == NO) {
    if (type == INT) {
      *((int *) ptr) = valuei;
    } else if (type == REAL) {
      *((real *) ptr) = valuer;
    } else if (type == STRING) {
      strcpy (ptr, deflt);
    }
  }
  VariableIndex++;
}

void
ReadVariables(filename)
char *filename;
{
  char            nm[300], s[350],stringval[290];
  char           *s1;
  double	  temp;
  real            valuer;
  int             i, found, valuei, success, type;
  int            *ptri;
  real           *ptrr;
  FILE           *input;

  InitVariables();
  input = fopen(filename, "r");
  if (input == NULL) {
      mastererr ("Unable to read '%s'. Program stopped.\n",filename);
    prs_exit(1);
  }
  mastererr ("Reading parameters file '%s'.\n", filename);
  while (fgets(s, 349, input) != NULL) {
    success = sscanf(s, "%s ", nm);
    if ((nm[0] != '#') && (success == 1)) {	/* # begins a comment
						 * line */
      s1 = s + strlen(nm);
      sscanf(s1 + strspn(s1, "\t :=>_"), "%lf", &temp);
      sscanf(s1 + strspn(s1, "\t :=>_"), "%289s ", stringval);
      valuer = (real) temp;
      valuei = (int) temp;
      for (i = 0; i < strlen(nm); i++) {
	nm[i] = (char) toupper(nm[i]);
      }
      found = NO;
      for (i = 0; i < VariableIndex; i++) {
	if (strcmp(nm, VariableSet[i].name) == 0) {
	  if (VariableSet[i].read == YES) {
	    mastererr("Warning : %s defined more than once.\n", nm);
	  }
	  found = YES;
	  VariableSet[i].read = YES;
	  ptri = (int *) (VariableSet[i].variable);
	  ptrr = (real *) (VariableSet[i].variable);
	  if (VariableSet[i].type == INT) {
	    *ptri = valuei;
	  } else if (VariableSet[i].type == REAL) {
	    *ptrr = valuer;
	  } else if (VariableSet[i].type == STRING) {
	    strcpy (VariableSet[i].variable, stringval);
	  }
	}
      }
      if (found == NO) {
	mastererr("Warning : variable %s defined but non-existent in code.\n", nm);
      }
    }
  }

  found = NO;
  for (i = 0; i < VariableIndex; i++) {
    if ((VariableSet[i].read == NO) && (VariableSet[i].necessary == YES)) {
      if (found == NO) {
	mastererr("Fatal error : undefined mandatory variable(s):\n");
	found = YES;
      }
      mastererr("%s\n", VariableSet[i].name);
    }
    if (found == YES)
      prs_exit(1);

  }
  found = NO;
  for (i = 0; i < VariableIndex; i++) {
    if (VariableSet[i].read == NO) {
      if (found == NO) {
	mastererr("Secondary variables omitted :\n");
	found = YES;
      }
      if ((type = VariableSet[i].type) == REAL)
	mastererr("%s ;\t Default Value : %.5g\n", VariableSet[i].name, *((real *) VariableSet[i].variable));
      if (type == INT)
	mastererr("%s ;\t Default Value : %d\n", VariableSet[i].name, *((int *) VariableSet[i].variable));
      if (type == STRING)
	mastererr("%s ;\t Default Value : %s\n", VariableSet[i].name, VariableSet[i].variable);
    }
  }
  if ((*ADVLABEL == 'y') || (*ADVLABEL == 'Y')) AdvecteLabel = YES;
  if ((*OUTERSOURCEMASS == 'y') || (*OUTERSOURCEMASS == 'Y')) OuterSourceMass = YES;
  if ((*TRANSPORT == 's') || (*TRANSPORT == 'S')) FastTransport = NO;
  if ((*OPENINNERBOUNDARY == 'O') || (*OPENINNERBOUNDARY == 'o')) OpenInner = YES;
  if ((*OPENINNERBOUNDARY == 'N') || (*OPENINNERBOUNDARY == 'n')) NonReflecting = YES;
  if ((*GRIDSPACING == 'L') || (*GRIDSPACING == 'l')) LogGrid = YES;
  if ((*DISK == 'N') || (*DISK == 'n')) IsDisk = NO;
  if ((*FRAME == 'C') || (*FRAME == 'c')) Corotating = YES;
  if ((*FRAME == 'G') || (*FRAME == 'g')) {
    Corotating = YES;
    GuidingCenter = YES;
  }
  if ((*WRITEVELOCITY == 'N') || (*WRITEVELOCITY == 'n')) Write_Velocity = NO;
  if ((*WRITEDENSITY == 'N') || (*WRITEDENSITY == 'n')) Write_Density = NO;
  if ((*INDIRECTTERM == 'N') || (*INDIRECTTERM == 'n')) Indirect_Term = NO;
  if ((*EXCLUDEHILL == 'Y') || (*EXCLUDEHILL == 'y')) ExcludeHill = YES;
  if ((ALPHAVISCOSITY != 0.0) && (VISCOSITY != 0.0)) {
    mastererr ("You cannot use at the same time\n");
    mastererr ("VISCOSITY and ALPHAVISCOSITY.\n");
    mastererr ("Edit the parameter file so as to remove\n");
    mastererr ("one of these variables and run again.\n");
    prs_exit (1);
  }
  if (ALPHAVISCOSITY != 0.0) {
    ViscosityAlpha = YES;
    masterprint ("Viscosity is of alpha type\n");
  }
  if ((THICKNESSSMOOTHING != 0.0) && (ROCHESMOOTHING != 0.0)) {
    mastererr ("You cannot use at the same time\n");
    mastererr ("`ThicknessSmoothing' and `RocheSmoothing'.\n");
    mastererr ("Edit the parameter file so as to remove\n");
    mastererr ("one of these variables and run again.\n");
    prs_exit (1);
  }
  if ((THICKNESSSMOOTHING <= 0.0) && (ROCHESMOOTHING <= 0.0)) {
    mastererr ("A non-vanishing potential smoothing length is required.\n");
    mastererr ("Please use either of the following variables:\n");
    mastererr ("`ThicknessSmoothing' *or* `RocheSmoothing'.\n");
    mastererr ("before launching the run again.\n");
    prs_exit (1);
  }
  if (ROCHESMOOTHING != 0.0) {
    RocheSmoothing = YES;
    masterprint ("Planet potential smoothing scales with their Hill sphere.\n");
  }
  if (OverridesOutputdir == YES) {
    sprintf (OUTPUTDIR, "%s", NewOutputdir);
  }
				/* Add a trailing slash to OUTPUTDIR if needed */
  if (*(OUTPUTDIR+strlen(OUTPUTDIR)-1) != '/')
    strcat (OUTPUTDIR, "/");
}

void PrintUsage (execname)
char *execname;
{
  mastererr("Usage : %s [-abcdeimnptvz] [-(0-9)] [-s number] [-f scaling] parameters file\n", execname);
  mastererr("\n-a : Monitor mass and angular momentum at each timestep\n");
  mastererr("-b : Adjust azimuthal velocity to impose strict centrifugal balance at t=0\n");
  mastererr("-c : Sloppy CFL condition (checked at each DT, not at each timestep)\n");
  mastererr("-d : Print some debugging information on 'stdout' at each timestep\n");
  mastererr("-e : Activate EU test problem torque file output\n");
  mastererr("-f : Scale density array by 'scaling'. Useful to increase/decrease\n");
  mastererr("     disk surface density after a restart, for instance.            \n");
  mastererr("-i : tabulate Sigma profile as given by restart files\n");
  mastererr("-m : Merge output files from different CPUs\n");
  mastererr("-n : Disable simulation. The program just reads parameters file\n");
  mastererr("-o : Overrides output directory of input file.\n");
  mastererr("-p : Give profiling information at each time step\n");
  mastererr("-s : Restart simulation, taking #'number' files as initial conditions\n");
  mastererr("-t : Monitor CPU time usage at each time step\n");
  mastererr("-v : Verbose mode. Tells everything about parameters file\n");
  mastererr("-z : fake sequential built when evaluating sums on HD meshes\n");
  mastererr("-(0-9) : only write initial (or restart) HD meshes,\n");
  mastererr("     proceed to the next nth output and exit\n");
  mastererr("     This option must stand alone on one switch (-va -4 is legal, -v4a is not)\n");
  prs_exit (1);
}

real TellNbOrbits (time)
real time;
{
  return time/2.0/PI*sqrt(G*1.0/1.0/1.0/1.0);
}

real TellNbOutputs (time)
real time;
{
  return (time/DT/NINTERM);
}

void TellEverything () {
  real temp;
  if (!CPU_Master) return;
  printf ("\nDisc properties:\n");
  printf ("----------------\n");
  printf ("Inner Radius          : %g\n", RMIN);
  printf ("Outer Radius          : %g\n", RMAX);
  printf ("Aspect Ratio          : %g\n", ASPECTRATIO);
  printf ("VKep at inner edge    : %.3g\n", sqrt(G*1.0*(1.-0.0)/RMIN));
  printf ("VKep at outer edge    : %.3g\n", sqrt(G*1.0/RMAX));
  temp=SIGMA0*PI*(RMAX*RMAX-RMIN*RMIN);	/* correct this and what follows... */
  printf ("Disk Mass             : %g\n", temp);
  temp=SIGMA0*PI*(1.0-RMIN*RMIN);
  printf ("Mass inner to r=1.0  : %g \n", temp);
  temp=SIGMA0*PI*(RMAX*RMAX-1.0);
  printf ("Mass outer to r=1.0  : %g \n", temp);
  printf ("Travelling time for acoustic density waves :\n");
  temp = 2.0/3.0/ASPECTRATIO*(pow(RMAX,1.5)-pow(RMIN,1.5));
  printf (" * From Rmin to Rmax  : %.2g = %.2f orbits ~ %.1f outputs\n", temp, TellNbOrbits(temp), TellNbOutputs(temp));
  temp = 2.0/3.0/ASPECTRATIO*(pow(RMAX,1.5)-pow(1.0,1.5));
  printf (" * From r=1.0 to Rmax: %.2g = %.2f orbits ~ %.1f outputs\n", temp, TellNbOrbits(temp), TellNbOutputs(temp));
  temp = 2.0/3.0/ASPECTRATIO*(pow(1.0,1.5)-pow(RMIN,1.5));
  printf (" * From r=1.0 to Rmin: %.2g = %.2f orbits ~ %.1f outputs\n", temp, TellNbOrbits(temp), TellNbOutputs(temp));
  temp = 2.0*PI*sqrt(RMIN*RMIN*RMIN/G/1.0);
  printf ("Orbital time at Rmin  : %.3g ~ %.2f outputs\n", temp, TellNbOutputs(temp));
  temp = 2.0*PI*sqrt(RMAX*RMAX*RMAX/G/1.0);
  printf ("Orbital time at Rmax  : %.3g ~ %.2f outputs\n", temp, TellNbOutputs(temp));
  printf ("Sound speed :\n");
  printf (" * At unit radius     : %.3g\n", ASPECTRATIO*sqrt(G*1.0));
  printf (" * At outer edge      : %.3g\n", ASPECTRATIO*sqrt(G*1.0/RMAX));
  printf (" * At inner edge      : %.3g\n", ASPECTRATIO*sqrt(G*1.0/RMIN));
  printf ("\nGrid properties:\n");
  printf ("----------------\n");
  printf ("Number of rings       : %d\n", NRAD);
  printf ("Number of sectors     : %d\n", NSEC);
  printf ("Total cells           : %d\n", NRAD*NSEC);
  printf ("\nOutputs properties:\n");
  printf ("-------------------\n");
  printf ("Time increment between outputs : %.3f = %.3f orbits\n", NINTERM*DT, TellNbOrbits(NINTERM*DT));
  printf ("At each output #i, the following files are written:\n");
  printf ("gasdens[i].dat       : %d bytes\n",(int)(NRAD*NSEC*sizeof(real)));
  printf ("gasvrad[i].dat       : %d bytes\n",(int)(NRAD*NSEC*sizeof(real)));
  printf ("gasvtheta[i].dat     : %d bytes\n",(int)(NRAD*NSEC*sizeof(real)));
  if (AdvecteLabel == YES)
    printf ("gaslabel[i].dat      : %d bytes\n",(int)(NRAD*NSEC*sizeof(real)));
  printf ("There will be in total %d outputs\n", NTOT/NINTERM);
  printf ("(which correspond to an elapsed time = %.3f or to %.2f orbits)\n", NTOT*DT, TellNbOrbits(NTOT*DT));
  temp =3.0*NRAD*NSEC*sizeof(real);
  if (AdvecteLabel == YES)
    temp *= 4.0/3.0;
  temp *= (real)NTOT/(real)NINTERM;
  temp /= 1024.0*1024.0;
  printf ("So the code will produce ~%.2f Mbytes of data\n", temp);
  printf ("Check (eg by issuing a 'df' command) that you have enough disk space,\n");
  printf ("otherwise you will get a system full and the code will stop.\n");
  fflush (stdout);
}

void GiveTimeInfo (number)
int number;
{
  struct tms buffer;
  real total, last, mean, totalu;
  Current = times (&buffer);
  CurrentUser = buffer.tms_utime;
  if (FirstStep == YES) {
    First = Current;
    FirstUser = CurrentUser;
    fprintf (stderr, "Time counters initialized\n");
    FirstStep = NO;
    Ticks = sysconf (_SC_CLK_TCK);
  }
  else {
    total = (real)(Current - First)/Ticks;
    totalu= (real)(CurrentUser-FirstUser)/Ticks;
    last  = (real)(CurrentUser - PreceedingUser)/Ticks;
    number -= begin_i/NINTERM;
    mean  = totalu / number;
    fprintf (stderr, "Total Real Time elapsed    : %.3f s\n", total);
    fprintf (stderr, "Total CPU Time of process  : %.3f s (%.1f %%)\n", totalu, 100.*totalu/total);
    fprintf (stderr, "CPU Time since last time step : %.3f s\n", last);
    fprintf (stderr, "Mean CPU Time between time steps : %.3f s\n", mean);
    fprintf (stderr, "CPU Load on last time step : %.1f %% \n", (real)(CurrentUser-PreceedingUser)/(real)(Current-Preceeding)*100.);

  }	
  PreceedingUser = CurrentUser;
  Preceeding = Current;
}

void InitSpecificTime (profiling, process_name, title)
boolean profiling;
TimeProcess *process_name;
char *title;
{
  struct tms buffer;
  if (profiling == NO) return;
  Ticks = sysconf (_SC_CLK_TCK);
  times (&buffer);
  process_name->clicks = buffer.tms_utime;
  strcpy (process_name->name, title);
}

void GiveSpecificTime (profiling, process_name)
boolean profiling;
TimeProcess process_name;
{
  struct tms buffer;
  long ticks;
  real t;
  if (profiling == NO) return;
  Ticks = sysconf (_SC_CLK_TCK);
  times (&buffer);
  ticks = buffer.tms_utime - process_name.clicks;
  t = (real)ticks / (real)Ticks;
  fprintf (stderr, "Time spent in %s : %.3f s\n", process_name.name, t);
}

