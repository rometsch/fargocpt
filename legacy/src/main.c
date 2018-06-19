/** \file main.c

Main file of the distribution. Manages the call to initialization
functions, then the main loop.

*/

#include "fargo.h"

boolean         Restart = NO, OpenInner = NO;
int             begin_i = 0, NbRestart = 0;
static int      InnerOutputCounter=0, StillWriteOneOutput;

extern real     LostMass;
extern boolean  Corotating;
real            ScalingFactor = 1.0;

int
main(argc, argv)
int argc;
char *argv[];
{
  PolarGrid   *gas_density, *gas_v_rad, *gas_v_theta, *gas_label;
  int          i;
  boolean      disable = NO, TimeInfo = NO, Profiling = NO;
  boolean      TimeToWrite, verbose = NO;
  TimeProcess  t_Hydro;
  char         ParameterFile[256];
  PlanetarySystem *sys;

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &CPU_Rank);
  MPI_Comm_size (MPI_COMM_WORLD, &CPU_Number);
  CPU_Master = (CPU_Rank == 0 ? 1 : 0);
  setfpe ();			/* Control behavior for floating point
				   exceptions trapping */
  if (argc == 1) PrintUsage (argv[0]);
  strcpy (ParameterFile, "");
  for (i = 1; i < argc; i++) {
    if (*(argv[i]) == '-') {
      if (strspn (argv[i], "-secndovtpfamzib0123456789") != strlen (argv[i]))
	PrintUsage (argv[0]);
      if (strchr (argv[i], 'n'))
	disable = YES;
      if (strchr (argv[i], 'e'))
	Stockholm = YES;
      if (strchr (argv[i], 'v'))
	verbose = YES;
      if (strchr (argv[i], 't'))
	TimeInfo = YES;
      if (strchr (argv[i], 'c'))
	SloppyCFL = YES;
      if (strchr (argv[i], 'p'))
	Profiling = YES;
      if (strchr (argv[i], 'd'))
	debug = YES;
      if (strchr (argv[i], 'b'))
	CentrifugalBalance = YES;
      if (strchr (argv[i], 'm'))
	Merge = YES;
      if (strchr (argv[i], 'a'))
	MonitorIntegral = YES;
      if (strchr (argv[i], 'z'))
	FakeSequential = YES;
      if (strchr (argv[i], 'i'))
	StoreSigma = YES;
      if (strchr (argv[i], '0'))
	OnlyInit = YES;
      if ((argv[i][1] >= '1') && (argv[i][1] <= '9')) {
	GotoNextOutput = YES;
	StillWriteOneOutput = (int)(argv[i][1]-'0');
      }
      if (strchr (argv[i], 's')) {
	Restart = YES;
	i++;
	NbRestart = atoi(argv[i]);
	if ((NbRestart < 0)) {
	  masterprint ("Incorrect restart number\n");
	  PrintUsage (argv[0]);
	}
      }
      if (strchr (argv[i], 'o')) {
	OverridesOutputdir = YES;
	i++;
	sprintf (NewOutputdir, "%s", argv[i]);
      } else {
	if (strchr (argv[i], 'f')) {
	  i++;
	  ScalingFactor = atof(argv[i]);
	  masterprint ("Scaling factor = %g\n", ScalingFactor);
	  if ((ScalingFactor <= 0)) {
	    masterprint ("Incorrect scaling factor\n");
	    PrintUsage (argv[0]);
	  }
	}
      }
    }
    else strcpy (ParameterFile, argv[i]);
  }
  if ((StoreSigma) && !(Restart)) {
    mastererr ("You cannot use tabulated surface density\n");
    mastererr ("in a non-restart run. Aborted\n");
    prs_exit (0);
  }
  if (ParameterFile[0] == 0) PrintUsage (argv[0]);
  ReadVariables (ParameterFile);
  sys = InitPlanetarySystem (PLANETCONFIG);
  ListPlanets (sys);
  SplitDomain ();
  if (verbose == YES) 
    TellEverything ();
  if (disable == YES)
    prs_exit (0);
  MakeDir (OUTPUTDIR);
  DumpSources (argc, argv);
  masterprint ("Allocating arrays...");
  fflush (stdout);
  gas_density        = CreatePolarGrid(NRAD, NSEC, "dens");
  gas_v_rad          = CreatePolarGrid(NRAD, NSEC, "vrad");
  gas_v_theta        = CreatePolarGrid(NRAD, NSEC, "vtheta");
  gas_label          = CreatePolarGrid(NRAD, NSEC, "label");
  masterprint ("done.\n");
  OmegaFrame = OMEGAFRAME;
  if (Corotating == YES) OmegaFrame = GetPsysInfo (sys, FREQUENCY);
  Initialization (gas_density, gas_v_rad, gas_v_theta, gas_label);
  InitComputeAccel ();
  if (Restart == YES) {
    begin_i         = NbRestart * NINTERM;
    RestartPlanetarySystem (NbRestart, sys);
    LostMass = GetfromPlanetFile (NbRestart, 7, 0); /* 0 refers to planet #0 */
    PhysicalTime  = GetfromPlanetFile (NbRestart, 8, 0);
    OmegaFrame  = GetfromPlanetFile (NbRestart, 9, 0);
  } else {			/* We initialize 'planet[i].dat' file */
    EmptyPlanetSystemFile (sys);
  }
  PhysicalTimeInitial = PhysicalTime;
  MultiplyPolarGridbyConstant (gas_density, ScalingFactor);
  for (i = begin_i; i <= NTOT; i++) {
    InnerOutputCounter++;
    if (InnerOutputCounter == 1) {
      InnerOutputCounter = 0;
      WriteBigPlanetSystemFile (sys, TimeStep);
      UpdateLog (sys, gas_density, PhysicalTime);
      if (Stockholm == YES)
	UpdateLogStockholm (sys, gas_density, TimeStep, PhysicalTime);
    }
    if (NINTERM * (TimeStep = (i / NINTERM)) == i)	/* Outputs are done here */ {
      TimeToWrite = YES;
      SendOutput (TimeStep, gas_density, gas_v_rad, gas_v_theta, gas_label);
      WritePlanetSystemFile (sys, TimeStep);
      if ((OnlyInit) || ((GotoNextOutput) && (!StillWriteOneOutput))) {
	MPI_Finalize();
	return 0;
      }
      StillWriteOneOutput--;
      if (TimeInfo == YES)	/* Time monitoring is done here */
	GiveTimeInfo (TimeStep);
    }
    else {
      TimeToWrite = NO;
    }
				/* Algorithm loop begins here */

				/***********************/
				/* Hydrodynamical Part */
				/***********************/
    InitSpecificTime (Profiling, &t_Hydro, "Eulerian Hydro algorithms");
    AlgoGas (gas_density, gas_v_rad, gas_v_theta, gas_label, sys);
    GiveSpecificTime (Profiling, t_Hydro);
    SolveOrbits (sys);
    if (MonitorIntegral == YES) {
      masterprint ("Gas Momentum   : %.18g\n", GasMomentum (gas_density, gas_v_theta));
      masterprint ("Gas total Mass : %.18g\n", GasTotalMass (gas_density));
    }
  }
  FreePlanetary (sys);
  MPI_Finalize ();
  return 0;
}
