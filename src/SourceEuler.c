/** \file SourceEuler.c 

Contains routines used by the hydrodynamical loop. More specifically,
it contains the main loop itself and all the source term substeps
(with the exception of the evaluation of the viscous force). The
transport substep is treated elsewhere. */

#include "fargo.h"

#define CFLSECURITY 0.5		/* Maximum fraction of zone size */
				/* swept in one timestep */

#define CVNR 1.41       	/* Shocks are spread over CVNR zones:       */
                                /* von Neumann-Richtmyer viscosity constant */
				/* Beware of misprint in Stone and Norman's */
				/* paper : use C2^2 instead of C2           */

static PolarGrid *TemperInt;
static PolarGrid *VradNew,   *VradInt;
static PolarGrid *VthetaNew, *VthetaInt;
static real timeCRASH;  
extern boolean Corotating;

static int AlreadyCrashed = 0, GasTimeStepsCFL;

extern int TimeStep;
extern boolean FastTransport, IsDisk;
Pair DiskOnPrimaryAcceleration;


boolean DetectCrash (array)
PolarGrid *array;
{
  int i, j, l, nr, ns;
  real *ptr;
  boolean bool = NO;
  nr = array->Nrad;
  ns = array->Nsec;
  ptr= array->Field;
#pragma omp parallel for private(j,l) shared(bool)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (ptr[l] < 0.0) 
	bool = YES;
    }
  }
  return bool;
}
 
void FillPolar1DArrays ()
{
  FILE *input, *output;
  int i,ii;
  real drrsep;
  float temporary;
  char InputName[256], OutputName[256];
  drrsep = (RMAX-RMIN)/(real)GLOBALNRAD;
  sprintf (InputName, "%s%s", OUTPUTDIR, "radii.dat");
  sprintf (OutputName, "%s%s", OUTPUTDIR, "used_rad.dat");
  input = fopen (InputName, "r");
  if (input == NULL) {
    mastererr ("Warning : no `radii.dat' file found. Using default.\n");
    if (LogGrid == YES) {
      for (i = 0; i <= GLOBALNRAD; i++) {
	Radii[i] = RMIN*exp((real)i/(real)GLOBALNRAD*log(RMAX/RMIN));
      }
    } else {
      for (i = 0; i <= GLOBALNRAD; i++) {
	Radii[i] = RMIN+drrsep*(real)(i);
      }
    }
  } else {
    mastererr ("Reading 'radii.dat' file.\n");
    for (i = 0; i <= GLOBALNRAD; i++) {
      fscanf (input, "%f", &temporary);
      Radii[i] = (real)temporary;
    }
  }
  for (i = 0; i < GLOBALNRAD; i++) {
    GlobalRmed[i] = 2.0/3.0*(Radii[i+1]*Radii[i+1]*Radii[i+1]-Radii[i]*Radii[i]*Radii[i]);
    GlobalRmed[i] = GlobalRmed[i] / (Radii[i+1]*Radii[i+1]-Radii[i]*Radii[i]);
  }
  for (i = 0; i < NRAD; i++) {
    ii = i+IMIN;
    Rinf[i] = Radii[ii];
    Rsup[i] = Radii[ii+1];
    Rmed[i] = 2.0/3.0*(Rsup[i]*Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]*Rinf[i]);
    Rmed[i] = Rmed[i] / (Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]);
    Surf[i] = PI*(Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i])/(real)NSEC;
    InvRmed[i] = 1.0/Rmed[i];
    InvSurf[i] = 1.0/Surf[i];
    InvDiffRsup[i] = 1.0/(Rsup[i]-Rinf[i]);
    InvRinf[i] = 1.0/Rinf[i];
  }
  Rinf[NRAD]=Radii[NRAD+IMIN];
  for (i = 1; i < NRAD; i++) {
    InvDiffRmed[i] = 1.0/(Rmed[i]-Rmed[i-1]);
  }
  if (CPU_Master) {
    output = fopen (OutputName, "w");
    if (output == NULL) {
      mastererr ("Can't write %s.\nProgram stopped.\n", OutputName);
      prs_exit (1);
    }
    for (i = 0; i <= GLOBALNRAD; i++) {
      fprintf (output, "%.18g\n", Radii[i]);
    }
    fclose (output);
  }
  if (input != NULL) fclose (input);
}

void InitEuler (Rho, Vr, Vt)
PolarGrid *Rho, *Vr, *Vt;
{
  FillPolar1DArrays ();
  FillSigma ();
  InitTransport ();
  InitViscosity ();
  RhoStar      = CreatePolarGrid(NRAD, NSEC, "RhoStar");
  RhoInt       = CreatePolarGrid(NRAD, NSEC, "RhoInt");
  VradNew      = CreatePolarGrid(NRAD, NSEC, "VradNew");
  VradInt      = CreatePolarGrid(NRAD, NSEC, "VradInt");
  VthetaNew    = CreatePolarGrid(NRAD, NSEC, "VthetaNew");
  VthetaInt    = CreatePolarGrid(NRAD, NSEC, "VthetaInt");
  TemperInt    = CreatePolarGrid(NRAD, NSEC, "TemperInt");
  Potential    = CreatePolarGrid(NRAD, NSEC, "Potential");
  InitGas (Rho, Vr, Vt);
}

real min2 (a,b)
real a,b;
{
  if (b < a) return b;
  return a;
}

real max2 (a,b)
real a,b;
{
  if (b > a) return b;
  return a;
}


void ActualiseGas (array, newarray)
PolarGrid *array, *newarray;
{
  int i,j,l,ns,nr;
  real *old, *new;
  nr = array->Nrad;
  ns = array->Nsec;
  old= array->Field;
  new= newarray->Field;
#pragma omp parallel for private(j,l)
  for (i = 0; i <= nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+ns*i;
      old[l] = new[l];
    }
  }
}


void AlgoGas (Rho,Vrad,Vtheta,Label,sys)
PolarGrid *Rho, *Vrad, *Vtheta, *Label;
PlanetarySystem *sys;
{
  real dt, dtemp=0.0;
  real OmegaNew, domega;
  int gastimestepcfl;
  boolean Crashed=NO;
  gastimestepcfl = 1;
  if (IsDisk == YES) {
    CommunicateBoundaries (Rho,Vrad,Vtheta,Label);
    if (SloppyCFL == YES) 
      gastimestepcfl = ConditionCFL (Vrad, Vtheta, DT-dtemp);
  }
  MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  dt = DT / (real)GasTimeStepsCFL;
  while (dtemp < 0.999999999*DT) {
    MassTaper = PhysicalTime/(MASSTAPER*2.0*M_PI);
    MassTaper = (MassTaper > 1.0 ? 1.0 : pow(sin(MassTaper*M_PI/2.0),2.0));
    if (IsDisk == YES) {
      CommunicateBoundaries (Rho,Vrad,Vtheta,Label);
      if (SloppyCFL == NO) {
	gastimestepcfl = 1;
	gastimestepcfl = ConditionCFL (Vrad, Vtheta, DT-dtemp);
	MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	dt = (DT-dtemp)/(real)GasTimeStepsCFL;
      }
      AccreteOntoPlanets (Rho, Vrad, Vtheta, dt, sys);
    }
    dtemp += dt;
    DiskOnPrimaryAcceleration.x = 0.0;
    DiskOnPrimaryAcceleration.y = 0.0;
    if (Corotating == YES) GetPsysInfo (sys, MARK);
    if (IsDisk == YES) {
      DiskOnPrimaryAcceleration   = ComputeAccel (Rho, 0.0, 0.0, 0.0, 0.0);
      FillForcesArrays (sys);
      AdvanceSystemFromDisk (Rho, sys, dt);
    }
    AdvanceSystemRK5 (sys, dt);
    if (Corotating == YES) {
      OmegaNew = GetPsysInfo(sys, GET) / dt;
      domega = OmegaNew-OmegaFrame;
      if (IsDisk == YES)
	CorrectVtheta (Vtheta, domega);
      OmegaFrame = OmegaNew;
    }
    RotatePsys (sys, OmegaFrame*dt);
    if (IsDisk == YES) {
      ApplyBoundaryCondition(Vrad, Vtheta, Rho, dt);
      Crashed = DetectCrash (Rho);  /* test for negative density values */
      if (Crashed == YES) {
	if (AlreadyCrashed == 0) {
	  timeCRASH=PhysicalTime;   /* if it appears to be the first crash */
	  fprintf (stdout,"\nCrash! at time %.12g\n", timeCRASH);
	  WriteDiskPolar (Rho, 999);    /* We write the HD arrays */
	  WriteDiskPolar (Vrad, 999);   /* in order to keep a track */
	  WriteDiskPolar (Vtheta, 999); /* of what happened */
	}
	AlreadyCrashed++;
	masterprint ("c");
      } else {
	masterprint (".");
      }
      fflush (stdout);
      SubStep1 (Vrad, Vtheta, Rho, dt);
      SubStep2 (Rho, dt);
      ActualiseGas (Vrad, VradNew);
      ActualiseGas (Vtheta, VthetaNew);
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, dt);
      Transport (Rho, Vrad, Vtheta, Label, dt);
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, dt);
    }
    PhysicalTime += dt;
  }
  masterprint ("\n");
}

void SubStep1 (Vrad, Vtheta, Rho, dt)
PolarGrid *Vrad, *Vtheta, *Rho;
real dt;
{
  int i, j, l, lim, ljm, ljp, nr, ns;
  real *vrad, *vtheta, *rho;
  real *vradint, *vthetaint;
  real gradp, gradphi, vt2, dxtheta, cs2, cs2m;
  real invdxtheta;
  real supp_torque=0.0;		/* for imposed disk drift */
  real *Pot;
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  rho = Rho->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  vradint = VradInt->Field;
  vthetaint = VthetaInt->Field;
  Pot = Potential->Field;
				/* In this substep we take into account     */
				/* the source part of Euler equations       */
				/* (i.e. the R.H.S. in classical notation). */
#pragma omp parallel private(j,l,lim,ljm,ljp,cs2,cs2m,dxtheta,vt2,gradp,gradphi,invdxtheta,supp_torque)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	cs2 = SOUNDSPEED[i]*SOUNDSPEED[i];
	cs2m= SOUNDSPEED[i-1]*SOUNDSPEED[i-1];
	gradp = (cs2*rho[l]-cs2m*rho[lim])*2.0/(rho[l]+rho[lim])*InvDiffRmed[i];
	gradphi = (Pot[l]-Pot[lim])*InvDiffRmed[i];
	vt2 = vtheta[l]+vtheta[ljp]+vtheta[lim]+vtheta[ljp-ns];
	vt2 = vt2/4.0+Rinf[i]*OmegaFrame;
	vt2 = vt2*vt2;
	vradint[l] = vrad[l]+dt*(-gradp-gradphi+vt2*InvRinf[i]);
      }
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      supp_torque = IMPOSEDDISKDRIFT*.5*pow(Rmed[i],-2.5+SIGMASLOPE);
      dxtheta = 2.0*PI/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	cs2m = cs2 = SOUNDSPEED[i]*SOUNDSPEED[i];
	gradp = cs2*(rho[l]-rho[ljm])*2.0/(rho[l]+rho[ljm])*invdxtheta;
	gradphi = (Pot[l]-Pot[ljm])*invdxtheta;   
	vthetaint[l] = vtheta[l]-dt*(gradp+gradphi);
	vthetaint[l] += dt*supp_torque;
      }
    }
  }
  ViscousTerms (VradInt, VthetaInt, Rho, dt);
}

void SubStep2 (Rho, dt)
PolarGrid *Rho;
real dt;
{
  int i, j, l, lim, lip, ljm, ljp, nr, ns;
  real *vrad, *vtheta, *rho;
  real *vradnew, *vthetanew, *qt, *qr;
  real dxtheta, invdxtheta;
  real dv;
  nr = VradInt->Nrad;
  ns = VradInt->Nsec;
  rho = Rho->Field;
  vrad = VradInt->Field;
  vtheta = VthetaInt->Field;
  qr = RhoInt->Field;
  vradnew = VradNew->Field;
  vthetanew = VthetaNew->Field;
  qt = TemperInt->Field;
#pragma omp parallel for private(j,dxtheta,l,lim,lip,ljm,ljp,dv)
  for (i = 0; i < nr; i++) {
    dxtheta = 2.0*PI/(real)ns*Rmed[i];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lim = l-ns;
      lip = l+ns;
      ljm = l-1;
      if (j == 0) ljm = i*ns+ns-1;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      dv = vrad[lip]-vrad[l];
      if (dv < 0.0)
        qr[l] = CVNR*CVNR*rho[l]*dv*dv;
      else 
        qr[l] = 0.0; 
      dv = vtheta[ljp]-vtheta[l];
      if (dv < 0.0)
        qt[l] = CVNR*CVNR*rho[l]*dv*dv;
      else
	qt[l] = 0.0;
    }
  }
#pragma omp parallel private(l,lim,lip,ljm,ljp,j,dxtheta,invdxtheta)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	lip = l+ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	vradnew[l] = vrad[l]-dt*2.0/(rho[l]+rho[lim])*(qr[l]-qr[lim])*InvDiffRmed[i];
      }
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      dxtheta = 2.0*PI/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	lip = l+ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	vthetanew[l] = vtheta[l]-dt*2.0/(rho[l]+rho[ljm])*(qt[l]-qt[ljm])*invdxtheta;
      }
    }
  }
}
		   		   
int ConditionCFL (Vrad, Vtheta, deltaT)
PolarGrid *Vrad, *Vtheta;
real deltaT;
{
  static real Vresidual[MAX1D], Vmoy[MAX1D];
  int i, j, l, ns, nr, lip, ljp;
  real invdt1, invdt2, invdt3, invdt4, invdt5=1e-30, cs, newdt, dt;
  int ideb=0, jdeb=0;
  real itdbg1=0., itdbg2=0., itdbg3=0., itdbg4=0., itdbg5=0., mdtdbg=0.;
  /* debugging variables */
  real *vt, *vr, dxrad, dxtheta, dvr, dvt, viscr=0., visct=0.;
  ns = Vtheta->Nsec;
  nr = Vtheta->Nrad;
  vt = Vtheta->Field;
  vr = Vrad->Field;
  newdt = 1e30;
  for (i = 0; i < nr; i++) {
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*2.0*PI/(real)ns;
    Vmoy[i] = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      Vmoy[i] += vt[l];
    }
    Vmoy[i] /= (real)ns;
  }
  for (i = One_or_active; i < Max_or_active; i++) {
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*2.0*PI/(real)ns;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (FastTransport == YES)
	Vresidual[j] = vt[l]-Vmoy[i];  /* FARGO algorithm */
      else
	Vresidual[j] = vt[l];	       /* Standard algorithm */
    }
    Vresidual[ns]=Vresidual[0];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      cs = SOUNDSPEED[i];
      invdt1 = cs/(min2(dxrad,dxtheta));
      invdt2 = fabs(vr[l])/dxrad;
      invdt3 = fabs(Vresidual[j])/dxtheta;
      dvr = vr[lip]-vr[l];
      dvt = vt[ljp]-vt[l];
      if (dvr >= 0.0) dvr = 1e-10;
      else dvr = -dvr;
      if (dvt >= 0.0) dvt = 1e-10;
      else dvt = -dvt;
      invdt4 = max2(dvr/dxrad,dvt/dxtheta);
      invdt4*= 4.0*CVNR*CVNR;
      if ( ViscosityAlpha || (VISCOSITY != 0.0) )
	invdt5 = FViscosity(Rmed[i])*4./pow(min2(dxrad,dxtheta),2);
      dt = CFLSECURITY/sqrt(invdt1*invdt1+invdt2*invdt2+invdt3*invdt3+\
			    invdt4*invdt4+invdt5*invdt5);
      if (dt < newdt) {
	newdt = dt;
	if (debug == YES) {
	  ideb = i;
	  jdeb = j;
	  itdbg1 = 1.0/invdt1; itdbg2=1.0/invdt2;
	  itdbg3=1.0/invdt3; itdbg4=1.0/invdt4;
	  itdbg5=1.0/invdt5;
	  mdtdbg = newdt;
	  viscr = dxrad/dvr/4.0/CVNR/CVNR;     
	  visct = dxtheta/dvt/4.0/CVNR/CVNR;
	}
      }  
    }
  }
  for (i = Zero_or_active; i < MaxMO_or_active; i++) {
    dt = 2.0*PI*CFLSECURITY/(real)NSEC/fabs(Vmoy[i]*InvRmed[i]-Vmoy[i+1]*InvRmed[i+1]);
    if (dt < newdt) newdt = dt;
  }
  if (deltaT < newdt) newdt = deltaT;
  if (debug == YES) {
    printf ("Timestep control information for CPU %d: \n", CPU_Rank);
    printf ("Most restrictive cell at i=%d and j=%d\n", ideb, jdeb);
    printf ("located at radius Rmed         : %g\n", Rmed[ideb]);
    printf ("Sound speed limit              : %g\n", itdbg1);
    printf ("Radial motion limit            : %g\n", itdbg2);
    printf ("Residual circular motion limit : %g\n", itdbg3);
    printf ("Artificial viscosity limit     : %g\n", itdbg4);
    printf ("   Arise from r with limit     : %g\n", viscr);
    printf ("   and from theta with limit   : %g\n", visct);
    printf ("Physical viscosity limit       : %g\n", itdbg5);
    printf ("Limit time step for this cell  : %g\n", mdtdbg);
    printf ("Limit time step adopted        : %g\n", newdt);
    if (newdt < mdtdbg) {
      printf ("Discrepancy arise either from shear.\n");
      printf ("or from the imposed DT interval.\n");
    }
  }
  return (int)(ceil(deltaT/newdt));
}
