/** \file Pframeforce.c

Functions that evaluate the %force between the planets and the disk.
The FillForcesArrays() function is ill-named: it rather fill an array
of the potential, that is later use to derive the force acting on the
disk at every zone.  The name of this file is due to the fact that we
work in the frame centered on the primary (which is therefore not
inertial). Old versions of fargo also featured the file Gframeforce.c,
in which we worked in the frame centered on the center of gravity of
the system.  The present file also contains the functions necessary to
update the planets' position and velocities (taking into account, or
not, the indirect term, ie the term arising from the fact that the
frame is not inertial), as well as the function that initializes the
hydrodynamics fields with analytic prescription.

*/

#include "fargo.h"

extern boolean AllowAccretion, Corotating, Indirect_Term;
extern Pair DiskOnPrimaryAcceleration;
static Pair IndirectTerm;
static real q0[MAX1D], q1[MAX1D], PlanetMasses[MAX1D];
static real vt_int[MAX1D], vt_cent[MAX1D];

void ComputeIndirectTerm () {
  IndirectTerm.x = -DiskOnPrimaryAcceleration.x;
  IndirectTerm.y = -DiskOnPrimaryAcceleration.y; 
  if (Indirect_Term == NO) {
    IndirectTerm.x = 0.0;
    IndirectTerm.y = 0.0;
  }
}
				/* Below : work in non-rotating frame */
				/* centered on the primary */
void FillForcesArrays (sys)
PlanetarySystem *sys;
{
  int i,ii,j,l,nr,ns,k,NbPlanets;
  real x, y, angle, distance, distancesmooth, iplanet;
  real xplanet, yplanet, RRoche,smooth, mplanet, frac;
  real PlanetDistance, *Pot, pot, smoothing, cs;
  real InvPlanetDistance3, InvDistance;
  Pot= Potential->Field;
  nr = Potential->Nrad;
  ns = Potential->Nsec;
  NbPlanets = sys->nb;
  ComputeIndirectTerm();
#pragma omp parallel for
  for (i = 0; i < (nr+1)*ns; i++) Pot[i] = 0.0;
  for (k = 0; k < NbPlanets; k++) {
    xplanet = sys->x[k];
    yplanet = sys->y[k];
    mplanet = sys->mass[k]*MassTaper;
    PlanetDistance = sqrt(xplanet*xplanet+yplanet*yplanet);
    InvPlanetDistance3 =  1.0/PlanetDistance/PlanetDistance/PlanetDistance;
    RRoche = PlanetDistance*pow((1.0/3.0*mplanet),1.0/3.0);
    if (RocheSmoothing) {
      smoothing = RRoche*ROCHESMOOTHING;
    } else {
      iplanet = GetGlobalIFrac (PlanetDistance);
      frac = iplanet-floor(iplanet);
      ii = (int)iplanet;
      cs = GLOBAL_SOUNDSPEED[ii]*(1.0-frac)+\
	GLOBAL_SOUNDSPEED[ii+1]*frac;
      smoothing = cs * PlanetDistance * sqrt(PlanetDistance) * THICKNESSSMOOTHING;
    }
    smooth = smoothing*smoothing;
#pragma omp parallel for private(InvDistance,j,l,angle,x,y,distance,distancesmooth,pot)
    for (i = 0; i < nr; i++) {
      InvDistance = 1.0/Rmed[i];
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	angle = (real)j/(real)ns*2.0*PI;
	x = Rmed[i]*cos(angle);
	y = Rmed[i]*sin(angle);
	distance = (x-xplanet)*(x-xplanet)+(y-yplanet)*(y-yplanet);
	distancesmooth = sqrt(distance+smooth);
	pot = -G*mplanet/distancesmooth;
	if (Indirect_Term == YES)
	  pot += G*mplanet*InvPlanetDistance3*(x*xplanet+y*yplanet); /* Indirect term from planet  */
	Pot[l] += pot;
      }
    }
  }
#pragma omp parallel for private(InvDistance,j,l,angle,x,y,pot)
  for (i = 0; i < nr; i++) {
    InvDistance = 1.0/Rmed[i];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      angle = (real)j/(real)ns*2.0*PI;
      x = Rmed[i]*cos(angle);
      y = Rmed[i]*sin(angle);
      pot = -G*1.0*InvDistance;  /*    Central Mass is 1 */
      pot -= IndirectTerm.x*x + IndirectTerm.y*y; /* Indirect term from disk only */
      Pot[l] += pot;	
    }
  }
} 

void AdvanceSystemFromDisk (Rho, sys, dt)
PlanetarySystem *sys;
PolarGrid *Rho;
real dt;		       
{
  int NbPlanets, k, ii;
  Pair gamma;
  real x,y,r,m, iplanet, frac, cs, smoothing;
  NbPlanets = sys->nb;
  for (k = 0; k < NbPlanets; k++) {
    if (sys->FeelDisk[k] == YES) {
      m=sys->mass[k];
      x=sys->x[k];
      y=sys->y[k];
      r=sqrt(x*x+y*y);
      if (RocheSmoothing) {
	smoothing = r*pow(m/3.0,1./3.)*ROCHESMOOTHING;
      } else {
	iplanet = GetGlobalIFrac (r);
	frac = iplanet-floor(iplanet);
	ii = (int)iplanet;
	cs = GLOBAL_SOUNDSPEED[ii]*(1.0-frac)+\
	  GLOBAL_SOUNDSPEED[ii+1]*frac;
	smoothing = cs * r * sqrt(r) * THICKNESSSMOOTHING;
      }
      gamma = ComputeAccel (Rho, x, y, smoothing, m);
      sys->vx[k] += dt * gamma.x;
      sys->vy[k] += dt * gamma.y;
      sys->vx[k] += dt * IndirectTerm.x;
      sys->vy[k] += dt * IndirectTerm.y;
    }
  }
}

void AdvanceSystemRK5 (sys, dt)
PlanetarySystem *sys;
real dt;
{
  int i, n;
  boolean *feelothers;
  real dtheta, omega, rdot, x, y, r, new_r, vx, vy, theta, denom;
  n = sys->nb;
  for (i = 0; i < n; i++) {
    q0[i] = sys->x[i];
    q0[i+n] = sys->y[i];
    q0[i+2*n] = sys->vx[i];
    q0[i+3*n] = sys->vy[i];
    PlanetMasses[i] = sys->mass[i];
  }
  feelothers = sys->FeelOthers;
  RungeKunta (q0, dt, PlanetMasses, q1, n, feelothers);
  for (i = 1-(PhysicalTime >= RELEASEDATE); i < sys->nb; i++) {
    sys->x[i] = q1[i];
    sys->y[i] = q1[i+n];
    sys->vx[i] = q1[i+2*n];
    sys->vy[i] = q1[i+3*n];
  }
  if (PhysicalTime < RELEASEDATE) {
    x = sys->x[0];
    y = sys->y[0];
    r = sqrt(x*x+y*y);
    theta = atan2(y,x);
    rdot = (RELEASERADIUS-r)/(RELEASEDATE-PhysicalTime);
    omega = sqrt((1.+sys->mass[0])/r/r/r);
    new_r = r + rdot*dt;
    denom = r-new_r;
    if (denom != 0.0) {
      dtheta = 2.*dt*r*omega/denom*(sqrt(r/new_r)-1.);
    } else {
      dtheta = omega*dt;
    }
    vx = rdot;
    vy = new_r*sqrt((1.+sys->mass[0])/new_r/new_r/new_r);
    sys->x[0] = new_r*cos(dtheta+theta);
    sys->y[0] = new_r*sin(dtheta+theta);
    sys->vx[0]= vx*cos(dtheta+theta)-vy*sin(dtheta+theta); 
    sys->vy[0]= vx*sin(dtheta+theta)+vy*cos(dtheta+theta); 
  }
}

void SolveOrbits (sys)
PlanetarySystem *sys;
{
  int i, n;
  real x,y,vx,vy;
  n = sys->nb;
  for (i = 0; i < n; i++) {
    x = sys->x[i];
    y = sys->y[i];
    vx = sys->vx[i];
    vy = sys->vy[i];
    FindOrbitalElements (x,y,vx,vy,1.0+sys->mass[i],i);
  }
} 

real ConstructSequence (u, v, n)
     real *u, *v;
     int n;
{
  int i;
  real lapl=0.0;
  for (i = 1; i < n; i++)
    u[i] = 2.0*v[i]-u[i-1];
  for (i = 1; i < n-1; i++) {
    lapl += fabs(u[i+1]+u[i-1]-2.0*u[i]);
  }
  return lapl;
}

void
InitGas (Rho, Vr, Vt)
PolarGrid *Rho, *Vr, *Vt;
{
  int i, j, l, nr, ns;
  real *dens, *vr, *vt;
  float temporary;
  FILE *CS;
  char csfile[512];
  real  r, rg, omega, ri;
  real viscosity, t1, t2, r1, r2;
  dens= Rho->Field;
  vr  = Vr->Field;
  vt  = Vt->Field;
  nr  = Rho->Nrad;
  ns  = Rho->Nsec;
  sprintf (csfile, "%s%s", OUTPUTDIR, "soundspeed.dat");
  CS = fopen (csfile, "r");
  if (CS == NULL) {
    for (i = 0; i < nr; i++) {
      SOUNDSPEED[i] = AspectRatio(Rmed[i]) * sqrt(G*1.0/Rmed[i]) * pow(Rmed[i], FLARINGINDEX);
    }
    for (i = 0; i < GLOBALNRAD; i++) {
      GLOBAL_SOUNDSPEED[i] = AspectRatio(GlobalRmed[i]) * sqrt(G*1.0/GlobalRmed[i]) * pow(GlobalRmed[i], FLARINGINDEX);
    }
  } else {
    masterprint ("Reading soundspeed.dat file\n");
    for (i = 0; i < GLOBALNRAD; i++) {
      fscanf (CS, "%f", &temporary);
      GLOBAL_SOUNDSPEED[i] = (real)temporary;
    }
    for (i = 0; i < nr; i++) {
      SOUNDSPEED[i] = GLOBAL_SOUNDSPEED[i+IMIN];
    }
  }
  for (i = 1; i < GLOBALNRAD; i++) {
    vt_int[i]=(GLOBAL_SOUNDSPEED[i]*GLOBAL_SOUNDSPEED[i]*Sigma(GlobalRmed[i])-\
	       GLOBAL_SOUNDSPEED[i-1]*GLOBAL_SOUNDSPEED[i-1]*Sigma(GlobalRmed[i-1]))/\
      (.5*(Sigma(GlobalRmed[i])+Sigma(GlobalRmed[i-1])))/(GlobalRmed[i]-GlobalRmed[i-1])+\
      G*(1.0/GlobalRmed[i-1]-1.0/GlobalRmed[i])/(GlobalRmed[i]-GlobalRmed[i-1]);
    vt_int[i] = sqrt(vt_int[i]*Radii[i])-Radii[i]*OmegaFrame;
  }
  t1 = vt_cent[0] = vt_int[1]+.75*(vt_int[1]-vt_int[2]);
  r1 = ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
  vt_cent[0] += .25*(vt_int[1]-vt_int[2]);
  t2 = vt_cent[0];
  r2 = ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
  t1 = t1-r1/(r2-r1)*(t2-t1);
  vt_cent[0] = t1;
  ConstructSequence (vt_cent, vt_int, GLOBALNRAD);
  vt_cent[GLOBALNRAD]=vt_cent[GLOBALNRAD-1];
  for (i = 0; i <= nr; i++) {
    if (i == nr) {
      r = Rmed[nr-1];
      ri= Rinf[nr-1];
    }
    else {
      r = Rmed[i];
      ri= Rinf[i];
    }
    viscosity = FViscosity (r);
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      rg = r;
      omega = sqrt(G*1.0/rg/rg/rg);
      vt[l] = omega*r*\
	sqrt(1.0-pow(ASPECTRATIO,2.0)*\
	     pow(r,2.0*FLARINGINDEX)*\
	     (1.+SIGMASLOPE-2.0*FLARINGINDEX));
      vt[l] -= OmegaFrame*r;
      if (CentrifugalBalance)
	vt[l] = vt_cent[i+IMIN];
      if (i == nr) 
	vr[l] = 0.0;
      else {
	vr[l] = IMPOSEDDISKDRIFT*SIGMA0/SigmaInf[i]/ri;
	if (ViscosityAlpha) {
	  vr[l] -= 3.0*viscosity/r*(-SIGMASLOPE+2.0*FLARINGINDEX+1.0);
	} else {
	  vr[l] -= 3.0*viscosity/r*(-SIGMASLOPE+.5);
	}
      }
      dens[l] = SigmaMed[i];
    }
  }
  for (j = 0; j < ns; j++)
    vr[j] = vr[j+ns*nr] = 0.0;
}
