/** \file Planet.c

Accretion of disk material onto the planets, and solver of planetary
orbital elements.  The prescription used for the accretion is the one
designed by W. Kley.

*/

#include "fargo.h"

void AccreteOntoPlanets (Rho, Vrad, Vtheta, dt, sys)
real dt;
PolarGrid *Rho, *Vrad, *Vtheta;
PlanetarySystem *sys;
{
  real RRoche, Rplanet, distance, dx, dy, deltaM, angle, temp;
  int i_min,i_max, j_min, j_max, i, j, l, jf, ns, nr, lip, ljp, k;
  real Xplanet, Yplanet, Mplanet, VXplanet, VYplanet;
  real facc, facc1, facc2, frac1, frac2; /* We adopt the same notations as W. Kley */
  real *dens, *abs, *ord, *vrad, *vtheta;
  real PxPlanet, PyPlanet, vrcell, vtcell, vxcell, vycell, xc, yc;
  real dPxPlanet, dPyPlanet, dMplanet;
  nr     = Rho->Nrad;
  ns     = Rho->Nsec;
  dens   = Rho->Field;
  abs    = CellAbscissa->Field;
  ord    = CellOrdinate->Field;
  vrad   = Vrad->Field;
  vtheta = Vtheta->Field;
  for (k=0; k < sys->nb; k++) {
    if (sys->acc[k] > 1e-10) {
      dMplanet = dPxPlanet = dPyPlanet = 0.0;
				/* Hereafter : initialization of W. Kley's parameters */
      facc = dt*(sys->acc[k]);
      facc1 = 1.0/3.0*facc;
      facc2 = 2.0/3.0*facc;
      frac1 = 0.75;
      frac2 = 0.45;
				/* W. Kley's parameters initialization finished */
      Xplanet = sys->x[k];
      Yplanet = sys->y[k];
      VXplanet = sys->vx[k];
      VYplanet = sys->vy[k];
      Mplanet = sys->mass[k];
      Rplanet = sqrt(Xplanet*Xplanet+Yplanet*Yplanet);
      RRoche = pow((1.0/3.0*Mplanet),(1.0/3.0))*Rplanet; /* Central mass is 1.0 */
      i_min=0;
      i_max=nr-1;
      while ((Rsup[i_min] < Rplanet-RRoche) && (i_min < nr)) i_min++;
      while ((Rinf[i_max] > Rplanet+RRoche) && (i_max > 0)) i_max--;
      angle = atan2 (Yplanet, Xplanet);
      j_min =(int)((real)ns/2.0/PI*(angle - 2.0*RRoche/Rplanet));
      j_max =(int)((real)ns/2.0/PI*(angle + 2.0*RRoche/Rplanet));
      PxPlanet = Mplanet*VXplanet;
      PyPlanet = Mplanet*VYplanet;
#pragma omp parallel for private(j,jf,vrcell,vtcell,vxcell,vycell,l,lip,ljp,xc,yc,dx,dy,distance,deltaM) shared(dPxPlanet, dPyPlanet, dMplanet)
      for (i = i_min; i <= i_max; i++) {
	for (j = j_min; j <= j_max; j++) {
	  jf = j;
	  while (jf <  0)  jf += ns;
	  while (jf >= ns) jf -= ns;
	  l   = jf+i*ns;
	  lip = l+ns;
	  ljp = l+1;
	  if (jf == ns-1) ljp = i*ns;
	  xc = abs[l];
	  yc = ord[l];
	  dx = Xplanet-xc;
	  dy = Yplanet-yc;
	  distance = sqrt(dx*dx+dy*dy);
	  vtcell=0.5*(vtheta[l]+vtheta[ljp])+Rmed[i]*OmegaFrame;
	  vrcell=0.5*(vrad[l]+vrad[lip]);
	  vxcell=(vrcell*xc-vtcell*yc)/Rmed[i];
	  vycell=(vrcell*yc+vtcell*xc)/Rmed[i];
	  if (distance < frac1*RRoche) {
	    deltaM = facc1*dens[l]*Surf[i];
	    if (i < Zero_or_active) deltaM = 0.0;
	    if (i >= Max_or_active) deltaM = 0.0;
	    dens[l] *= (1.0 - facc1);
#pragma omp atomic
	    dPxPlanet    += deltaM*vxcell;
#pragma omp atomic
	    dPyPlanet    += deltaM*vycell;
#pragma omp atomic
	    dMplanet     += deltaM;
	  }
	  if (distance < frac2*RRoche) {
	    deltaM = facc2*dens[l]*Surf[i];
	    if (i < Zero_or_active) deltaM = 0.0;
	    if (i >= Max_or_active) deltaM = 0.0;
	    dens[l] *= (1.0 - facc2);
#pragma omp atomic
	    dPxPlanet    += deltaM*vxcell;
#pragma omp atomic
	    dPyPlanet    += deltaM*vycell;
#pragma omp atomic
	    dMplanet     += deltaM;
	  }
	}
      }
      MPI_Allreduce (&dMplanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      dMplanet = temp;
      MPI_Allreduce (&dPxPlanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      dPxPlanet = temp;
      MPI_Allreduce (&dPyPlanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      dPyPlanet = temp;
      PxPlanet += dPxPlanet;
      PyPlanet += dPyPlanet;
      Mplanet  += dMplanet;
      if (sys->FeelDisk[k] == YES) {
	sys->vx[k] = PxPlanet/Mplanet;
	sys->vy[k] = PyPlanet/Mplanet;
      }
      sys->mass[k] = Mplanet;
    }
  }
}


void FindOrbitalElements (x,y,vx,vy,m,n)
real x,y,vx,vy,m;
int n;
{
  real Ax, Ay, e, d, h, a, E, M, V;
  real PerihelionPA;
  FILE *output;
  char name[256];
  if (CPU_Rank != CPU_Number-1) return;
  sprintf (name, "%sorbit%d.dat", OUTPUTDIR, n);
  output = fopenp (name, "a");
  h = x*vy-y*vx;
  d = sqrt(x*x+y*y);
  Ax = x*vy*vy-y*vx*vy-G*m*x/d;
  Ay = y*vx*vx-x*vx*vy-G*m*y/d;
  e = sqrt(Ax*Ax+Ay*Ay)/G/m;
  a = h*h/G/m/(1-e*e);
  if (e != 0.0) {
    E = acos((1.0-d/a)/e);
  } else {
    E = 0.0;
  }
  if ((x*y*(vy*vy-vx*vx)+vx*vy*(x*x-y*y)) < 0) E= -E;
  M = E-e*sin(E);
  if (e != 0.0) {
    V = acos ((a*(1.0-e*e)/d-1.0)/e);
  } else {
    V = 0.0;
  }
  if (E < 0.0) V = -V;
  if (e != 0.0) {
    PerihelionPA=atan2(Ay,Ax);
  } else {
    PerihelionPA=atan2(y,x);
  }
  fprintf (output, "%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n", PhysicalTime, e, a, M, V,\
	   PerihelionPA);
  fclose (output);
}
 
