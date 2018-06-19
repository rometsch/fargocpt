/** \file Force.c

Contains the function used to evaluate the %force due to disk, and the
function that writes the 'tqwk' log files. Although the planet mass is
given as an argument to ComputeForce(), this mass is used only to
specify the distance cutoff in the case of the Hill sphere avoidance.
The force returned is a specific force. It has therefore the dimension
of an acceleration (LT^-2).
*/

#include "fargo.h"

extern boolean OpenInner, NonReflecting;
extern Pair DiskOnPrimaryAcceleration;

Force ComputeForce (Rho, x, y, rsmoothing, mass)
PolarGrid *Rho;
real x, y, rsmoothing, mass;
{
  int i, j, l, ns;
  real localforce[8]={0.,0.,0.,0.,0.,0.,0.,0.}, globalforce[8]={0.,0.,0.,0.,0.,0.,0.,0.};
  real xc, yc, cellmass, dx, dy, distance, dist2, rh, a;
  real InvDist3, fxi, fyi, fxhi, fyhi, fxo, fyo, fxho, fyho, hill_cut;
  real *dens, *abs, *ord;
  Force Force;
  ns = Rho->Nsec;
  dens = Rho->Field;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  fxi = fyi = fxhi = fyhi = fxo = fyo = fxho = fyho = 0.0;
  a = sqrt(x*x+y*y);
  rh = pow(mass/3.0, 1./3.)*a+1e-15;
  if (FakeSequential && (CPU_Rank > 0)) {
    MPI_Recv (&globalforce, 8, MPI_DOUBLE, CPU_Rank-1, 27, MPI_COMM_WORLD, &fargostat);
    fxi = globalforce[0];
    fyi = globalforce[1];
    fxhi= globalforce[2];
    fyhi= globalforce[3];
    fxo = globalforce[4];
    fyo = globalforce[5];
    fxho= globalforce[6];
    fyho= globalforce[7];
  }
#pragma omp parallel for private(j,hill_cut,cellmass,l,xc,yc,dist2,distance,InvDist3,dx,dy) shared(fxi,fyi,fxhi,fyhi,fxo,fyo,fxho,fyho)
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      xc = abs[l];
      yc = ord[l];
      cellmass = Surf[i]*dens[l];
      dx = xc-x;
      dy = yc-y;
      dist2 = dx*dx+dy*dy;
      hill_cut = 1.-exp(-dist2/(rh*rh));
      dist2 += rsmoothing*rsmoothing;
      distance = sqrt(dist2);
      InvDist3 = 1.0/dist2/distance;
      if (Rmed[i] < a) {
#pragma omp atomic
	fxi += G*cellmass*dx*InvDist3;
#pragma omp atomic
	fyi += G*cellmass*dy*InvDist3;
#pragma omp atomic
	fxhi+= G*cellmass*dx*InvDist3*hill_cut;
#pragma omp atomic
	fyhi+= G*cellmass*dy*InvDist3*hill_cut;
      } else {
#pragma omp atomic
	fxo += G*cellmass*dx*InvDist3;
#pragma omp atomic
	fyo += G*cellmass*dy*InvDist3;
#pragma omp atomic
	fxho+= G*cellmass*dx*InvDist3*hill_cut;
#pragma omp atomic
	fyho+= G*cellmass*dy*InvDist3*hill_cut;
      }
    }
  }
  if (FakeSequential) {
    globalforce[0]=fxi;
    globalforce[1]=fyi;
    globalforce[2]=fxhi;
    globalforce[3]=fyhi;
    globalforce[4]=fxo;
    globalforce[5]=fyo;
    globalforce[6]=fxho;
    globalforce[7]=fyho;
    if (CPU_Rank < CPU_Number-1)
      MPI_Send (&globalforce, 8, MPI_DOUBLE, CPU_Rank+1, 27, MPI_COMM_WORLD);
  } else {
    localforce[0] = fxi;
    localforce[1] = fyi;
    localforce[2] = fxhi;
    localforce[3] = fyhi;
    localforce[4] = fxo;
    localforce[5] = fyo;
    localforce[6] = fxho;
    localforce[7] = fyho;
    MPI_Allreduce (&localforce, &globalforce, 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  if (FakeSequential)
    MPI_Bcast (&globalforce, 8, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
  Force.fx_inner    = globalforce[0];
  Force.fy_inner    = globalforce[1];
  Force.fx_ex_inner = globalforce[2];
  Force.fy_ex_inner = globalforce[3];
  Force.fx_outer    = globalforce[4];
  Force.fy_outer    = globalforce[5];
  Force.fx_ex_outer = globalforce[6];
  Force.fy_ex_outer = globalforce[7];
  return Force;
}

void UpdateLog (psys, Rho, outputnb, time)
     PolarGrid *Rho;
     PlanetarySystem *psys;
     int outputnb;
     real time;
{
  int i, ii, nb;
  real x, y, r, m, vx, vy, smoothing, iplanet, cs, frac;
  Force fc;
  FILE *out;
  char filename[MAX1D];
  nb = psys->nb;
  for (i = 0; i < nb; i++) {
    x = psys->x[i];
    y = psys->y[i];
    vx = psys->vx[i];
    vy = psys->vy[i];
    r = sqrt(x*x+y*y);
    m = psys->mass[i];
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
    fc = ComputeForce (Rho, x, y, smoothing, m);
    if (CPU_Rank == CPU_Number-1) {
      sprintf (filename, "%stqwk%d.dat", OUTPUTDIR, i);
      out = fopenp (filename, "a");
      fprintf (out, "%d\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\t%.18g\n", outputnb,\
	       x*fc.fy_inner-y*fc.fx_inner,\
	       x*fc.fy_outer-y*fc.fx_outer,\
	       x*fc.fy_ex_inner-y*fc.fx_ex_inner,\
	       x*fc.fy_ex_outer-y*fc.fx_ex_outer,\
	       vx*fc.fx_inner+vy*fc.fy_inner,\
	       vx*fc.fx_outer+vy*fc.fy_outer,\
	       vx*fc.fx_ex_inner+vy*fc.fy_ex_inner,\
	       vx*fc.fx_ex_outer+vy*fc.fy_ex_outer, time);
      fclose (out);
    }
  }
}
