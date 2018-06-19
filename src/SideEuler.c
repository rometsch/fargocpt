/** \file SideEuler.c

Total mass and angular momentum monitoring, and boundary conditions.
In addition, this file contains a few low-level functions that
manipulate PolarGrid 's or initialize the forces evaluation.

*/

#include "fargo.h"

extern boolean OpenInner, NonReflecting, OuterSourceMass;
extern Pair DiskOnPrimaryAcceleration;

real GasTotalMass (array)
PolarGrid *array;
{
   int i,j,ns;
   real *density, total = 0.0, fulltotal=0.0;
   ns = array->Nsec;
   density = array->Field;
   if (FakeSequential && (CPU_Rank > 0)) 
     MPI_Recv (&total, 1, MPI_DOUBLE, CPU_Rank-1, 0, MPI_COMM_WORLD, &fargostat);
   for (i = Zero_or_active; i < Max_or_active; i++) {
     for (j = 0; j < ns; j++) {
       total += Surf[i]*density[j+i*ns];
     }
   }
   if (FakeSequential) {
     if (CPU_Rank < CPU_Number-1)
       MPI_Send (&total, 1, MPI_DOUBLE, CPU_Rank+1, 0, MPI_COMM_WORLD);
   }
   else
     MPI_Allreduce (&total, &fulltotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if (FakeSequential) {
     MPI_Bcast (&total, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
     fulltotal = total;
   }
   return fulltotal;
}

real GasMomentum (Density, Vtheta)
PolarGrid *Density, *Vtheta;
{
   int i,j,ns;
   real *density, *vtheta, total = 0.0, fulltotal=0.0;
   ns = Density->Nsec;
   density = Density->Field;
   vtheta = Vtheta->Field;
   if (FakeSequential && (CPU_Rank > 0)) 
     MPI_Recv (&total, 1, MPI_DOUBLE, CPU_Rank-1, 2, MPI_COMM_WORLD, &fargostat);
   for (i = Zero_or_active; i < Max_or_active; i++) {
     for (j = 1; j < ns; j++) {
       total += Surf[i]*(density[j+i*ns]+density[j-1+i*ns])*Rmed[i]*(vtheta[j+i*ns]+OmegaFrame*Rmed[i]);
     }
     total += Surf[i]*(density[i*ns]+density[i*ns+ns-1])*Rmed[i]*(vtheta[i*ns]+OmegaFrame*Rmed[i]);
   }
   if (FakeSequential) {
     if (CPU_Rank < CPU_Number-1)
       MPI_Send (&total, 1, MPI_DOUBLE, CPU_Rank+1, 2, MPI_COMM_WORLD);
   }
   else
     MPI_Allreduce (&total, &fulltotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if (FakeSequential) {
     MPI_Bcast (&total, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
     fulltotal = total;
   }
   return 0.5*fulltotal;
}

void DivisePolarGrid (Num, Denom, Res)
PolarGrid *Num, *Denom, *Res;
{
  int i,j,l,nr,ns;
  real *num, *denom, *res;
  num = Num->Field;
  denom=Denom->Field;
  res = Res->Field;
  ns = Res->Nrad;
  nr = Res->Nsec;
#pragma omp parallel for private(j,l)
  for (i = 0; i <= nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+ns*i;
      res[l] = num[l]/(denom[l]+1e-20);
    }
  }
}

void InitComputeAccel ()
{
  int i, j, l, nr, ns;
  real *abs, *ord;
  CellAbscissa = CreatePolarGrid (NRAD,NSEC,"abscissa");
  CellOrdinate = CreatePolarGrid (NRAD,NSEC,"ordinate");
  nr = CellAbscissa->Nrad;
  ns = CellAbscissa->Nsec;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
       l = j+i*ns;
       abs[l] = Rmed[i] * cos(2.0*PI*(real)j/(real)ns);
       ord[l] = Rmed[i] * sin(2.0*PI*(real)j/(real)ns);
    }
  }
}
  
Pair ComputeAccel (Rho, x, y, rsmoothing, mass)
PolarGrid *Rho;
real x, y, rsmoothing, mass;
{
  Pair acceleration;
  Force force;
  force = ComputeForce (Rho, x, y, rsmoothing, mass);
  if (ExcludeHill) {
    acceleration.x = force.fx_ex_inner+force.fx_ex_outer;
    acceleration.y = force.fy_ex_inner+force.fy_ex_outer;
  } else {
    acceleration.x = force.fx_inner+force.fx_outer;
    acceleration.y = force.fy_inner+force.fy_outer;
  }
  return acceleration;
}

void OpenBoundary (Vrad, Rho)
PolarGrid *Vrad, *Rho;
{
  int i,j,l,ns;
  real *rho, *vr;
  if (CPU_Rank != 0) return;
  ns = Rho->Nsec;
  rho = Rho->Field;
  vr  = Vrad->Field;
  i = 1;
#pragma omp parallel for private(l)
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    rho[l-ns] = rho[l];		/* copy first ring into ghost ring */
    if ((vr[l+ns] > 0.0) || (rho[l] < SigmaMed[0]))
      vr[l] = 0.0; /* we just allow outflow [inwards] */
    else
      vr[l] = vr[l+ns];
  }
}

void NonReflectingBoundary (Vrad, Rho)
PolarGrid *Vrad, *Rho;
{
  int i,j,l,ns,nr,jp,lp,i_angle;
  real *rho, *vr;
  real dangle, mean, vr_med;
  ns = Rho->Nsec;
  nr = Rho->Nrad;
  rho = Rho->Field;
  vr  = Vrad->Field;
  if (CPU_Rank == 0) {
    i = 1;			/* The expression below should be refined */
    /* We need to know the orbital frequency of the nearest planet */
    dangle = (pow(Rinf[1],-1.5)-1.0)/(.5*(SOUNDSPEED[0]+SOUNDSPEED[1]));
    dangle *= (Rmed[1]-Rmed[0]);
    i_angle = (int)(dangle/2.0/PI*(real)NSEC+.5);
#pragma omp parallel for private(l,jp,lp,vr_med)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      jp = j+i_angle;
      if (jp >= ns) jp -= ns;
      if (jp < 0) jp += ns;
      lp = jp;
      rho[lp] = rho[l];		/* copy first ring into ghost ring */
      vr_med = -SOUNDSPEED[i]*(rho[l]-SigmaMed[1])/SigmaMed[1];
      vr[l] = 2.*vr_med-vr[l+ns];
    }
    mean = 0.0;
    for (j = 0; j < ns; j++) {
      mean += rho[j];
    }
    mean /= (real)ns;
    for (j = 0; j < ns; j++) {
      rho[j] += SigmaMed[0]-mean;
    }
  }
  if (CPU_Rank == CPU_Number-1) {
    i = nr-1;			/* The expression below should be refined */
    /* We need to know the orbital frequency of the nearest planet */
    dangle = (pow(Rinf[nr-2],-1.5)-1.0)/(.5*(SOUNDSPEED[nr-1]+SOUNDSPEED[nr-2]));
    dangle *= (Rmed[nr-1]-Rmed[nr-2]);
    i_angle = (int)(dangle/2.0/PI*(real)NSEC+.5);
#pragma omp parallel for private(l,jp,lp,vr_med)
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      jp = j-i_angle;
      if (jp >= ns) jp -= ns;
      if (jp < 0) jp += ns;
      lp = jp+(i-1)*ns;
      rho[l] = rho[lp];		/* copy first ring into ghost ring */
      vr_med = SOUNDSPEED[i]*(rho[l-ns]-SigmaMed[nr-2])/SigmaMed[nr-2];
      vr[l] = 2.*vr_med-vr[l-ns];
    }
    mean = 0.0;
    for (j = 0; j < ns; j++) {
      mean += rho[j+ns*(nr-1)];
    }
    mean /= (real)ns;
    for (j = 0; j < ns; j++) {
      rho[j+(nr-1)*ns] += SigmaMed[nr-1]-mean;
    }
  }
}

void ApplyOuterSourceMass (Rho, Vrad)
     PolarGrid *Rho, *Vrad;
{
  int i,j,l,nr,ns;
  real *rho, average_rho = 0.0, *vr, penul_vr;
  if (CPU_Rank != CPU_Number-1) return;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho= Rho->Field;
  vr = Vrad->Field;
  i = nr-1;
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    average_rho += rho[l];
  }
  average_rho /= (real)ns;
  average_rho = SigmaMed[nr-1]-average_rho;
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    rho[l] += average_rho;
  }
  i = nr-1;
  penul_vr = IMPOSEDDISKDRIFT*pow((Rinf[nr-1]/1.0),-SIGMASLOPE);
  for (j = 0; j < ns; j++) {
    l = j+i*ns;
    vr[l] = penul_vr;
  }
}


void ApplyBoundaryCondition (Vrad, Vtheta, Rho, dt)
PolarGrid *Vrad, *Vtheta, *Rho;
real dt;
{
  if (Stockholm == YES) {
    StockholmBoundary (Vrad, Vtheta, Rho, dt);
    return;
  }
  if (OpenInner == YES) OpenBoundary (Vrad, Rho);
  if (NonReflecting == YES) NonReflectingBoundary (Vrad, Rho);
  if (OuterSourceMass == YES) ApplyOuterSourceMass (Rho, Vrad);
}

void CorrectVtheta (vtheta, domega)
PolarGrid *vtheta;
real domega;
{
  int i,j,l,nr,ns;
  real *vt;
  nr = vtheta->Nrad;
  ns = vtheta->Nsec;
  vt = vtheta->Field;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      vt[l] -= domega*Rmed[i];
    }
  }
}
 
