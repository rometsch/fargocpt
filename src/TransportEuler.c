/** \file TransportEuler.c

Functions that handle the transport substep of a hydrodynamical time
step.  The FARGO algorithm is implemented here. The transport is
performed in a manner similar to what is done for the ZEUS code (Stone
& Norman, 1992), except for the momenta transport (we define a left
and right momentum for each zone, which we declare zone centered; we
then transport then normally, and deduce the new velocity in each zone
by a proper averaging).

*/

#include "fargo.h"

extern real OmegaFrame;

static real VMed[MAX1D];
static int Nshift[MAX1D];
static real *TempShift;

static real *dq;
static PolarGrid *RadMomP, *RadMomM, *ThetaMomP, *ThetaMomM, *ExtLabel, *VthetaRes;
static PolarGrid *Work, *QRStar, *Elongations;

extern int TimeStep;
extern boolean OpenInner, FastTransport;

real LostMass = 0.0;

void Transport (Rho, Vrad, Vtheta, Label, dt)
PolarGrid *Rho, *Vrad, *Vtheta, *Label;
real dt;
{
  ComputeLRMomenta (Rho, Vrad, Vtheta);
  if (AdvecteLabel == YES)
     ComputeExtQty (Rho, Label, ExtLabel);
				/* No-Alternate Directionnal Splitting */
  OneWindRad (Rho, Vrad, dt);
  OneWindTheta (Rho, Vtheta, dt);
  ComputeVelocities (Rho, Vrad, Vtheta);
  if (AdvecteLabel == YES) 
     ComputeSpeQty (Rho, Label, ExtLabel);
}

void OneWindRad (Rho, Vrad, dt)
PolarGrid *Rho, *Vrad;
real dt;
{
  ComputeStarRad (Rho, Vrad, RhoStar, dt);
  ActualiseGas (RhoInt, Rho);
  VanLeerRadial (Vrad, RadMomP, dt); 
  VanLeerRadial (Vrad, RadMomM, dt);  
  VanLeerRadial (Vrad, ThetaMomP, dt);   
  VanLeerRadial (Vrad, ThetaMomM, dt);
  if (AdvecteLabel == YES)
    VanLeerRadial (Vrad, ExtLabel, dt);
  LostMass += VanLeerRadial (Vrad, Rho, dt); /* MUST be the last line */
}

	
			/* Hereafter are the new specific procedures to the fast algorithm */

void ComputeResiduals (Vtheta, dt)
PolarGrid *Vtheta;
real dt;
{
  int i,j,l,nr,ns;
  long nitemp;
  real *vt, moy, *vtr, Ntilde, Nround;
  real invdt, dpinvns, vtemp;
  nr = Vtheta->Nrad;
  ns = Vtheta->Nsec;
  vt = Vtheta->Field;
  vtr= VthetaRes->Field;
  invdt = 1.0/dt;
  dpinvns=2.0*PI/(real)ns;
#pragma omp parallel for private(j,l, moy, Ntilde, Nround, nitemp, vtemp)
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l=j+i*ns;
      moy += vt[l];
    }
    VMed[i] = moy/(real)ns;
    Ntilde = VMed[i]*InvRmed[i]*dt*(real)ns/2.0/PI;
    Nround = floor(Ntilde+0.5);
    nitemp = (long)Nround;
    Nshift[i] = (long)nitemp;
    if (FastTransport == YES) {
      vtemp = (Ntilde-Nround)*Rmed[i]*invdt*dpinvns;
      for (j = 0; j < ns; j++) {
	l=j+i*ns;
	vtr[l] = vt[l]-VMed[i];
	vt[l] = vtemp; 
      }
    } else {
      for (j = 0; j < ns; j++) {
	l=j+i*ns;
	vtr[l] = vt[l];
	vt[l] = 0.0; 
      }
    }
  }
}

void AdvectSHIFT (array)
PolarGrid *array;
{
  int i,j,ji,l,li,nr,ns;
  real *val;
  val = array->Field;
  nr  = array->Nrad;
  ns  = array->Nsec;
#pragma omp parallel for private (j,ji,l,li) default(none) shared(nr,ns,Nshift, val, TempShift)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      ji = j-Nshift[i];
      while (ji < 0) ji += ns;
      while (ji >= ns) ji -= ns;
      li= ji+i*ns;
      l=j + i*ns;
      TempShift[l]=val[li];
    }
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      val[l] = TempShift[l];
    }
  }
}

void OneWindTheta (Rho, Vtheta, dt)
PolarGrid *Rho, *Vtheta;
real dt;
{
  ComputeResiduals (Vtheta, dt);   /* Constant residual is in Vtheta from now on */
  QuantitiesAdvection (Rho, VthetaRes, dt);
  if (FastTransport == YES) { /* Useless in standard transport as Vtheta is zero */
    QuantitiesAdvection (Rho, Vtheta, dt); /* Uniform Transport here */
    AdvectSHIFT (RadMomP);
    AdvectSHIFT (RadMomM);
    AdvectSHIFT (ThetaMomP);
    AdvectSHIFT (ThetaMomM);
    if (AdvecteLabel == YES)
      AdvectSHIFT (ExtLabel);
    AdvectSHIFT (Rho);
  }
}
				/* End of new specific procedures to the fast algorithm */

void QuantitiesAdvection (Rho, Vtheta, dt)
PolarGrid *Rho, *Vtheta;
real dt;
{
  ComputeStarTheta (Rho, Vtheta, RhoStar, dt);
  ActualiseGas (RhoInt, Rho);
  VanLeerTheta (Vtheta, RadMomP, dt);
  VanLeerTheta (Vtheta, RadMomM, dt);    
  VanLeerTheta (Vtheta, ThetaMomP, dt);    
  VanLeerTheta (Vtheta, ThetaMomM, dt); 
  if (AdvecteLabel == YES)
    VanLeerTheta (Vtheta, ExtLabel, dt); 
  VanLeerTheta (Vtheta, Rho, dt); /* MUST be the last line */
}

void ComputeExtQty (Rho, Label, ExtLabel)
PolarGrid *Rho, *Label, *ExtLabel;
{
  int i,j,l,nr,ns;
  real *rho, *lab, *extlab;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho= Rho->Field;
  lab= Label->Field;
  extlab= ExtLabel->Field;
#pragma omp parallel for private(j,l)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      extlab[l] = rho[l]*lab[l]; /* compressive flow  if line commentarized
				    extlab[l] = lab[l]; */
    }
  }
}

void ComputeSpeQty (Rho, Label, ExtLabel)
PolarGrid *Rho, *Label, *ExtLabel;
{
  int i,j,l,nr,ns;
  real *rho, *lab, *extlab;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho= Rho->Field;
  lab= Label->Field;
  extlab= ExtLabel->Field;
#pragma omp parallel for private(j,l)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lab[l] = extlab[l]/rho[l]; /* compressive flow if line commentarized
				    lab[l] = extlab[l]; */
    }
  }
}

void InitTransport () 
{
  RadMomP      = CreatePolarGrid(NRAD, NSEC, "RadMomP");
  RadMomM      = CreatePolarGrid(NRAD, NSEC, "RadMomM");
  ThetaMomP    = CreatePolarGrid(NRAD, NSEC, "ThetaMomP");
  ThetaMomM    = CreatePolarGrid(NRAD, NSEC, "ThetaMomM");
  Work         = CreatePolarGrid(NRAD, NSEC, "WorkGrid");
  QRStar       = CreatePolarGrid(NRAD, NSEC, "QRStar");
  ExtLabel     = CreatePolarGrid(NRAD, NSEC, "ExtLabel");
  VthetaRes    = CreatePolarGrid(NRAD, NSEC, "VThetaRes");
  Elongations  = CreatePolarGrid(NRAD, NSEC, "Elongations");
  TempShift    = (real *)malloc(NRAD*NSEC*sizeof(real));
  dq           = (real *)malloc(NRAD*NSEC*sizeof(real));
}

void ComputeStarRad (Qbase, Vrad, QStar, dt)
PolarGrid *Qbase, *Vrad, *QStar;
real dt;
{
  int i,j,l,lq,lip,lim,nr,ns;
  real *qb, *qs, *vr;
  real dqp, dqm;
  nr = Qbase->Nrad;
  ns = Qbase->Nsec;
  qb = Qbase->Field;
  qs = QStar->Field;
  vr = Vrad->Field;
#pragma omp parallel for private (i,l,lq,lip,lim,dqm,dqp)
  for (j = 0; j < ns; j++) {
    for (i = 0; i < nr; i++) {
      l = j+i*ns;
      lq= i+j*nr;
      lip = l+ns;
      lim = l-ns;
      if ((i == 0) || (i == nr-1)) dq[lq] = 0.0;
      else {
	dqm = (qb[l]-qb[lim])*InvDiffRmed[i];
	dqp = (qb[lip]-qb[l])*InvDiffRmed[i+1];
	if (dqp * dqm > 0.0)
	  dq[lq] = 2.0*dqp*dqm/(dqp+dqm);
	else
	  dq[lq] = 0.0;
      }
    }
    for (i = 0; i < nr; i++) {
      l = j+i*ns;
      lq= i+j*nr;
      lip = l+ns;
      lim = l-ns;
      if (vr[l] > 0.0)
	qs[l] = qb[lim]+(Rmed[i]-Rmed[i-1]-vr[l]*dt)*0.5*dq[lq-1];
      else
	qs[l] = qb[l]-(Rmed[i+1]-Rmed[i]+vr[l]*dt)*0.5*dq[lq];
    }
    qs[j] = qs[j+ns*nr] = 0.0;
  }
}
	
void ComputeStarTheta (Qbase, Vtheta, QStar, dt)
PolarGrid *Qbase, *Vtheta, *QStar;
real dt;
{
  int i,j,l,ljp,ljm,jm,nr,ns;
  real *qb, *qs, *vt;
  real dqp, dqm,dxtheta,ksi,invdxtheta;
  nr = Qbase->Nrad;
  ns = Qbase->Nsec;
  qb = Qbase->Field;
  qs = QStar->Field;
  vt = Vtheta->Field;
#pragma omp parallel for private (dxtheta,invdxtheta,l,j,ljp,ljm,dqm,dqp,jm,ksi)
  for (i = 0; i < nr; i++) {
    dxtheta = 2.0*PI/(real)ns*Rmed[i];
    invdxtheta = 1.0/dxtheta;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      ljp = l+1;
      ljm = l-1;
      if (j == 0) ljm = i*ns+ns-1;
      if (j == ns-1) ljp = i*ns;
      dqm = (qb[l]-qb[ljm]);
      dqp = (qb[ljp]-qb[l]);
      if (dqp * dqm > 0.0)
	dq[l] = dqp*dqm/(dqp+dqm)*invdxtheta;
      else
	dq[l] = 0.0;
    }
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      jm = j-1; 
      if (j == 0) jm = ns-1;
      ljm = jm+i*ns;
      ksi=vt[l]*dt;
      if (ksi > 0.0)
	qs[l] = qb[ljm]+(dxtheta-ksi)*dq[ljm];
      else
	qs[l] = qb[l]-(dxtheta+ksi)*dq[l];
    }
  }
}
	
void ComputeLRMomenta (Rho, Vrad, Vtheta)
PolarGrid *Rho, *Vrad, *Vtheta;
{
  int i,j,l,lip,ljp,nr,ns;
  real *vr,*vt,*rho;
  real *rp, *rm, *tp, *tm;
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  rho= Rho->Field;
  vr = Vrad->Field;
  vt = Vtheta->Field;
  rp = RadMomP->Field;
  rm = RadMomM->Field;
  tp = ThetaMomP->Field;
  tm = ThetaMomM->Field;
#pragma omp parallel for private(j,l,lip,ljp)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1)
	ljp = i*ns;
      rp[l] = rho[l]*vr[lip];
      rm[l] = rho[l]*vr[l];
      tp[l] = rho[l]*(vt[ljp]+Rmed[i]*OmegaFrame)*Rmed[i]; /* it is the angular momentum */
      tm[l] = rho[l]*(vt[l]+Rmed[i]*OmegaFrame)*Rmed[i];
    }
  }
}


void ComputeVelocities (Rho, Vrad, Vtheta)
PolarGrid *Rho, *Vrad, *Vtheta;
{
  int i,j,l,lim,ljm,nr,ns;
  real *vr,*vt,*rho;
  real *rp, *rm, *tp, *tm;
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  rho= Rho->Field;
  vr = Vrad->Field;
  vt = Vtheta->Field;
  rp = RadMomP->Field;
  rm = RadMomM->Field;
  tp = ThetaMomP->Field;
  tm = ThetaMomM->Field;
#pragma omp parallel for private(j,l,lim,ljm)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lim = l-ns;
      ljm = l-1;
      if (j == 0)
	ljm = i*ns+ns-1;
      if (i == 0)
	vr[l] = 0.0;
      else		
	vr[l] = (rp[lim]+rm[l])/(rho[l]+rho[lim]+1e-20);
      vt[l] = (tp[ljm]+tm[l])/(rho[l]+rho[ljm]+1e-15)/Rmed[i]-Rmed[i]*OmegaFrame;
				/* It was the angular momentum */
    }
  }
}


real VanLeerRadial (Vrad, Qbase, dt)
PolarGrid *Vrad, *Qbase;
real dt;
{
  int i,j,nr,ns,l,lip;
  real *qrs, *rhos, *vr, *qb;
  real dtheta, varq;
  real LostByDisk=0.0;
  DivisePolarGrid (Qbase, RhoInt, Work);
  ComputeStarRad (Work, Vrad, QRStar, dt);
  nr = Qbase->Nrad;
  ns = Qbase->Nsec;
  qrs  = QRStar->Field;
  rhos = RhoStar->Field;
  vr = Vrad->Field;
  qb=  Qbase->Field;
  dtheta = 2.0*PI/(real)ns;
#pragma omp parallel for private(j,l,lip,varq)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l=j+i*ns;
      lip=l+ns;
      varq =dt*dtheta*Rinf[i]*qrs[l]*rhos[l]*vr[l];
      varq-=dt*dtheta*Rsup[i]*qrs[lip]*rhos[lip]*vr[lip];
      qb[l] += varq*InvSurf[i];
      if ((i == 0) && (OpenInner == YES))
#pragma omp atomic
	LostByDisk += varq;
    }
  }
  return LostByDisk;
}
      
void VanLeerTheta (Vtheta, Qbase, dt)
PolarGrid *Vtheta, *Qbase;
real dt;
{
  int i,j,nr,ns,l,ljp;
  real *qrs, *rhos, *vt, *qb;
  real dxrad, varq, invsurf;
  DivisePolarGrid (Qbase, RhoInt, Work);
  ComputeStarTheta (Work, Vtheta, QRStar, dt);
  nr = Qbase->Nrad;
  ns = Qbase->Nsec;
  qrs  = QRStar->Field;
  rhos = RhoStar->Field;
  vt = Vtheta->Field;
  qb=  Qbase->Field;
#pragma omp parallel for private(dxrad,invsurf,j,l,ljp,varq)
  for (i = 0; i < nr; i++) {
    dxrad = (Rsup[i]-Rinf[i])*dt;
    invsurf = 1.0/Surf[i];
    for (j = 0; j < ns; j++) {
      l=j+i*ns;
      ljp=l+1;
      if (j == ns-1) ljp=i*ns;
      varq  = dxrad*qrs[l]*rhos[l]*vt[l];
      varq -= dxrad*qrs[ljp]*rhos[ljp]*vt[ljp];
      qb[l] += varq*invsurf;
    }
  }
}


 
