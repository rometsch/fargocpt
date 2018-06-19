/** \file Psys.c

Contains the functions that set up the planetary system configuration.
In addition, the last two functions allow to track the first planet
(number 0) of the planetary system, in order to perform a calculation
in the frame corotating either with this planet or with its
guiding-center.

*/

#include "fargo.h"

static real Xplanet, Yplanet;

extern boolean GuidingCenter;

int FindNumberOfPlanets (filename)
char *filename;
{
  FILE *input;
  char s[512];
  int Counter=0;
  input = fopen (filename, "r");
  if (input == NULL) {
    fprintf (stderr, "Error : can't find '%s'.\n", filename);
    prs_exit (1);
  }
  while (fgets(s, 510, input) != NULL) {
    if (isalpha(s[0]))
      Counter++;
  }
  fclose (input);
  return Counter;
}

PlanetarySystem *AllocPlanetSystem (nb)
     int nb;
{
  real *mass, *x, *y, *vx, *vy, *acc;
  boolean *feeldisk, *feelothers;
  int i;
  PlanetarySystem *sys;
  sys  = (PlanetarySystem *)malloc (sizeof(PlanetarySystem));
  if (sys == NULL) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  x    = (real *)malloc (sizeof(real)*(nb+1));
  y    = (real *)malloc (sizeof(real)*(nb+1));
  vy   = (real *)malloc (sizeof(real)*(nb+1));
  vx   = (real *)malloc (sizeof(real)*(nb+1));
  mass = (real *)malloc (sizeof(real)*(nb+1));
  acc  = (real *)malloc (sizeof(real)*(nb+1));
  if ((x == NULL) || (y == NULL) || (vx == NULL) || (vy == NULL) || (acc == NULL) || (mass == NULL)) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  feeldisk   = (boolean *)malloc (sizeof(real)*(nb+1));
  feelothers = (boolean *)malloc (sizeof(real)*(nb+1));
  if ((feeldisk == NULL) || (feelothers == NULL)) {
    fprintf (stderr, "Not enough memory.\n");
    prs_exit (1);
  }
  sys->x = x;
  sys->y = y;
  sys->vx= vx;
  sys->vy= vy;
  sys->acc=acc;
  sys->mass = mass;
  sys->FeelDisk = feeldisk;
  sys->FeelOthers = feelothers;
  for (i = 0; i < nb; i++) {
    x[i] = y[i] = vx[i] = vy[i] = mass[i] = acc[i] = 0.0;
    feeldisk[i] = feelothers[i] = YES;
  }
  return sys;
}

void FreePlanetary (sys)
PlanetarySystem *sys;
{
  free (sys->x);
  free (sys->vx);
  free (sys->y);
  free (sys->vy);
  free (sys->mass);
  free (sys->acc);
  free (sys->FeelOthers);
  free (sys->FeelDisk);
  free (sys);
}

PlanetarySystem *InitPlanetarySystem (filename)
char *filename;
{
  FILE *input;
  char s[512], nm[512], test1[512], test2[512], *s1;
  PlanetarySystem *sys;
  int i=0, nb;
  float mass, dist, accret;
  boolean feeldis, feelothers;
  nb = FindNumberOfPlanets (filename);
  if (CPU_Master)
    printf ("%d planet(s) found.\n", nb);
  sys = AllocPlanetSystem (nb);
  input = fopen (filename, "r");
  sys->nb = nb;
  while (fgets(s, 510, input) != NULL) {
    sscanf(s, "%s ", nm);
    if (isalpha(s[0])) {
      s1 = s + strlen(nm);
      sscanf(s1 + strspn(s1, "\t :=>_"), "%f %f %f %s %s", &dist, &mass, &accret, test1, test2);
      sys->mass[i] = (real)mass;
      feeldis = feelothers = YES;
      if (tolower(*test1) == 'n') feeldis = NO;
      if (tolower(*test2) == 'n') feelothers = NO;
      sys->x[i] = (real)dist*(1.0+ECCENTRICITY);
      sys->y[i] = 0.0;
      sys->vy[i] = (real)sqrt((1.0+mass)/dist)*\
	sqrt(1.0-ECCENTRICITY*ECCENTRICITY)/(1.0+ECCENTRICITY);
      sys->vx[i] = -0.0000000001*sys->vy[i];
      sys->acc[i] = accret;
      sys->FeelDisk[i] = feeldis;
      sys->FeelOthers[i] = feelothers;
      i++;
    }
  }
  return sys;
}

void ListPlanets (sys)
PlanetarySystem *sys;
{
  int nb;
  int i;
  nb = sys->nb;
  if (!CPU_Master) return;
  for (i = 0; i < nb; i++) {
    printf ("Planet number %d\n", i);
    printf ("---------------\n");
    printf ("x = %.10f\ty = %.10f\n", sys->x[i],sys->y[i]);
    printf ("vx = %.10f\tvy = %.10f\n", sys->vx[i],sys->vy[i]);
    if (sys->acc[i] == 0.0)
      printf ("Non-accreting.\n");
    else
      printf ("accretion time = %.10f\n", 1.0/(sys->acc[i]));
    if (sys->FeelDisk[i] == YES) {
      printf ("Feels the disk potential\n");
    } else {
      printf ("Doesn't feel the disk potential\n");
    }
    if (sys->FeelOthers[i] == YES) {
      printf ("Feels the other planets potential\n");
    } else {
      printf ("Doesn't feel the other planets potential\n");
    }
    printf ("\n");
  }
}

real GetPsysInfo (sys, action)
PlanetarySystem *sys;
boolean action;
{
  real d1,d2,cross;
  real x,y, vx, vy, m, h, d, Ax, Ay, e, a, E, M;
  real xc, yc, vxc, vyc, omega;
  real arg, PerihelionPA;
  xc = x = sys->x[0];
  yc = y = sys->y[0];
  vxc = vx= sys->vx[0];
  vyc = vy= sys->vy[0];
  m = sys->mass[0]+1.;
  h = x*vy-y*vx;
  d = sqrt(x*x+y*y);
  Ax = x*vy*vy-y*vx*vy-G*m*x/d;
  Ay = y*vx*vx-x*vx*vy-G*m*y/d;
  e = sqrt(Ax*Ax+Ay*Ay)/m;
  a = h*h/G/m/(1.-e*e);
  if (e == 0.0) {
    arg = 1.0;
  } else {
    arg = (1.0-d/a)/e;
  }
  if (fabs(arg) >= 1.0) 
    E = PI*(1.-arg/fabs(arg))/2.;
  else
    E = acos((1.0-d/a)/e);
  if ((x*y*(vy*vy-vx*vx)+vx*vy*(x*x-y*y)) < 0) E= -E;
  M = E-e*sin(E);
  omega = sqrt(m/a/a/a);
  PerihelionPA=atan2(Ay,Ax);
  if (GuidingCenter == YES) {
    xc = a*cos(M+PerihelionPA);
    yc = a*sin(M+PerihelionPA);
    vxc = -a*omega*sin(M+PerihelionPA);
    vyc =  a*omega*cos(M+PerihelionPA);
  } 
  if (e < 1e-8) {
    xc = x;
    yc = y;
    vxc = vx;
    vyc = vy;
  }
  switch (action) {
  case MARK: 
    Xplanet = xc;
    Yplanet = yc;
    return 0.;
    break;
  case GET:
    x = xc;
    y = yc;
    vx = vxc;
    vy = vyc;
    d2 = sqrt(x*x+y*y);
    d1 = sqrt(Xplanet*Xplanet+Yplanet*Yplanet);
    cross = Xplanet*y-x*Yplanet;
    Xplanet = x;
    Yplanet = y;
    return asin(cross/(d1*d2));
    break;
  case FREQUENCY:
    return omega;
    break;
  }
  return 0.0;
}

void RotatePsys (sys, angle)	/* Rotate by angle '-angle' */
PlanetarySystem *sys;
real angle;
{
  int nb;
  int i;
  real sint, cost, xt, yt;
  nb = sys->nb;
  sint = sin(angle);
  cost = cos(angle);
  for (i = 0; i < nb; i++) {
    xt = sys->x[i];
    yt = sys->y[i];
    sys->x[i] = xt*cost+yt*sint;
    sys->y[i] = -xt*sint+yt*cost;
    xt = sys->vx[i];
    yt = sys->vy[i];
    sys->vx[i] = xt*cost+yt*sint;
    sys->vy[i] = -xt*sint+yt*cost;
  }
}
 
