/** \file split.c

Split (radially) the mesh among the different processors.  A simple
Round Robin algorithm is used to achieve a proper load balancing,
assuming the cluster to be homogeneous.  A ring described by a given
process has its 5 (CPUOVERLAP) innermost zones that are ghost zones
which are filled by communications with the previous inner process
(unless it is the innermost process itself), and its 5 (CPUOVERLAP)
outermost zones that are ghost zones which are filled by
communications with the next outer process (unless it is the outermost
process itself). The "active" part of the submesh described by a given
process is the part of this mesh that excludes the ghost zones.  For
each process, the (local) radial index of the first active ring is
Zero_or_active, and the (local) radial index of the last active ring
is Max_or_active. MaxMO_or_active is Max_or_active for all processes,
except the very last one (the outermost) for which it is
Max_or_active-1 (MO stands for 'Minus One').

*/


#include "fargo.h"

void SplitDomain () {
  int remainder;
  int size_low, size_high;
  GLOBALNRAD = NRAD;
  size_low = NRAD/CPU_Number;
  size_high = size_low+1;
  remainder = NRAD % CPU_Number;
  if (size_low < 2*CPUOVERLAP) {
    mastererr("The number of processes is too large\n");
    mastererr("or the mesh is radially too narrow.\n");
    prs_exit(1);
  }
  if (CPU_Rank < remainder) {
    IMIN = size_high*CPU_Rank;
    IMAX = IMIN+size_high-1;
  } else {
    IMIN = size_high*remainder+(CPU_Rank-remainder)*size_low;
    IMAX = IMIN+size_low-1;
  }
  if (CPU_Rank > 0) IMIN -= CPUOVERLAP;
  if (CPU_Rank < CPU_Number-1) IMAX +=CPUOVERLAP;
  NRAD = IMAX-IMIN+1;
  Zero_or_active = CPUOVERLAP * (CPU_Rank > 0 ? 1 : 0);
  One_or_active = 1+(CPUOVERLAP-1) * (CPU_Rank > 0 ? 1 : 0);
  Max_or_active = NRAD-CPUOVERLAP*(CPU_Rank < CPU_Number-1 ? 1 : 0);
  MaxMO_or_active = NRAD-1-(CPUOVERLAP-1)*(CPU_Rank < CPU_Number-1 ? 1 : 0);
  if (debug == YES) {
    masterprint ("%d = %d * %d + %d\n", GLOBALNRAD, CPU_Number, size_low, remainder);
    printf ("IMIN process %d : %d\n", CPU_Rank, IMIN);
    printf ("IMAX process %d : %d\n", CPU_Rank, IMAX);
    printf ("NRAD process %d : %d\n", CPU_Rank, NRAD);
    printf ("Zero_or_active process %d : %d\n", CPU_Rank, Zero_or_active);
    printf ("One_or_active process %d : %d\n", CPU_Rank, One_or_active);
    printf ("Max_or_active process %d : %d\n", CPU_Rank, Max_or_active);
    printf ("MaxMO_or_active process %d : %d\n", CPU_Rank, MaxMO_or_active);
    printf ("GLOB process %d : %d\n", CPU_Rank, GLOBALNRAD);
  }
}
      
