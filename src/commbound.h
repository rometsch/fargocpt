#ifndef _COMMBOUND_H_
#define _COMMBOUND_H_

#include "types.h"

void AllocateBoundaryCommunicationBuffers();
void CommunicateBoundaries(t_polargrid* Density, t_polargrid* Vrad, t_polargrid* Vtheta, t_polargrid* Energy);

#endif
