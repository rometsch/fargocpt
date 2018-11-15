#ifndef COMMBOUND_H
#define COMMBOUND_H

#include "types.h"

void AllocateBoundaryCommunicationBuffers();
void CommunicateBoundaries(t_polargrid* Density, t_polargrid* Vrad, t_polargrid* Vtheta, t_polargrid* Energy);

#endif // COMMBOUND_H
