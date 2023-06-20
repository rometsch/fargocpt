#pragma once

#include "types.h"
#include "data.h"

void DeallocateBoundaryCommunicationBuffers();
void AllocateBoundaryCommunicationBuffers();
void CommunicateBoundaries(t_polargrid *Density, t_polargrid *Vrad,
			   t_polargrid *Vtheta, t_polargrid *Energy);

void CommunicateBoundariesAll(t_data& data);
void CommunicateBoundariesAllInitial(t_data& data);