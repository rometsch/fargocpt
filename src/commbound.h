#pragma once

#include "data.h"

void DeallocateBoundaryCommunicationBuffers();
void AllocateBoundaryCommunicationBuffers();
void CommunicateBoundaries(t_polargrid *Density, t_polargrid *Vrad,
			   t_polargrid *Vazi, t_polargrid *Energy);

void CommunicateBoundariesAll(t_data& data);
void CommunicateBoundariesAllInitial(t_data& data);