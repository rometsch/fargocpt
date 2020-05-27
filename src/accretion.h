#pragma once

#include "data.h"

namespace accretion
{
	bool AccreteOntoSinglePlanet(t_data &data, double dt);
	void AccreteOntoPlanets(t_data &data, double dt);
} // namespace accretion
