#pragma once

#include "data.h"

void CheckAngularMomentumConservation(t_data &data);
void divise_polargrid(const t_polargrid& num, const t_polargrid& denom,
			  t_polargrid &result);
void InitCellCenterCoordinates();
void FreeCellCenterCoordinates();

void correct_v_azimuthal(t_polargrid &v_azimuthal, double dOmega);
