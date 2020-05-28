#ifndef SIDEEULER_H
#define SEIEEULER_H

#include "data.h"
#include "types.h"

void CheckAngularMomentumConservation(t_data &data);
void divise_polargrid(t_polargrid &num, t_polargrid &denom,
		      t_polargrid &result);
void InitCellCenterCoordinates();
void FreeCellCenterCoordinates();

void NonReflectingBoundary_inner(t_data &data, t_polargrid *VRadial,
				 t_polargrid *Density, t_polargrid *Energy);
void NonReflectingBoundary_outer(t_data &data, t_polargrid *VRadial,
				 t_polargrid *Density, t_polargrid *Energy);

void EvanescentBoundary(t_data &data, double step);

void ApplyOuterSourceMass(t_polargrid *Density, PolarGrid *VRadial);
void ApplySubKeplerianBoundary(t_polargrid &v_azimuthal);

void correct_v_azimuthal(t_polargrid &v_azimuthal, double dOmega);

#endif // SEIEEULER_H
