#ifndef SIDEEULER_H
#define SIDEEULER_H

#include "data.h"
#include "types.h"

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
void ApplySubKeplerianBoundaryInner(t_polargrid &v_azimuthal);
void ApplySubKeplerianBoundaryOuter(t_polargrid &v_azimuthal,
				    const bool did_sg);

void correct_v_azimuthal(t_polargrid &v_azimuthal, double dOmega);

#endif // SIDEEULER_H
