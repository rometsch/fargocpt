#ifndef TRANSPORTEULER_H
#define TRANSPORTEULER_H

#include "data.h"
#include "types.h"

void Transport(t_data &data, PolarGrid* Density, PolarGrid* VRadial, PolarGrid* VAzimuthal, PolarGrid* Energy, double dt);
void OneWindRad(t_data &data, PolarGrid* Density, PolarGrid* VRadial, PolarGrid* Energy, double dt);
void ComputeThetaElongations(PolarGrid* VAzimuthal, double dt);

void compute_average_azimuthal_velocity(t_polargrid &v_azimuthal);
void compute_residual_velocity(t_polargrid &v_azimuthal);

void ComputeConstantResidual(PolarGrid* VAzimuthal, double dt);
void AdvectSHIFT(t_polargrid &array);
void OneWindTheta(t_data &data, PolarGrid* Density, PolarGrid* VAzimuthal, PolarGrid* Energy, double dt);
void QuantitiesAdvection(t_data &data, PolarGrid* Density, PolarGrid* VAzimuthal, PolarGrid* Energy, double dt);
void InitTransport();
void FreeTransport();
void compute_star_radial(t_polargrid* Qbase, t_polargrid* VRadial, t_polargrid* QStar, double dt);
void ComputeStarTheta(PolarGrid* Qbase, PolarGrid* VAzimuthal, PolarGrid* QStar, double dt);

void compute_momenta_from_velocities(t_polargrid &density, t_polargrid &v_radial, t_polargrid &v_azimuthal);
void compute_velocities_from_momenta(t_polargrid &density, t_polargrid &v_radial, t_polargrid &v_azimuthal);

void VanLeerRadial(t_data &data, PolarGrid* VRadial, PolarGrid* Qbase, double dt);
void VanLeerTheta(t_data &data, PolarGrid* VAzimuthal, PolarGrid* Qbase, double dt);
void boundary_layer_mass_influx(PolarGrid* QStar, PolarGrid* VRadial);


#endif // TRANSPORTEULER_H
