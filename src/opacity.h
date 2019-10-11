#ifndef OPAC_H
#define OPAC_H

namespace opacity
{

double opacity(double density, double temperature);
double lin(double density, double temperature);
double bell(double density, double temperature);
double zhu(double density, double temperature);
double kramers(double density, double temperature);

}; // namespace opacity

#endif // OPAC_H
