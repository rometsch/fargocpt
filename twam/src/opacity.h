#ifndef _OPAC_H_
#define _OPAC_H_

namespace opacity {

double opacity(double density, double temperature);
double lin(double density, double temperature);
double bell(double density, double temperature);
double zhu(double density, double temperature);
double kramers(double density, double temperature);

};

#endif
