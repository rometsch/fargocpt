#include "Theo.h"

/**
	\param Density
*/
void RefillSigma(t_polargrid* Density)
{
	unsigned int nRadial, nAzimuthal, cell;
	double mean;

	for (nRadial = 0; nRadial < Density->Nrad; ++nRadial) {
		mean = 0.0;
		for (nAzimuthal = 0; nAzimuthal < Density->Nsec; ++nAzimuthal) {
			cell = nAzimuthal+nRadial*Density->Nsec;
			mean += Density->Field[cell];
		}
		mean /= (double)(Density->Nsec);
		SigmaMed[nRadial] = mean;
	}

	SigmaInf[0] = SigmaMed[0];

	for (nRadial = 1; nRadial < Density->Nrad; ++nRadial) {
		SigmaInf[nRadial] = (SigmaMed[nRadial-1]*(Rmed[nRadial]-Rinf[nRadial])+SigmaMed[nRadial]*(Rinf[nRadial]-Rmed[nRadial-1]))/(Rmed[nRadial]-Rmed[nRadial-1]);
	}
}

/**

*/
void RefillEnergy(t_polargrid* Energy)
{
	unsigned int nRadial, nAzimuthal, cell;
	double mean;

	for (nRadial = 0; nRadial < Energy->Nrad; ++nRadial) {
		mean = 0.0;
    		for (nAzimuthal = 0; nAzimuthal < Energy->Nsec; ++nAzimuthal) {
			cell = nAzimuthal+nRadial*Energy->Nsec;
			mean += Energy->Field[cell];
		}
		mean /= (double)Energy->Nsec;
		EnergyMed[nRadial] = mean;
	}
}

double calculate_omega_kepler(double r)
{
	return sqrt(constants::G * hydro_center_mass / (r * r * r));
}
