/***************************************************************************
 *   Copyright (C) 2011 by Tobias Mueller                                  *
 *   Tobias_Mueller@twam.info                                              *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANYs WARRANTY; without even the implied warranty of       *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>

// g in cgs
const double G = 6.6742867e-8;
const double AU = 149.60e11;

int main(int argc, char* argv[]) {
	// check parameter count
	if (argc < 5) {
		fprintf(stderr, "Not enough parameters!\n");
		fprintf(stderr, "Usage: %s <n_radial> <surface density 1D file> <soundspeed 1D file> <v_azimthal 1D file>\n",argv[0]);
		return EXIT_FAILURE;
	}

	unsigned int N_radial = atoi(argv[1]);

	// read radius & surface density 
	FILE *file_sigma;
	double *radius = new double[N_radial];
	double *sigma = new double[N_radial];

	if ((file_sigma = fopen(argv[2], "r")) == NULL) {
		fprintf(stderr, "Could not open file '%s'!\n", argv[2]);
		return EXIT_FAILURE;		
	}

	for (unsigned int line = 0; line < N_radial; ++line) {
		// read radius
		if (fread(&radius[line], sizeof(double), 1, file_sigma) == 0) {
			fprintf(stderr, "File '%s' is too short!\n", argv[2]);
			fclose(file_sigma);
			return EXIT_FAILURE;
		}
		radius[line]*=AU;

		// read mean
		if (fread(&sigma[line], sizeof(double), 1, file_sigma) == 0) {
			fprintf(stderr, "File '%s' is too short!\n", argv[2]);
			fclose(file_sigma);
			return EXIT_FAILURE;
		}

		// skip min/max
		fseek(file_sigma, 2*sizeof(double), SEEK_CUR);
	}
	fclose(file_sigma);

	// read soundspeed
	FILE *file_cs;
	double *cs = new double[N_radial];

	if ((file_cs = fopen(argv[3], "r")) == NULL) {
		fprintf(stderr, "Could not open file '%s'!\n", argv[3]);
		return EXIT_FAILURE;		
	}

	for (unsigned int line = 0; line < N_radial; ++line) {
		// skip radius
		fseek(file_cs, 1*sizeof(double), SEEK_CUR);

		// read mean
		if (fread(&cs[line], sizeof(double), 1, file_cs) == 0) {
			fprintf(stderr, "File '%s' is too short!\n", argv[3]);
			fclose(file_cs);
			return EXIT_FAILURE;
		}

		// skip min/max
		fseek(file_cs, 2*sizeof(double), SEEK_CUR);
		//printf("%u %lf %lg\n", line, radius[line], sigma[line]);
	}
	fclose(file_cs);
	
	// read v_azimuthal
	FILE *file_v;
	double *v = new double[N_radial];

	if ((file_v = fopen(argv[4], "r")) == NULL) {
		fprintf(stderr, "Could not open file '%s'!\n", argv[4]);
		return EXIT_FAILURE;		
	}

	for (unsigned int line = 0; line < N_radial; ++line) {
		// skip radius
		fseek(file_v, 1*sizeof(double), SEEK_CUR);

		// read mean
		if (fread(&v[line], sizeof(double), 1, file_v) == 0) {
			fprintf(stderr, "File '%s' is too short!\n", argv[4]);
			fclose(file_v);
			return EXIT_FAILURE;
		}

		// skip min/max
		fseek(file_v, 2*sizeof(double), SEEK_CUR);
	}
	fclose(file_v);

	// calculate lambda
	double *lambda = new double[N_radial];

	for (unsigned int i = 0; i < N_radial; ++i) {
		lambda[i] = sqrt(M_PI*sqrt(2.0*M_PI)*pow(cs[i],3.0)*radius[i]/(G*sigma[i]*v[i]));
		//printf("%u %lg\n", i, lambda[i]);
	}
	
	// calculate J
	double *J = new double[N_radial];

	printf("# radius J (see Truelove et al., 1997)\n");
	for (unsigned int i = 0; i < N_radial-1; ++i) {
		J[i] = (radius[i+1]-radius[i])/lambda[i];
		printf("%lg %lg\n", radius[i]/AU, J[i]);
	}

	delete [] J;
	delete [] lambda;
	delete [] cs;
	delete [] radius;
	delete [] sigma;

	return EXIT_SUCCESS;
}
