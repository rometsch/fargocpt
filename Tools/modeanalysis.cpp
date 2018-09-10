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

int main(int argc, char* argv[]) {
	FILE *fd;

	// check parameter count
	if (argc < 9) {
		fprintf(stderr, "Not enough parameters!\n");
		fprintf(stderr, "Usage: %s filename n_radial n_azimuthal loggrid rmin rmax max_mode modus\n",argv[0]);
		fprintf(stderr, "  loggrid should be 'l' for logarithmic grid\n");
		fprintf(stderr, "  mode should be 'm' for mass weighted average, 'r' for radial dependency of modes\n");
		return EXIT_FAILURE;
	}

	unsigned int N_radial = atoi(argv[2]);
	unsigned int N_azimuthal = atoi(argv[3]);
	double r_min = atof(argv[5]);
	double r_max = atof(argv[6]);
	unsigned int grid_size = N_radial*N_azimuthal;

	fprintf(stdout, "# Using %ux%u grid.\n",N_radial, N_azimuthal);
	fprintf(stdout, "# r_min = %f\n",r_min);
	fprintf(stdout, "# r_max = %f\n",r_max);

	double* radii = new double[(N_radial+1)*sizeof(double)];

	if (argv[4][0] == 'l')  {
		fprintf(stdout, "# Using logarithmic grid.\n");
		for (unsigned int n_radial = 0; n_radial <= N_radial; ++n_radial) {
			radii[n_radial] = r_min*exp((double)(n_radial-1.0)/(double)(N_radial-2.0)*log(r_max/r_min));
		}
	} else {
		fprintf(stdout, "# Using equidistant grid.\n");
		double interval = (r_max-r_min)/(double)(N_radial-2.0);
		for (unsigned int n_radial = 0; n_radial <= N_radial; ++n_radial) {
			radii[n_radial] = r_min+interval*(double)(n_radial-1.0);
		}
	}

	unsigned int max_mode = atoi(argv[7]);
	printf("# max_mode = %i\n", max_mode);

	struct stat filestatus;
	if (stat(argv[1], &filestatus) == -1) {
		fprintf(stderr, "Could net read file size of '%s'!\n", argv[1]);
		return EXIT_FAILURE;
	}

	if ((unsigned int)filestatus.st_size != N_radial * N_azimuthal * sizeof(double)) {
		fprintf(stderr, "File '%s' has a size of %u bytes, but should be %lu bytes long!\n", argv[1], (unsigned int)filestatus.st_size, N_radial * N_azimuthal * sizeof(double));	
		return EXIT_FAILURE;
	}

	double* grid = new double[grid_size];

	if ((fd = fopen(argv[1],"r")) == NULL) {
		fprintf(stderr, "Could not open file '%s'!\n", argv[1]);
	}

	if (fread(grid, sizeof(double), grid_size, fd) == 0) {
		fprintf(stderr, "Error while reading '%s'!\n", argv[1]);
	}
	fclose(fd);

	double* amplitude = new double[(max_mode+1)*N_radial];
	double* phase = new double[(max_mode+1)*N_radial];

	// Don't treat ghost cells
	for (unsigned int n_radial = 1; n_radial < N_radial-1; ++n_radial) {
		for (unsigned int mode = 0; mode <= max_mode; ++mode) {
			double a = 0.0;
			double b = 0.0;

			// Don't treat ghost cells
			for (unsigned int n_azimuthal = 0; n_azimuthal < N_azimuthal; ++n_azimuthal) {
				a += 2.0*M_PI/(double)N_azimuthal * grid[n_radial*N_azimuthal+n_azimuthal] * cos(2.0*M_PI/(double)N_azimuthal*(double)n_azimuthal*(double)mode);
				b += 2.0*M_PI/(double)N_azimuthal * grid[n_radial*N_azimuthal+n_azimuthal] * sin(2.0*M_PI/(double)N_azimuthal*(double)n_azimuthal*(double)mode);
			}

			amplitude[n_radial*(max_mode+1)+mode] = sqrt(a*a+b*b);
			phase[n_radial*(max_mode+1)+mode] = atan2(b,a);
			if (phase[n_radial*(max_mode+1)+mode] < 0)
				phase[n_radial*(max_mode+1)+mode] += 2.0*M_PI;
		}
	}

	if (argv[8][0] == 'm')  {
		double* surface = new double[N_radial];

		for (unsigned int n_radial = 0; n_radial < N_radial; ++n_radial) {
			double rinf = radii[n_radial];
			double rsup = radii[n_radial+1];
			surface[n_radial] = M_PI*(rsup*rsup-rinf*rinf)/(double)N_azimuthal;
		}

		// mass weighted average
		double total_mass = 0;
		for (unsigned int n_radial = 1; n_radial < N_radial-1; ++n_radial) {
			for (unsigned int n_azimuthal = 0; n_azimuthal < N_azimuthal; ++n_azimuthal) {
				total_mass += surface[n_radial]*grid[n_radial*N_azimuthal+n_azimuthal];
			}
		}


		printf("# total mass: %g \n", total_mass);

		// print legend
		printf("# mode\tvalue\n");

		for (unsigned int mode = 1; mode <= max_mode; ++mode) {
			double mean_amplitude = 0;

			for (unsigned int n_radial = 1; n_radial < N_radial-1; ++n_radial) {
				double ring_mass = 0;
				for (unsigned int n_azimuthal = 0; n_azimuthal < N_azimuthal; ++n_azimuthal) {
					ring_mass += surface[n_radial]*grid[n_radial*N_azimuthal+n_azimuthal];
				}
				mean_amplitude += ring_mass * amplitude[n_radial*(max_mode+1)+mode];
			}
			mean_amplitude /= total_mass;

			printf("%u\t%.20g\n",mode,mean_amplitude);
		}

		delete [] surface;
	} else if (argv[8][0] == 'r') {
		// print legend
		printf("# radius");
		for (unsigned int mode = 0; mode <= max_mode; ++mode) {
			printf("\tmode %u",mode);
		}
		printf("\n");

		for (unsigned int n_radial = 1; n_radial <= N_radial-1; ++n_radial) {
			// print radii
			printf("%.20g\t",radii[n_radial]);

			// normalize amplitudes
			double total = 0;
			for (unsigned int mode = 0; mode <= max_mode; ++mode) {
				total += pow(amplitude[n_radial*(max_mode+1)+mode],2);
			}
			total = sqrt(total);

			// print mode values
			for (unsigned int mode = 0; mode <= max_mode; ++mode) {
				printf("%.20g\t",amplitude[n_radial*(max_mode+1)+mode]/total);
			}
			printf("\n");
		}
	} else {
		printf("# radii\tmodes\n");
		for (unsigned int n_radial = 1; n_radial < N_radial-1; ++n_radial) {
			printf("%.20g\t", radii[n_radial]);
			for (unsigned int mode = 1; mode <= max_mode; ++mode) {
				printf("%.20g\t",amplitude[n_radial*(max_mode+1)+mode]);
			}
			printf("\n");
		}
	}

	delete [] phase;
	delete [] amplitude;
	delete [] grid;
	delete [] radii;

	return EXIT_SUCCESS;
}
