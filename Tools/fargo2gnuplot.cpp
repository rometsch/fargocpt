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
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>

int main(int argc, char* argv[]) {
	FILE *fd;

	// check parameter count
	if (argc < 10) {
		fprintf(stderr, "Not enough parameters!\n");
		fprintf(stderr, "Usage: %s input_filename n_radial n_azimuthal loggrid rmin rmax output_filename out_n_radial out_n_azimuthal\n",argv[0]);
		fprintf(stderr, "  loggrid should be 'l' for logarithmic grid\n");
		return EXIT_FAILURE;
	}

	unsigned int N_radial = atoi(argv[2]);
	unsigned int N_azimuthal = atoi(argv[3]);
	double r_min = atof(argv[5]);
	double r_max = atof(argv[6]);
	unsigned int grid_size = N_radial*N_azimuthal;
	unsigned int N_radial_out = atoi(argv[8]);
	unsigned int N_azimuthal_out = atoi(argv[9]);

	if (N_radial_out > N_radial)
		N_radial_out = N_radial;
	if (N_azimuthal_out > N_azimuthal)
		N_azimuthal_out = N_azimuthal;

	fprintf(stdout, "# Input file '%s' has %ux%u %s grid with r_min = %f and r_max = %f.\n",argv[1], N_radial, N_azimuthal, (argv[4][0] == 'l') ? "logarithmic" : "arithmetic" ,r_min, r_max);
	fprintf(stdout, "# Output file '%s' has %ux%u grid.\n", argv[7], N_radial_out, N_azimuthal_out);

	double* radii = new double[(N_radial+1)*sizeof(double)];

	if (argv[4][0] == 'l')  {
		for (unsigned int n_radial = 0; n_radial <= N_radial; ++n_radial) {
			radii[n_radial] = r_min*exp((double)(n_radial-1.0)/(double)(N_radial-2.0)*log(r_max/r_min));
		}
	} else {
		double interval = (r_max-r_min)/(double)(N_radial-2.0);
		for (unsigned int n_radial = 0; n_radial <= N_radial; ++n_radial) {
			radii[n_radial] = r_min+interval*(double)(n_radial-1.0);
		}
	}

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

	if ((fd = fopen(argv[7],"w")) == NULL) {
		fprintf(stderr, "Could not open file '%s'!\n", argv[7]);
	}

	if ((N_radial_out != N_radial) || (N_azimuthal_out != N_azimuthal)) {
		double dataset[3];

		for (unsigned int n_radial = 0; n_radial < N_radial_out; ++n_radial) {
			// calculate new radius
			if (argv[4][0] == 'l') {
				dataset[1] = r_min*exp((double)(n_radial-1.0)/(double)(N_radial_out-2.0)*log(r_max/r_min));
			} else {
				double interval = (r_max-r_min)/(double)(N_radial_out-1.0);
				dataset[1] = r_min+interval*(double)(n_radial);
			}

			unsigned i = 0;
			if (n_radial < N_radial_out-1) {
				while (radii[i]<=dataset[1]) ++i;
			} else {
				i = N_radial-1;
			}

			for (unsigned int n_azimuthal = 0; n_azimuthal < N_azimuthal_out; ++n_azimuthal) {
				// calculate phi
				dataset[0] = 2*M_PI*(double)n_azimuthal/(double)N_azimuthal_out;

				unsigned int j = n_azimuthal*N_azimuthal/N_azimuthal_out;

				dataset[2] = grid[i*N_azimuthal+j];
				fwrite(dataset, sizeof(dataset), 1, fd);
			}

			// write first again
			dataset[0] = 0;
			dataset[2] =  grid[i*N_azimuthal+0];
			fwrite(dataset, sizeof(dataset), 1, fd);
		}
	} else {
		double dataset[3];

		for (unsigned int n_radial = 0; n_radial < N_radial_out; ++n_radial) {
			dataset[1] = radii[n_radial];
			for (unsigned int n_azimuthal = 0; n_azimuthal < N_azimuthal_out; ++n_azimuthal) {
				dataset[0] = 2*M_PI*(double)n_azimuthal/(double)N_azimuthal_out;
				dataset[2] = grid[n_radial*N_azimuthal+n_azimuthal];
				fwrite(dataset, sizeof(dataset), 1, fd);
			}
			dataset[0] = 2*M_PI*(double)0/(double)N_azimuthal_out;
			dataset[2] = grid[n_radial*N_azimuthal+0];
			fwrite(dataset, sizeof(dataset), 1, fd);			
		}
	}

	fprintf(stdout, "splot \"%s\" binary record=%ux%u format=\"%%double\" u ($2*cos($1)):($2*sin($1)):3\n", argv[7], N_azimuthal_out+1, N_radial_out);

	fclose(fd);

	delete [] grid;
	delete [] radii;

	return EXIT_SUCCESS;
}
