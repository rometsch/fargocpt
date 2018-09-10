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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char* argv[])
{
	FILE *fd;
	unsigned int N_radial, N_azimuthal, grid_size;

	// check parameter count
	if (argc < 8) {
		fprintf(stderr, "Not enough parameters!\n");
		fprintf(stderr, "Usage: %s filename n_radial n_azimuthal loggrid rmin rmax mode [mode_options]\n",argv[0]);
		return EXIT_FAILURE;
	}

	N_radial = atoi(argv[2]);
	N_azimuthal = atoi(argv[3]);
	double r_min = atof(argv[5]);
	double r_max = atof(argv[6]);
	grid_size = N_radial*N_azimuthal;

	fprintf(stdout, "Using %ux%u grid.\n",N_radial, N_azimuthal);
	fprintf(stdout, "r_min = %f\n",r_min);
	fprintf(stdout, "r_max = %f\n",r_max);

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

	double* grid = new double[grid_size];

	switch (argv[7][0]) {
		case 'p': {
			fprintf(stdout, "Using mode 'Gauss-Point'.\n");
			if (argc < 11) {
				fprintf(stderr, "Not enough parameters!\n");
				fprintf(stderr, "mode_options are: r_0 phi_0 sigma y_0\n");
				return EXIT_FAILURE;
			}

			double r_0 = atof(argv[8]);
			double phi_0 = atof(argv[9]);
			double sigma = atof(argv[10]);
			double y_0 = atof(argv[11]);

			fprintf(stdout, "r_0 = %f\n",r_0);
			fprintf(stdout, "phi_0 = %f\n",phi_0);
			fprintf(stdout, "sigma = %f\n",sigma);
			fprintf(stdout, "y_0 = %f\n",y_0);
			for (unsigned int n_radial = 0; n_radial < N_radial; ++n_radial) {
				for (unsigned int n_azimuthal = 0; n_azimuthal < N_azimuthal; ++n_azimuthal) {
					unsigned int cell = (n_radial*N_azimuthal)+n_azimuthal;
					double r = radii[n_radial];
					double phi = (double)n_azimuthal*(2.0*M_PI/(double)N_azimuthal);
					grid[cell] = y_0 *exp(-0.5*pow( (  sqrt( pow(r*cos(phi)-r_0*cos(phi_0),2) + pow(r*sin(phi)-r_0*sin(phi_0),2)  )  )/sigma,2));
				}
			}
			}
			break;
		case 'g': {
			fprintf(stdout,"Using mode 'Gauss'.\n");
			if (argc < 10) {
				fprintf(stderr, "Not enough parameters!\n");
				fprintf(stderr, "mode_options are: x_0 sigma y_0\n");
				return EXIT_FAILURE;
			}

			double x_0 = atof(argv[8]);
			double sigma = atof(argv[9]);
			double y_0 = atof(argv[10]);

			fprintf(stdout, "x_0 = %f\n",x_0);
			fprintf(stdout, "sigma = %f\n",sigma);
			fprintf(stdout, "y_0 = %f\n",y_0);
			for (unsigned int n_radial = 0; n_radial < N_radial; ++n_radial) {
				for (unsigned int n_azimuthal = 0; n_azimuthal < N_azimuthal; ++n_azimuthal) {
					unsigned int cell = (n_radial*N_azimuthal)+n_azimuthal;
					grid[cell] = y_0 *exp(-0.5*pow((radii[n_radial]-x_0)/sigma,2));
				}
			}
			}
			break;
		case 'm': {
			fprintf(stdout,"Using mode 'mode'.\n");
			if (argc < 10) {
				fprintf(stderr, "Not enough parameters!\n");
				fprintf(stderr, "mode_options are: exponent mode amplitude\n");
				return EXIT_FAILURE;
			}

			double exp = atof(argv[8]);
			double mode = atof(argv[9]);
			double amplitude = atof(argv[10]);

			fprintf(stdout, "exp = %f\n", exp);
			fprintf(stdout, "mode = %f\n", mode);
			fprintf(stdout, "amplitude = %f\n", amplitude);

			for (unsigned int n_radial = 0; n_radial < N_radial; ++n_radial) {
				for (unsigned int n_azimuthal = 0; n_azimuthal < N_azimuthal; ++n_azimuthal) {
					unsigned int cell = (n_radial*N_azimuthal)+n_azimuthal;
					grid[cell] = pow(radii[n_radial],-exp) + sin(mode*(double)n_azimuthal/(double)N_azimuthal*2.0*M_PI)*amplitude;
				}
			}
			}
			break;
		case 'a': {
			fprintf(stdout,"Using axis 'mode'\n");
			if (argc < 7) {
				fprintf(stderr, "Not enough parameters!\n");
				fprintf(stderr, "mode_options are: none\n");
				return EXIT_FAILURE;
			}

			for (unsigned int n_radial = 0; n_radial < N_radial; ++n_radial) {
				for (unsigned int n_azimuthal = 0; n_azimuthal < N_azimuthal; ++n_azimuthal) {
					unsigned int cell = (n_radial*N_azimuthal)+n_azimuthal;
					grid[cell] = + cos((double)n_azimuthal/(double)N_azimuthal*2.0*M_PI);
				}
			}

			}
			break;
		case 'j': {
			fprintf(stdout,"Using jump 'mode'\n");
			if (argc < 7) {
				fprintf(stderr, "Not enough parameters!\n");
				fprintf(stderr, "mode_options are: none\n");
				return EXIT_FAILURE;
			}

			double Sigma0 = 100.0;
			double MU = 2.35;
			double ADIABATICINDEX = 1.4;
			double R = 8.31447e+07;
			double minT = 10;
			double maxT = 100;

			for (unsigned int n_radial = 0; n_radial < N_radial; ++n_radial) {
				for (unsigned int n_azimuthal = 0; n_azimuthal < N_azimuthal; ++n_azimuthal) {
					unsigned int cell = (n_radial*N_azimuthal)+n_azimuthal;
//					grid[cell] = (n_radial > N_radial/2 ? minT : (n_radial == N_radial/2 ? (minT+maxT)/2.0: maxT))*Sigma0/(ADIABATICINDEX-1.0)/MU*R; //* 1.545511414728764e+20;
					grid[cell] = (n_azimuthal > N_azimuthal/2 ? minT : (n_azimuthal == N_azimuthal/2 ? (minT+maxT)/2.0: maxT))*Sigma0/(ADIABATICINDEX-1.0)/MU*R; //* 1.545511414728764e+20;
				}
			}

			}

			break;

		case 's': {
			fprintf(stdout, "Using mode 'shocktube'\n");

			if (argc < 9) {
				fprintf(stderr, "Not enough parameters!\n");
				fprintf(stderr, "mode_options are: <left value> <right value>\n");
				return EXIT_FAILURE;
			}

			double left = atof(argv[8]);
			double right = atof(argv[9]);

			printf("Left value: %lf\n", left);
			printf("Right value: %lf\n", right);

			for (unsigned int n_radial = 0; n_radial < N_radial; ++n_radial) {
				for (unsigned int n_azimuthal = 0; n_azimuthal < N_azimuthal; ++n_azimuthal) {
					unsigned int cell = (n_radial*N_azimuthal)+n_azimuthal;
					grid[cell] = n_radial < N_radial/2 ? left : right;
				}
			}

			break;

		}

		default:
			fprintf(stderr,"Unknown mode: '%s'\n", argv[7]);
			break;
	}

	// open file
	fd = fopen(argv[1], "w");
	if (fd == NULL) {
		fprintf(stderr, "Error while opening '%s'\n", argv[1]);
		return EXIT_FAILURE;
	}

	fprintf(stdout, "Writing to '%s'.\n", argv[1]);

	fwrite(grid, grid_size, sizeof(double), fd);

	// close file
	fclose (fd);

	delete [] grid;
	delete [] radii;

	return EXIT_SUCCESS;
}
