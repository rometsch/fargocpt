#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
	FILE *fd;

	if (argc < 6) {
		fprintf(stderr, "Not enough parameters!\n");
		fprintf(stderr, "Usage: %s filename n_radial n_azimuthal profile_type cell_number\n", argv[0]);
		return EXIT_FAILURE;
	}

	unsigned int N_radial = atoi(argv[2]);
	unsigned int N_azimuthal = atoi(argv[3]);
	unsigned int cell_number = atoi(argv[5]);

	// read in
	unsigned int grid_size = N_radial*N_azimuthal;

        fprintf(stdout, "# Using %ux%u grid.\n",N_radial, N_azimuthal);

	double* grid = new double[grid_size];

        if ((fd = fopen(argv[1],"r")) == NULL) {
                fprintf(stderr, "Could not open file '%s'!\n", argv[1]);
        }

        if (fread(grid, sizeof(double), grid_size, fd) == 0) {
                fprintf(stderr, "Error while reading '%s'!\n", argv[1]);
        }
        fclose(fd);

	// azimuthal profile
	if (argv[4][0] == 'a') {
		for (unsigned int n_azimuthal = 0; n_azimuthal < N_azimuthal; n_azimuthal++) {
			fprintf(stdout, "%u %.16lg\n", n_azimuthal, grid[cell_number*N_azimuthal+n_azimuthal]);
		}

	}

	return EXIT_SUCCESS;
}
