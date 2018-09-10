#include <stdio.h>
#include <stdlib.h>

#include "../src/particle.h"

int main(int argc, char** argv) {
	if (argc < 4) {
		fprintf(stderr, "Not enough parameters!\n");
		fprintf(stderr, "Usage: %s folder planet_id outputfile\n", argv[0]);
		return EXIT_FAILURE;
	}

	unsigned int id = atoi(argv[2]);
	unsigned int timestep = 0;

	FILE *fd_output = fopen(argv[3], "w");
	if (fd_output == NULL) {
		fprintf(stderr, "Could not open outputfile '%s'!\n", argv[3]);
		return EXIT_FAILURE;
	}

	fprintf(fd_output, "# Particle id = %u\n", id);

	while (1) {
		char *filename;
		if (asprintf(&filename, "%s/particles%i.dat",argv[1],timestep) <0) {
		}

		FILE *fd_input = fopen(filename, "r");
		free(filename);
		if (fd_input == NULL) {
			break;
		}

		// find offset
		t_particle particle;
		bool found = false;

		while (1) {
			fread(&particle, sizeof(t_particle), 1, fd_input);

			if (feof(fd_input))
				break;

			if (particle.id == id) {
				found = true;
				break;
			}
		}

		if (!found) {
			break;
		}

		fprintf(fd_output, "%lf %lf %lf %lf %lf\n", particle.x, particle.y, particle.vx, particle.vy, particle.mass);

		fclose(fd_input);
		timestep++;
	}

	fclose(fd_output);

	printf("Parsed %u files.\n", timestep);

	return EXIT_SUCCESS;
}
