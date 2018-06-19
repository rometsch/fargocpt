#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
	FILE *fd;
	double buffer[4];

	// check parameter count
	if (argc != 2) {
		fprintf(stderr, "Not enough parameters!\n");
	}

	// open file
	fd = fopen(argv[1], "r");
	if (fd == NULL) {
		fprintf(stderr, "Error while opening '%s'\n", argv[1]);
		return EXIT_FAILURE;
	}

	while (feof(fd) == 0) {
		if (fread(buffer, sizeof(double), 4, fd) == 0) {
			fprintf(stderr, "Error while reading '%s'\n", argv[1]);
			goto clean_up;
		}
		printf("%20g\t%20g\n",buffer[0], buffer[1]);
	}

clean_up:
	// close file
	fclose (fd);

	return EXIT_SUCCESS;
}
