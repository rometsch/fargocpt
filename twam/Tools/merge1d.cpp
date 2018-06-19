#include <stdio.h>
#include <stdlib.h>



int main(int argc, char** argv) {
//	FILE *fd;

        // check parameter count
        if (argc < 3) {
                fprintf(stderr, "Not enough parameters!\n");
                fprintf(stderr, "Usage: %s n_radial output_filename input_filenames\n",argv[0]);
                return EXIT_FAILURE;
        }

	unsigned int number_of_input_files = argc - 3;
	unsigned int N = atoi(argv[1]);

	// create file handles
	FILE **fd_input = (FILE**)calloc(number_of_input_files, sizeof(FILE*));
	FILE *fd_output = NULL;
	unsigned int *columns = (unsigned int*)calloc(number_of_input_files, sizeof(unsigned int));

	printf("Using %u cells in radial direction.\n", N);
	printf("Output file is '%s'.\n", argv[2]);

	// get size of input files
	for (unsigned int file = 0; file < number_of_input_files; ++file) {
		fd_input[file] = fopen(argv[3+file], "r");
		if (fd_input[file] == NULL) {
			fprintf(stderr, "Could not open file '%s'\n", argv[3+file]);
			goto close_files;
		}

		fseek(fd_input[file], 0L, SEEK_END);
		unsigned int size = ftell(fd_input[file]);
		fseek(fd_input[file], 0L, SEEK_SET);

		columns[file] = size/N/sizeof(double);

		printf("Input file '%s' has %u columns.\n", argv[3+file], columns[file]);
	}

	fd_output = fopen(argv[2], "w");
	if (fd_output == NULL) {
		fprintf(stderr, "Could not open file '%s'\n", argv[2]);
		goto close_files;
	}

	for (unsigned int n = 0; n < N; ++n) {
		for (unsigned int file = 0; file < number_of_input_files; ++file) {
			for (unsigned int col = 0; col < columns[file]; ++col) {
				double data;

				if (fread(&data, sizeof(double), 1, fd_input[file]) < 1) {
					fprintf(stderr, "Could not read from '%s'\n", argv[3+file]);
					goto close_files;
				}

				if (fwrite(&data, sizeof(double), 1, fd_output) < 1) {
					fprintf(stderr, "Could not write to '%s'\n", argv[2]);
					goto close_files;
				}
			}
		}
	}

close_files:
	// close all files
	for (unsigned int file = 0; file < number_of_input_files; ++file) {
		if (fd_input[file] != NULL) {
			fclose(fd_input[file]);
		}
	}

	if (fd_output != NULL) {
		fclose(fd_output);
	}

	free(fd_input);

	return EXIT_SUCCESS;
}
