#include <stdio.h>
#include <stdlib.h>

#include <hdf5_read.h>

int main (int argcount, char** argvektor) {

	// basic variables
	double* config;

	// check whether user provided a filename and timestep to compute
	if (argcount != 3) {
		fprintf(stderr,    "Please provide a filename (1) and timestep (2) (number, \"first\" or \"last\") to read.\n
							Program will now terminate. \n");
		return EXIT_SUCCESS;
	}

	// read the expected configuration from file
	config = read_config(argvektor[1], argvektor[2]);

}