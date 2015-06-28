#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "structs.h"
#include "hdf5.h"
#include "picture.h"

/*----------------------------------------------------------------------------------------------------------------------------*/
int main(int argcount, char** argvector) {

	// basic variables
	char file[1024]; 
	char outfile[1024];
	char buf[32];

	int steps, counter;
	int stepsize;
	int init_image_num;

	double perc;

	// check whether the right amount of arguments is given
	if (argcount != 4) {
		fprintf(stderr, "Please provide a filename to read from, how many steps to skip (one is none) and the initial images number, process will terminate.\n");
		return EXIT_FAILURE;
	}

	// copy the filename, stepsize and the first images number
	strncpy(file, argvector[1], 1024);
	stepsize  		= strtol(argvector[2], NULL, 10);
	init_image_num	= strtol(argvector[3], NULL, 10);

	// user output
	fprintf(stderr, "Reading file: %s\n", file);

	// read from the given file, check whether successful
	struct parameters *param = hdf5_init(file);

	if (param == NULL) {
		return EXIT_FAILURE;
	}

	// copy the number of steps we have
	steps = param->steps;

	// initiate the SDL libraries, check if successful
	if (initiate(param) != EXIT_SUCCESS) 
		return EXIT_FAILURE;

	// user output
	fprintf(stderr, "Computing...\n");

	// initiate counter
	counter = init_image_num;

	// iterate over all steps, compute filenames and make pictures
	for (int i=0; i<=steps; i+=stepsize) {
		// compute percentage of completed iterations
		perc = 100.*i/steps;

		// output percentage to user
		fprintf(stderr, "Progress: \t\t[");

		for (int j=0; j<floor(perc); j++) {
			fprintf(stderr, "=");
		}

		if (perc != 100)
			fprintf(stderr, ">");

		for (int j=floor(perc)+1; j<100; j++) {
			fprintf(stderr, ".");
		}

		fprintf(stderr, "] \t%.1lf%%\r", perc);

		// compute the next filename
		strncpy(outfile, "/home/ayra/.images/image_", 984);

		sprintf(buf, "%d", counter);
		strncat(outfile, buf, 32);
		strncat(outfile, ".bmp\0", 5);

		// update position and psi4 array to the next step
		hdf5_read(i);

		// build the picture
		draw_picture(i, outfile);

		// increase counter
		counter++;
	}

	// user output
	fprintf(stderr, "\nDone!\n\n");
	fprintf(stdout, "%d\n", counter);

	// free memory used by SDL
	destroy();

	// free struct, return to caller
	free(param->positions);
	free(param);

	return EXIT_SUCCESS;
}