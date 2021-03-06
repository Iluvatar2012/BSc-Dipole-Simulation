#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "structs.h"
#include "hdf5.h"

/*----------------------------------------------------------------------------------------------------------------------------*/
int main(int argcount, char** argvec) {

	// basic variables
	char	infile[1024];
	char 	out_filename[1024];

	double 	timestep, t;
	int 	steps, step_size, N;

	// check if the right number of arguments ist given
	if (argcount != 5) {
		fprintf(stderr, "Please provide a filename, a filename to write to, the timestep of the simulation and the amount of steps to skip. \n");
		return EXIT_FAILURE;
	}

	// read the filenames, timestep of the simulation and steps between writes from the console
	strncpy(infile, argvec[1], 1024);
	strncpy(out_filename, argvec[2], 1024);

	timestep 	= atof(argvec[3]);
	step_size	= atoi(argvec[4]);


	// output to user
	fprintf(stdout, "Reading file...\n");
	fflush(stdout);

	// read data from file
	struct parameters *param = hdf5_init(infile);

	// output to user
	fprintf(stdout, "Finished initiating file, there are %d steps for %d particles.\n", param->steps, param->N);

	// open filestream with overwrite attribute and output the simulations basic data
	FILE *outfile = fopen(out_filename, "w+");

	fprintf(outfile, "%d # N\n", param->N);
	fprintf(outfile, "%d # steps\n\n", param->steps);

	// read the amount of particles and steps from the struct
	N 		= param->N;
	steps 	= param->steps;

	// print a mask explaining the significance of each value
	fprintf(outfile, "#time, x, y, disp_x, disp_y, psi_4, psi_6, laning\n\n");

	// iterate over all wanted times
	for (int i=0; i<=steps; i+=step_size) {

		// iterate current time
		t = i*timestep;

		// get the current step
		hdf5_read(i);

		// iterate over all particles
		for (int j=0; j<N; j++) {
			fprintf(outfile, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", t, param->positions[2*j], param->positions[2*j+1], \
					param->displacement[2*j], param->displacement[2*j+1], param->psi4[j], param->psi6[j], param->laning[j]);
		}
	}

	// close the stream and delete the struct
	fclose(outfile);

	free(param);

	// return to caller
	return EXIT_SUCCESS;
}