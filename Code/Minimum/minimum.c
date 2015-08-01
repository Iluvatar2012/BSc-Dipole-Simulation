#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "structs.h"
#include "hdf5.h"

/*----------------------------------------------------------------------------------------------------------------------------*/
int main(int argcount, char** argvec) {

	// basic variables
	char	infile[1024];
	int 	N;

	double* positions;

	double 	min = 1.0;


	// check if the right number of arguments ist given
	if (argcount != 2) {
		fprintf(stderr, "Please provide a filename to read from. \n");
		return EXIT_FAILURE;
	}

	// read the filenames and binsize
	strncpy(infile, argvec[1], 1024);

	// read data from file
	struct parameters *param = hdf5_init(infile);

	// calculate fraction of A particles and the cutoff sizes
	N 		= param->N;
	
	// copy the current positions
	positions = param->positions;

	// iterate over the last half of the steps
	for (int step=round(param->last_step/2.)+1; step<=param->last_step; step++) {

		// read the current step, increment step counter (done so nobody fucks up in determining the amount of steps read)
		hdf5_read(step);

		double loc_min = param->L_y;
		double dx, dy, r;

		for (int i=0; i<N; i++) {
			for(int j=0; j<N; j++) {

				if (i==j)
					continue;

				dx = positions[2*i] 	- positions[2*j];
				dy = positions[2*i+1] 	- positions[2*j+1];

				r  = sqrt(dx*dx + dy*dy);

				if (r < loc_min)
					loc_min = r;
			}
		}

		if (loc_min < min)
			min = loc_min;
		
	}

	fprintf(stdout, "%lf\n", min);

	free(positions);
	free(param);

	// return to caller
	return EXIT_SUCCESS;
}