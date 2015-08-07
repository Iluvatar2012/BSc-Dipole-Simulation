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
	char 	out_filename[1024];

	double 	binsize, binmin, binmax, bincen;
	int 	bins, N, N_A;

	int 	tot_count_A, tot_count_B;

	int* 	count_A;
	int*	count_B;
	double* positions;


	// check if the right number of arguments ist given
	if (argcount != 4) {
		fprintf(stderr, "Please provide a filename, a filename to write to and the amount of bins. \n");
		return EXIT_FAILURE;
	}

	// read the filenames and binsize
	strncpy(infile, argvec[1], 1024);
	strncpy(out_filename, argvec[2], 1024);

	bins = atoi(argvec[3]);

	// allocate memory for the counting array, initiate every value to zero
	count_A = malloc(bins*sizeof(int));
	count_B = malloc(bins*sizeof(int));

	for (int i=0; i<bins; i++) {
		count_A[i] = 0;
		count_B[i] = 0;
	}

	// read data from file
	struct parameters *param = hdf5_init(infile);

	// calculate fraction of A particles and the binsize
	N 		= param->N;
	N_A 	= round(N * param->X);
	binsize = param->L_y / bins;

	// copy the current positions
	positions = param->positions;

	// iterate over the last half of the steps
	for (int step=round(param->last_step*1/3.)+1; step<=param->last_step; step++) {

		// read the current step
		hdf5_read(step);

		// iterate over all bins
		for (int bin=0; bin<bins; bin++) {

			// calculate the current minimum and maximum binning values
			binmin = bin*binsize;
			binmax = (bin + 1)*binsize;

			// iterate over all A particles
			for (int i=0; i<N_A; i++) {
				// count whether the current particle is in the current bin
				if ((positions[2*i+1] >= binmin) && (positions[2*i+1] < binmax))
					count_A[bin]++;
			}	

			// iterate over all B particles
			for (int i=N_A; i<N; i++) {
				// count whether the current particle is in the current bin
				if ((positions[2*i+1] >= binmin) && (positions[2*i+1] < binmax))
					count_B[bin]++;
			}
		}
	}

	// open filestream with overwrite attribute and output the simulations basic data
	FILE *outfile = fopen(out_filename, "w+");

	// count the total number of counts through all bins
	tot_count_A = 0;
	tot_count_B = 0;

	for (int i=0; i<bins; i++) {
		tot_count_A += count_A[i];
		tot_count_B += count_B[i];
	}

	// iterate over all bins, output to user
	for (int i=0; i<bins; i++) {
		
		// calculate the center of the bin, output the binning values to file
		bincen = (0.5+i) * binsize - 0.5*param->L_y;
		fprintf(outfile, "%lf \t%lf \t%lf\n", bincen, count_A[i]*1./tot_count_A, count_B[i]*1./tot_count_B);
	}

	// close the stream and delete the struct
	fclose(outfile);

	free(positions);
	free(param);

	// return to caller
	return EXIT_SUCCESS;
}