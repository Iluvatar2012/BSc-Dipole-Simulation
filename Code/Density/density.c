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

	int 	N, N_A;
	int 	step_count;

	double	count_A_in, count_A_out;
	double 	count_B_in, count_B_out;
	double 	L_one, L_two;

	double* positions;


	// check if the right number of arguments ist given
	if (argcount != 3) {
		fprintf(stderr, "Please provide a filename to read from and one to append to. \n");
		return EXIT_FAILURE;
	}

	// read the filenames and binsize
	strncpy(infile, argvec[1], 1024);
	strncpy(out_filename, argvec[2], 1024);

	// initiate counting values
	count_A_in 	= 0;
	count_A_out = 0;
	count_B_in 	= 0;
	count_B_out = 0;

	step_count	= 0;

	// read data from file
	struct parameters *param = hdf5_init(infile);

	// calculate fraction of A particles and the cutoff sizes
	N 		= param->N;
	N_A 	= round(N * param->X);

	L_one = 0.25*param->L_y;
	L_two = 0.75*param->L_y;

	// There is no point in doing this analysis if we exclusively have A or B particles
	if ((N_A == 0) || (N_A == N))
		return EXIT_FAILURE;
	
	// copy the current positions
	positions = param->positions;

	// iterate over the last half of the steps
	for (int step=round(param->last_step*1/3.)+1; step<=param->last_step; step++) {

		// read the current step, increment step counter (done so nobody fucks up in determining the amount of steps read)
		hdf5_read(step);
		step_count++;

		// iterate over all A particles
		for (int i=0; i<N_A; i++) {
			// check whether the current particle is on the in- or outside
			if ((positions[2*i+1] >= L_one) && (positions[2*i+1] <= L_two))
				count_A_in++;
			else
				count_A_out++;
		}	

		// iterate over all B particles
		for (int i=N_A; i<N; i++) {
			// check whether the current particle is on the in- or outside
			if ((positions[2*i+1] >= L_one) && (positions[2*i+1] <= L_two))
				count_B_in++;
			else
				count_B_out++;
		}
	}

	// normalize the values to one step
	count_A_in 	/= step_count;
	count_A_out /= step_count;
	count_B_in 	/= step_count;
	count_B_out /= step_count;


	// open filestream with overwrite attribute and output the simulations basic data
	FILE *outfile = fopen(out_filename, "a+");

	fprintf(outfile, "%lf,\t%lf,\t%lf,\t%lf,\t%lf,\t%lf,\t%lf,\t%lf,\t%lf,\t%lf,\t%lf\n", 	param->gamma, param->G, param->X, \
																							count_A_in, count_A_out, count_B_in, count_B_out, \
																							count_A_in/N_A, count_A_out/N_A, count_B_in/(N-N_A), count_B_out/(N-N_A));

	// close the stream and delete the struct
	fclose(outfile);

	free(positions);
	free(param);

	// return to caller
	return EXIT_SUCCESS;
}