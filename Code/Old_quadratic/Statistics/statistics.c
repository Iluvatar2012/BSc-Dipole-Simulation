#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <math.h>

#include "structs.h"
#include "hdf5.h"
#include "utilities.h"


/*----------------------------------------------------------------------------------------------------------------------------*/
int main(int argcount, char** argvector) {

	// basic variables
	char 	file[1024];
	double 	t;
	int 	steps, N, bands;

	// check whether the right amount of arguments is given
	if (argcount != 2) {
		fprintf(stderr, "Please provide only a filename to read from, process will terminate.\n");
		return EXIT_FAILURE;
	}

	// copy the filename, we will need it
	strncpy(file, argvector[1], 1024);

	// read from the given file, check whether successful
	struct parameters *param = hdf5_read(file);

	if (param == NULL) {
		return EXIT_FAILURE;
	}

	// copy simulation attributes to variables
	N 		= param->N;
	steps 	= param->steps;
	bands	= ceil(2*sqrt(N/2.));

	// calculate mean values and standard deviations with step dependence
	struct results *res  = stat_analysis(param);

	// open new files for all parameters
	FILE* stat_lane = fopen("stat_lane", "w+");
	FILE* stat_psi4 = fopen("stat_psi4", "w+");
	FILE* stat_psi6 = fopen("stat_psi6", "w+");

	// print all parameters to respective files
	for (int i=0; i<=steps; i+=10) {
		t = i*1000*1e-6;


		for (int j=0; j<bands; j++) {
			fprintf(stat_lane, "%lf\t%lf\t%lf\n", t, res->y_lane[i*bands + j], res->max_lane[i*bands+j]);
		}
		
		// fprintf(stat_lane, "%lf\t%lf\t%lf\t%lf\t%lf\n", t, res->mean_lane[i], res->stddev_lane[i], res->max_lane[i], res->y_lane[i]);
		fprintf(stat_psi4, "%lf\t%lf\t%lf\n", t, res->mean_psi4[i], res->stddev_psi4[i]);
		fprintf(stat_psi6, "%lf\t%lf\t%lf\n", t, res->mean_psi6[i], res->stddev_psi6[i]);
	}

	// close all files
	fclose(stat_lane);
	fclose(stat_psi4);
	fclose(stat_psi6);

	// free memory
	free(param);
	free(res);

	// return to caller
	return EXIT_SUCCESS;
}