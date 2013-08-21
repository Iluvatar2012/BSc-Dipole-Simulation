/*
 * fileIO_v1.c
 * 		//TODO: update this
 *		This C-file, consisting only of the function write_file(), opens a file and writes given arrays in
 *		the format seen below into the file. If the file does not exist, a new file, with the name given in
 *		the functions variables, will be created. All data is added to the end of original file.
 *
 *
 *					read_file():
 *					 Reads the file, given as argument to the function and stores the given variables in
 *					 the static variables declared at the start of the program.
 *					WARNING:
 *					 The function requires a specific template, as seen in commentary just before the
 *					 function itself. The function checks whether the incoming file conforms to the specified
 *					 template and exits if it does not, so a premature closure of the program may be due to
 *					 a false file
 *
 *  	Last Changed: 	Aug 1, 2013
 *      Author: 		Aiko Bernehed
 */

#include "fileIO_v1.h"
#include "structs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*-------------------------------------------------------------------------------------------------------*/
int write_file(char *file, double *position, int timestep, int no_writeouts, int N) {

	// Open file, check whether operation was successful
	FILE *outfile = fopen(file, "a+");
	if (outfile == NULL) {
		fprintf(stderr, "The file %s could not be opened on timestep: %d\n", file, timestep);
		return EXIT_FAILURE;
	}

	// Upon first time setup, print N and amount of expected writeouts, no check is made, whether the data is actually written
	if (timestep == 0) {
		no_writeouts ++;
		fwrite(&N, sizeof(int), 1, outfile);
		fwrite(&no_writeouts, sizeof(int), 1, outfile);
	}

	// Write Parameters to file
	fwrite(&timestep, sizeof(int), 1, outfile);
	fwrite(position, sizeof(double), 2*N, outfile);


	// Close outgoing stream
	fclose(outfile);
	fprintf(stderr, "Completed writing file on timestep: %d\n", timestep);

	return EXIT_SUCCESS;
}

/*-------------------------------------------------------------------------------------------------------*/
/* Template for config file:
 * Number of particles (N)              : 1000
 * Energie (k_B*T)				        : 1
 * Dipole interaction relation (Gamma)  : 1
 * Shear rate (shear)                   : 1
 * Brownian Diffusion Time	            : 1
 * Brownian Diffusion		            : 1
 * Time difference (delta_t)            : 2e-7
 * Destination file (outfile)           : great_file_name.txt
 * Iterations                           : 100000
 * Number of Threads                    : 8
 * Writeouts                            : 10
 *
 * */
struct parameters *read_file(char* file) {
	// open configuration file, check whether operation was successful
	FILE *infile = fopen(file, "r");
	if(infile == NULL) {
		fprintf(stderr, "The file \"%s\" could not be found or opened.\n", file);
		return NULL;
	}

	// define a new struct which we will return
	struct parameters *param = malloc(sizeof(struct parameters));

	// variable for checking whether everything went well
	int check;

	// read input from document, adjust output location and return
	char temp[1024];

	if ((check = fscanf(infile, "Number of particles (N)              : %d\n",  (&(param->N)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Energie (k_B*T)				      : %lf\n", (&(param->kT)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Dipole interaction relation (Gamma)  : %lf\n", (&(param->Gamma)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Shear rate A (shear_A)               : %lf\n", (&(param->shear_A)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Shear rate B (shear_B)               : %lf\n", (&(param->shear_B)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Brownian Diffusion Time	          : %lf\n", (&(param->tau_B)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Brownian Diffusion		              : %lf\n", (&(param->D_Brown)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Time difference (delta_t)            : %lf\n", (&(param->timestep)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Destination file (outfile)           : %s\n",  temp)) < 1)
		return NULL;
	if ((check = fscanf(infile, "Iterations                           : %d\n",  (&(param->max_timesteps)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Number of Threads                    : %d\n",  (&(param->thread_number)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Writeouts                            : %d",    (&(param->no_writeouts)))) < 1)
		return NULL;

	// Concatenate current working directory with filename from config file, close the file stream
	strncat((param->outfile), temp, 1024-strlen((param->outfile)));
	fclose(infile);

	// return to caller
	return param;
}
/*-------------------------------------------------------------------------------------------------------*/
