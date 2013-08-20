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
int read_file(char* file, char* outfile, int* N, double* kT, double* Gamma, double* shear, \
		double* tau_B, double* D_Brown, double* timestep, int* max_timesteps, int* thread_number, \
		int* no_writeouts) {
	// open configuration file, check whether operation was successful
	FILE *infile = fopen(file, "r");
	if(infile == NULL) {
		fprintf(stderr, "The file \"%s\" could not be found or opened.\n", file);
		return EXIT_FAILURE;
	}

	// variable for checking whether everything went well
	int check;

	// read input from document, adjust output location and return
	char temp[1024];

	if ((check = fscanf(infile, "Number of particles (N)              : %d\n",  N)) < 1)
		return EXIT_FAILURE;
	if ((check = fscanf(infile, "Energie (k_B*T)				      : %lf\n", kT)) < 1)
		return EXIT_FAILURE;
	if ((check = fscanf(infile, "Dipole interaction relation (Gamma)  : %lf\n", Gamma)) < 1)
		return EXIT_FAILURE;
	if ((check = fscanf(infile, "Shear rate (shear)                   : %lf\n", shear)) < 1)
		return EXIT_FAILURE;
	if ((check = fscanf(infile, "Brownian Diffusion Time	          : %lf\n", tau_B)) < 1)
		return EXIT_FAILURE;
	if ((check = fscanf(infile, "Brownian Diffusion		              : %lf\n", D_Brown)) < 1)
		return EXIT_FAILURE;
	if ((check = fscanf(infile, "Time difference (delta_t)            : %lf\n", timestep)) < 1)
		return EXIT_FAILURE;
	if ((check = fscanf(infile, "Destination file (outfile)           : %s\n",  temp)) < 1)
		return EXIT_FAILURE;
	if ((check = fscanf(infile, "Iterations                           : %d\n",  max_timesteps)) < 1)
		return EXIT_FAILURE;
	if ((check = fscanf(infile, "Number of Threads                    : %d\n",  thread_number)) < 1)
		return EXIT_FAILURE;
	if ((check = fscanf(infile, "Writeouts                            : %d",    no_writeouts)) < 1)
		return EXIT_FAILURE;

	// Concatenate current working directory with filename from config file
	strncat(outfile, temp, 1024-strlen(outfile));
	fclose(infile);

	return EXIT_SUCCESS;
}
