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

// variables for initiating the read parameter
static char* file;


/*-------------------------------------------------------------------------------------------------------*/
int init_file(char* infile, int no_writeouts, int N) {
	file = infile;

	// Open file, check whether operation was successful
	FILE *outfile = fopen(file, "a+");
	if (outfile == NULL) {
		fprintf(stderr, "The file %s could not be opened on initialization.\n", file);
		return EXIT_FAILURE;
	}

	// Upon first time setup, print N and amount of expected writeouts, no check is made, whether the data is actually written
	no_writeouts ++;
	fwrite(&N, sizeof(int), 1, outfile);
	fwrite(&no_writeouts, sizeof(int), 1, outfile);

	// Close outgoing stream
	fclose(outfile);
	fprintf(stderr, "Completed initializing outgoing file.\n");

	return EXIT_SUCCESS;
}

/*-------------------------------------------------------------------------------------------------------*/
int write_file(double *position, int timestep, int N) {

	// Open file, check whether operation was successful
	FILE *outfile = fopen(file, "a+");
	if (outfile == NULL) {
		fprintf(stderr, "The file %s could not be opened on timestep: %d\n", file, timestep);
		return EXIT_FAILURE;
	}

	// Write Parameters to file, no check will be made whether the data is actually written
	fwrite(&timestep, sizeof(int), 1, outfile);
	fwrite(position, sizeof(double), 2*N, outfile);


	// Close outgoing stream
	fclose(outfile);
	fprintf(stderr, "Completed writing file on timestep: %d\n", timestep);

	return EXIT_SUCCESS;
}

/*-------------------------------------------------------------------------------------------------------*/
/* Template for config file:
 * Number of particles (N)              	: 200
 * Energie (k_B*T)				         	: 1
 * Dipole interaction relation (Gamma_A)  	: 50
 * Dipole interaction relation (Gamma_B)  	: 20
 * Shear rate A (shear_A)               	: 1000
 * Shear rate B (shear_B)               	: -1000
 * Brownian Diffusion Time (tau_B)       	: 1
 * Brownian Diffusion (D_Brown_A)	       	: 1
 * Time difference (delta_t)            	: 1e-6
 * Destination file (outfile)           	: /great_file_name.bin
 * Iterations                           	: 100000
 * Number of Threads                    	: 4
 * Writeouts                            	: 1000
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

	if ((check = fscanf(infile, "Number of particles (N)              	: %d\n",  (&(param->N)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Energie (k_B*T)				      	: %lf\n", (&(param->kT)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Dipole interaction relation (Gamma_A)	: %lf\n", (&(param->Gamma_A)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Dipole interaction relation (Gamma_B)	: %lf\n", (&(param->Gamma_B)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Shear rate A (shear_A)               	: %lf\n", (&(param->shear_A)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Shear rate B (shear_B)               	: %lf\n", (&(param->shear_B)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Brownian Diffusion Time (tau_B)       	: %lf\n", (&(param->tau_B)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Brownian Diffusion (D_Brown_A)	       	: %lf\n", (&(param->D_Brown_A)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Time difference (delta_t)            	: %lf\n", (&(param->timestep)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Destination file (outfile)           	: %s\n",  temp)) < 1)
		return NULL;
	if ((check = fscanf(infile, "Iterations                           	: %d\n",  (&(param->max_timesteps)))) < 1)
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
