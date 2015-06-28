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

#include "structs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*-------------------------------------------------------------------------------------------------------*/
/* Template for config file:
	Dipole interaction relation (Gamma_A)	: 10
	Particle dipole ratio (m)				: 1.0
	Shear rate (gamma)		               	: 20
	Particle Diffusion ratio (D_rat)		: 1.7
	Time difference (delta_t)            	: 0.000001
	Iteration time (tau)                  	: 10
	Writeout step                           : 1000
	Destination file (outfile)           	: Results/default.hdf5
 
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

	if ((check = fscanf(infile, "Dipole interaction relation (Gamma_A)	: %lf\n", (&(param->Gamma_A)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Particle dipole ratio (m)				: %lf\n", (&(param->m)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Shear rate (gamma)		               	: %lf\n", (&(param->gamma_shear)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Particle Diffusion ratio (D_rat)		: %lf\n", (&(param->D_rat)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Time difference (delta_t)            	: %lf\n", (&(param->timestep)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Iteration time (tau)                  	: %d\n",  (&(param->tau)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Writeout step                          : %d\n",    (&(param->write_step)))) < 1)
		return NULL;
	if ((check = fscanf(infile, "Destination file (outfile)           	: %s\n",  temp)) < 1)
		return NULL;

	// Concatenate current working directory with filename from config file, close the file stream
	strncpy((param->outfile), temp, 1024);
	fclose(infile);

	// return to caller
	return param;
}
/*-------------------------------------------------------------------------------------------------------*/
