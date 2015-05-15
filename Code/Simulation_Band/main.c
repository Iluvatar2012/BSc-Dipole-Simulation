/* main.c
 *		This program simulates the interaction of induced dipoles in an external magnetic field,
 *		confined to a 2-dimensional plane. All important features of the simulation can be changed
 *		via the config-files, giving the user full control of the system.
 *
 *		USE:		Simulation is started via the exec-file in folder Debug of this project. Parameters
 *					can be changed via the config-files in the folder Config_Files. A standard template
 *					is given in the commentary just before the method "read_file".
 *
 *		REQUIRES:	Common Libraries:
 *					stdio.h
 *					stdlib.h
 *					string.h
 *					math.h
 *					unistd.h
 *					pthread.h
 *
 *					Own Header files:
 *					fileIO_v1.h, fileIO_v1.c
 *					simulation.h
 *
 *		CONTENT:	init():
 *					 Initializes all remaining unset static variables and allocates memory for all arrays.
 *					 The function has return values EXIT_SUCCESS and EXIT_FAILURE, which should be regarded
 *					 due to memory failures.
 *
 *					simulation():
 *					 This is the actual part where an entire timestep is computed. It synchronizes with all
 *					 iterating threads, lets them compute a timestep and finally does some final analysis of the
 *					 step itself (should the current parameters be written to a file, does the verlet list need
 *					 updating, etc.).
 *
 *					iteration():
 *					 As the name suggests, this is the heart of the simulation. It computes, with periodic
 *					 boundary conditions in mind, which particle is when where with which speed and which
 *					 force acting upon it. Every couple of timesteps (varying according to particle positions)
 *					 the verlet list used as a neighbor list for all particles is updated and at regular
 *					 intervals, all position and speed data is written to a file. Several threads are started
 *					 from simulation() which then go on and execute this function until they are aborted.
 *
 *					update_verlet():
 *					 With periodic boundary conditions in mind, this method updates the verlet list.
 *					 Unfortunately several different arrays are needed, which makes the entire algorithm very
 *					 cost intensive with regard to memory. On the other hand this approach also reduces the
 *					 amount of time needed to iterate a single step, as the closest copy of a neighbor does not
 *					 need to be evaluated.
 *
 *					main():
 *					 Takes care of starting functions to read the supplied config files, initialize all data
 *					 and runs the simulation.
 *
 *      AUTHOR: 		Aiko Bernehed
 *
 */


/*-------------------------------------------------------------------------------------------------------*/
#define _GNU_SOURCE

// Standard header files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Own header files
#include "simulation.h"
#include "read_config.h"
#include "structs.h"
#include "hdf5_output.h"



/*-------------------------------------------------------------------------------------------------------*/
int main(int argcount, char** argvektor) {

	/* Read parameters from given file, default is "config.txt", check that input string fits char *infile.
	 * This way several different simulations with different starting parameters are possible. */
	char infile[1024];
	unsigned int length = 1024;

	double* init_positions;

	int sim_number;

	if(argcount == 4) {
		strncpy(infile, argvektor[1], length);
		sim_number = atoi(argvektor[2]);
		init_positions = read_configuration(argvektor[3]);

		// return with a failure if the specified array could not be read
		if (init_positions == NULL) 
			return EXIT_FAILURE;

	} else if (argcount == 3) {
		strncpy(infile, argvektor[1], length);
		sim_number = atoi(argvektor[2]);
		init_positions = NULL;
	} else if (argcount == 2) {
		strncpy(infile, argvektor[1], length);
		sim_number = 0;
		init_positions = NULL;
	} else {
		strncpy(infile,"Config_Files/default", length);		// default
		sim_number = 0;
		init_positions = NULL;
	}

	// get all parameters from file
	struct parameters *param = read_file(infile);

	// check whether the struct was read right, else return failure notice
	if (param == NULL) {
		fprintf(stderr, "Failure in processing configuration file, programme will exit.\n");
		return EXIT_FAILURE;
	}

	// copy the simulation number to the used struct
	param->sim_number = sim_number;

	// Initiate all variables, here a memory problem is most likely the issue at hand
	int success_check = init(param, init_positions);
	if(success_check != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

	// run simulation
	simulation();

	return EXIT_SUCCESS;
}
