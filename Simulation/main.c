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

// define basic variables
#define 	N 					1000
#define		thread_number		8

#define 	kT					1.0
#define		tau_B				1.0
#define		D_Brown_A			1.0



// Standard header files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Own header files
#include "simulation.h"
#include "read_config.h"
#include "structs.h"



/*-------------------------------------------------------------------------------------------------------*/
// Basic values of Simulation
static double timestep;
static int	  max_timesteps;
static int 	  write_step;
static int 	  sim_number;

// Particle properties
static double Gamma_A;
static double m;

// System properties
static double D_rat;
static double shear;

// Other miscellaneous system variables, especially for thread handling
static char 	outfile[1024];



/*-------------------------------------------------------------------------------------------------------*/
int read_struct (char* infile) {

	// get all parameters from file
	struct parameters *param = read_file(infile);

	// check whether the struct was read right, else return failure notice
	if (param == NULL)
		return EXIT_FAILURE;

	// get all the necessary variables from the struct
	strncpy(outfile, param->outfile, 1024);

	Gamma_A		= param->Gamma_A;
	m 			= param->m;
	shear		= param->shear_A;
	D_rat		= param->D_rat;

	timestep	  	= param->timestep;
	max_timesteps  	= param->max_timesteps;
	write_step		= param->write_step;

	// free the memory space needed
	free(param);

	// return to caller
	return EXIT_SUCCESS;
}



/*-------------------------------------------------------------------------------------------------------*/
struct sim_struct *make_sim_struct () {

	// create a new struct which will hold all relevant parameters for the simulation
	struct sim_struct *sim = malloc(sizeof(struct sim_struct));

	// copy simple variables to struct
	sim->struct_N 					= N;
	sim->struct_thread_number		= thread_number;

	sim->struct_kT					= kT;
	sim->struct_tau_B				= tau_B;
	sim->struct_D_Brown_A			= D_Brown_A;

	sim->Gamma_A			= Gamma_A;
	sim->m 					= m;
	sim->shear 				= shear;
	sim->D_rat				= D_rat;

	sim->timestep 			= timestep;
	sim->max_timesteps		= max_timesteps;
	sim->write_step			= write_step;
	sim->sim_number			= sim_number;

	// file to write to
	strncpy(sim->outfile, outfile, 1024);

	// return to caller
	return sim;
}




/*-------------------------------------------------------------------------------------------------------*/
int main(int argcount, char** argvektor) {

	/* Read parameters from given file, default is "config.txt", check that input string fits char *infile.
	 * This way several different simulations with different starting parameters are possible. */
	char infile[1024];
	unsigned int length = 1024;

	if(argcount == 3) {
		strncpy(infile, argvektor[1], length);
		sim_number = atoi(argvektor[2]);
	} else {
		strncpy(infile,"Config_Files/default", length);		// default
		sim_number = 0;
	}

	// pass arguments to file reader
	int success_check = read_struct(infile);

	// check whether file could be properly read
	if(success_check != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

	// create a new struct for the simulation
	struct sim_struct *sim = make_sim_struct();

	// Initiate all variables, here a memory problem is most likely the issue at hand
	success_check = init(sim);
	if(success_check != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

	// run simulation
	simulation();

	return EXIT_SUCCESS;
}
