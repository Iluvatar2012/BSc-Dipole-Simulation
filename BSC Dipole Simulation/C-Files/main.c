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
 *  	LAST CHANGED: 	Aug 1, 2013
 *      AUTHOR: 		Aiko Bernehed
 *
 */

/*-------------------------------------------------------------------------------------------------------*/
#define _GNU_SOURCE

// Standard header files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>

// Own header files
#include "simulation.h"
#include "fileIO_v1.h"
#include "extendedmath.h"
#include "structs.h"


/*-------------------------------------------------------------------------------------------------------*/
// Basic values of Simulation
static double cutoff;
static double cutoff_squared;
static double delta_t;
static double timestep;
static double box_x;
static int	  max_timesteps;
static int 	  no_writeouts;

// Particle properties
static double kT;
static double Gamma;

// Arrays for mechanical movement
static double* force;
static double* position;

// System properties
static int 	  N;
static double L;
static double Li;
static double shear;
static double D_Brown;
static double D_kT;
static double tau_B;
static double weigh_brown;

// Variables required for Verlet list creation
static int 	   N_Verlet;
static int*    verlet;
static double* verlet_max;
static double* verlet_distance;
static double  verlet_max_1;
static double  verlet_max_2;
static double  d_cutoff_verlet;

// Other miscellaneous system variables, especially for thread handling
static char outfile[1024];
static int* borders;
static int* numbers;
static int thread_number;
static pthread_t* threads;
static pthread_barrier_t barrier_main_one;
static pthread_barrier_t barrier_main_two;
static pthread_barrier_t barrier_internal;

/*-------------------------------------------------------------------------------------------------------*/
int read_struct (char* infile) {

	// get all parameters from file
	struct parameters *param = read_file(infile);

	// check whether the struct was read right, else return failure notice
	if (param == NULL)
		return EXIT_FAILURE;

	// get all the necessary variables from the struct
	strncat(outfile, param->outfile, 1024-strlen(outfile));

	N 		= param->N;
	kT		= param->kT;
	Gamma	= param->Gamma;
	shear	= param->shear;
	tau_B	= param->tau_B;
	D_Brown = param->D_Brown;

	timestep	  	= param->timestep;
	max_timesteps  	= param->max_timesteps;
	thread_number  	= param->thread_number;
	no_writeouts	= param->no_writeouts;

	// free the memory space needed
	free(param);

	// return to caller
	return EXIT_SUCCESS;
}

/*-------------------------------------------------------------------------------------------------------*/
int init(void) {
	// Initiate random number generator (0,1), use system time as seed
    time_t t;
    time(&t);
    srand((unsigned int)t);

	// compute Boxlength from a = sqrt(L²/(PI*N)) = 1, or, more inaccurately, narrow a as side of a box
	// so that: a = sqrt(L²/N) = 1
//	L = sqrt(PI*N);
	L 	= sqrt(N);
	Li 	= 1.0/L;

	// set values of remaining static variables
	delta_t				= tau_B * timestep;
	weigh_brown 		= sqrt(2.0 * D_Brown * delta_t);
	D_kT 				= D_Brown/kT;
	box_x 				= 0;
	cutoff 				= (L/2.0);
	cutoff_squared 		= cutoff*cutoff;
	d_cutoff_verlet 	= 0.16 * cutoff;	// equals 2.9/2.5-1, estimate for best runtime

	// compute the size of the Verlet list, add 20% as safety margin
	N_Verlet = N*PI*cutoff_squared/(L*L);
	N_Verlet *= 1.2;

	// allocate memory for fundamentally important arrays
	position = 			malloc(2*N*sizeof(double));
	force = 			malloc(2*N*sizeof(double));
	verlet = 			malloc(N*N_Verlet*sizeof(int));
	verlet_distance = 	malloc(2*N*sizeof(double));
	threads = 			malloc(thread_number*sizeof(pthread_t));
	borders = 			malloc((thread_number+1)*sizeof(int));
	numbers = 			malloc(thread_number*sizeof(int));
	verlet_max = 		malloc(2*thread_number*sizeof(double));

	if((position == NULL) || (force == NULL) || (verlet == NULL) || (verlet_distance == NULL) || \
			(threads == NULL) || (borders == NULL) || (numbers == NULL) || (verlet_max == NULL)) {
		fprintf(stderr, "Memory space for fundamentally important arrays could not be allocated.\n");
		return EXIT_FAILURE;
	}

	// compute initial particle positions, the particles shall be placed upon an evenly spaced grid
//	double init_x = 0;
//	double init_y = 0;

	for(int i=0; i<(2*N); i+=2) {
//		init_x += L/(sqrt(N));
//		if (init_x >= L) {
//			init_x = (init_x/L-(int)(init_x/L))*L;
//			init_y += L/sqrt(N);
//		}

//		position[i] 	= init_x;
//		position[i+1] 	= init_y;

		// try using random positions for all particles
		position[i]		= (rand()/(double)RAND_MAX)*L;
		position[i+1]	= (rand()/(double)RAND_MAX)*L;
	}

	// initiate Verlet list
	update_verlet();
	return EXIT_SUCCESS;
}


/*-------------------------------------------------------------------------------------------------------*/
static void *iteration (int *no) {
	// Define necessary variables
	double xi, xj, yi, yj;
	double dx, dy, sig_y;
	double r_squared;

	double temp_force;
	double temp;

	double u1, u2, g1, g2;

	int min = borders[*no];
	int max = borders[*no+1];

	int iterate;
	int j;

	while(1) {
		// iterate 2D-forces over all particles in Verlet-List
		for (int i=min; i<max; i++) {
			// get the position of the current particle
			xi = position[2*i];
			yi = position[2*i+1];

			// initiate forces for this round
			force[2*i]	 = 0;
			force[2*i+1] = 0;

			// get the amount of particles we have to iterate
			iterate = verlet[N_Verlet*(i+1)-1];

			// iterate to number given in last entry of i in Verlet-List
			for (int k=0; k<iterate; k+=1) {
				// get the next particle and it's position within the root box
				j  = verlet[N_Verlet*i+k];
				xj = position[2*j];
				yj = position[2*j+1];

				// Get the distance between both particles
				dx = xi - xj;
				dy = yi - yj;

				// alter dx, this accounts for Lees-Edwards conditions
				sig_y 	= dround(dy*Li);
				dx 		+= sig_y * box_x;

				// find images through altering dx and dy
				dx -= dround(dx*Li)*L;
				dy -= sig_y * L;

				// get square of distance and compute force in x and y direction
				r_squared = dx*dx + dy*dy;
				temp_force = 3*Gamma/(r_squared*r_squared*sqrt(r_squared)); // WARNING: this is F/r

				force[2*i] 		+= temp_force*dx; // equals F*x/r = F*cos(phi) (x-component)
				force[2*i+1]	+= temp_force*dy; // equals F*y/r = F*sin(phi) (y-component)
			}
		}

		// wait for all threads to finish iteration of forces, then continue with writing position data
		pthread_barrier_wait(&barrier_internal);

		for (int i=min; i<max; i++) {

			xi = position[2*i];
			yi = position[2*i+1];

			// use Box-Muller-method in order to obtain two independent standard normal values
			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

			// Calculate next positions and speeds for all particles
			position[2*i] 	+= D_kT*delta_t*(force[2*i]+(position[2*i+1]/L*shear)) + weigh_brown * g1;
			position[2*i+1] += D_kT*delta_t*force[2*i+1] + weigh_brown * g2;

			// update list of total displacement for each particle after last verlet-list update
			verlet_distance[2*i] 	+= (xi - position[2*i]);
			verlet_distance[2*i+1] 	+= (yi - position[2*i+1]);

			if ((temp = sqrt(verlet_distance[2*i]*verlet_distance[2*i] + verlet_distance[2*i+1]*verlet_distance[2*i+1])) > verlet_max[2*(*no)]) {
				verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
				verlet_max[2*(*no)] = temp;
			} else if (temp > verlet_max[2*(*no)+1]) {
				verlet_max[2*(*no)+1] = temp;
			}

			// Calculate x and y positions with periodic boundary conditions, save what row the particle travelled from
			position[2*i] 	-= floor(position[2*i]/L)*L;

			temp = (position[2*i+1] + L)/L;
			position[2*i+1] -= floor(position[2*i+1]/L)*L;

			// check whether the particle moved from one box to another, function: ((int)((y+L)/L)-1 gives the signum of the row traveled from (0 if still in the same row)
			position[2*i] -= (((int)(temp))-1)*box_x;

		}
		// Signal to main thread, that all threads have finished their iteration
		pthread_barrier_wait(&barrier_main_one);

		// thread is done with one iteration, it will now wait for a signal from the main thread to continue
		pthread_barrier_wait(&barrier_main_two);
	}

	return NULL;
}


/*-------------------------------------------------------------------------------------------------------*/
void simulation (void) {

	double temp;
	int timesteps;
	int ret_thread;

	// Initiate timestep-counter, write initial composition to file, compute borders for threads and set write-out interval
	timesteps = 0;

	if (write_file(outfile, position, timesteps, no_writeouts, N) == EXIT_FAILURE)
		exit(EXIT_FAILURE);

	timesteps ++;
	no_writeouts = max_timesteps / no_writeouts;

	for(int i=0; i<(thread_number); i++) {
		borders[i] = i*(int)(N/thread_number);
		numbers[i] = i;
	}
	borders[thread_number] = N;

	// Initiate Threads and barrier, catch problems, number of threads is given in config-file
	pthread_barrier_init(&barrier_main_one, NULL, thread_number+1);
	pthread_barrier_init(&barrier_main_two, NULL, thread_number+1);
	pthread_barrier_init(&barrier_internal, NULL, thread_number);

	for (int i=0; i<thread_number; i++) {
		ret_thread 	= pthread_create(&(threads[i]), NULL, (void*)&iteration, &(numbers[i]));

		if (ret_thread != 0) {
			fprintf(stderr, "Thread %d could not be created. \n", i);
		} else {
			fprintf(stderr, "Thread %d initialized...\n", i);
		}
	}

	while (timesteps <= max_timesteps) {
		// Synchronize all threads
//		fprintf(stderr, "Timestep: %d\n", timesteps);
		pthread_barrier_wait(&barrier_main_one);


		// Save maximum displacement of two particles, this determines when to update the Verlet-list
		for (int i=0; i<2*thread_number; i++) {
			if ((temp = verlet_max[i]) > verlet_max_1) {
				verlet_max_2 = verlet_max_1;
				verlet_max_1 = temp;
			} else if (temp > verlet_max_2) {
				verlet_max_2 = temp;
			}
		}

		// adjust orientation of upper and lower box row, check if Verlet list has to be updated
		box_x += shear*L*delta_t;
		box_x  = ((box_x/L)-(int)(box_x/L))*L;

		if ((verlet_max_1+verlet_max_2) > d_cutoff_verlet) {
			update_verlet();
			verlet_max_1 = 0;
			verlet_max_2 = 0;
		}

		// check whether parameters should be written to declared external file
		if ((timesteps%no_writeouts) == 0) {
			write_file(outfile, position, timesteps, no_writeouts, N);
		}

		// increase timesteps and continue waiting threads
		timesteps++;

		pthread_barrier_wait(&barrier_main_two);
	}

	// cancel all other threads, we're done here
	for (int i=0; i<thread_number; i++) {
		pthread_cancel(threads[i]);
	}

	// end of simulation, free all allocated memory space
	free(position);
	free(force);
	free(verlet);
	free(verlet_distance);
	free(threads);
	free(borders);
	free(numbers);

	fprintf(stderr, "exiting simulation\n");
}


/*-------------------------------------------------------------------------------------------------------*/
void update_verlet (void) {
	// define necessary variables: i, j and temporary position variables,
	double xi, yi, xj, yj;
	double dx, dy, sig_y;
	double r_squared;

	// count how many values there are in the list
	int k;

	// iterate over all particles, read positions
	for (int i=0; i<N; i++) {
		xi = position[2*i];
		yi = position[2*i+1];

		// save the amount of neighbors to particle i in k
		k = 0;

		for (int j=0; j<i; j++) {
			// get the position of particle j
			xj = position[2*j];
			yj = position[2*j+1];

			// get the distance between both particles
			dx = xi - xj;
			dy = yi - yj;

			// alter dx, this accounts for Lees-Edwards conditions
			sig_y 	= dround(dy*Li);
			dx 		+= sig_y * box_x;

			// find images through altering dx and dy
			dx -= dround(dx*Li)*L;
			dy -= sig_y * L;

			// find out squared distance and check against the verlet cutoff
			r_squared = dx*dx + dy*dy;

			// squaring is a strict monotonous function, thus we can check with the squares of the values (saves N*N sqrt()-calls)
			if(cutoff_squared >= r_squared) {
				// check whether theere is enough space in the verlet list
				if (k == N_Verlet -1) {
					fprintf(stderr, "Verlet list is too short!!\n");
					fprintf(stderr, "N: %d, N_Verlet: %d, k: %d\n", N, N_Verlet, k);
					exit(EXIT_FAILURE);
				}

				// add neighbor and signums into verlet and sign list
				verlet[N_Verlet*i+k] 	= j;
				k++;
			}
		}

		// ignore entry i=j
		for (int j=(i+1); j<N; j++) {
			// get the position of paritcle j
			xj = position[2*j];
			yj = position[2*j+1];

			// get the distance between both particles
			dx = xi - xj;
			dy = yi - yj;

			// alter dx, this accounts for Lees-Edwards conditions
			sig_y 	= dround(dy*Li);
			dx 		+= sig_y * box_x;

			// find images through altering dx and dy
			dx -= dround(dx*Li)*L;
			dy -= sig_y * L;

			// find out squared distance and check against the verlet cutoff
			r_squared = dx*dx + dy*dy;

			// squaring is a strict monotonous function, thus we can check with the squares of the values (saves N*N sqrt()-calls)
			if(cutoff_squared >= r_squared) {
				// check whether theere is enough space in the verlet list
				if (k == N_Verlet -1) {
					fprintf(stderr, "Verlet list is too short!!\n");
					fprintf(stderr, "N: %d, N_Verlet: %d, k: %d\n", N, N_Verlet, k);
					exit(EXIT_FAILURE);
				}

				// add neighbor and signums into verlet and sign list
				verlet[N_Verlet*i+k] 	= j;
				k++;
			}
		}

		// edit last entry of row, here we store how many neighbors particle i has, and reset distance counting vector
		verlet[N_Verlet*(i+1)-1] 		= k;
		verlet_distance[2*i]	= 0;
		verlet_distance[2*i+1] 	= 0;
	}

	// reset thread specific verlet counter
	for (int i=0; i<thread_number; i++) {
		verlet_max[2*i] 	= 0;
		verlet_max[2*i+1]	= 0;
	}

	fprintf(stderr, "Verlet-list updated...\n");
}


/*-------------------------------------------------------------------------------------------------------*/
int main(int argcount, char** argvektor) {

	/* Read parameters from given file, default is "config.txt", check that input string fits char *infile.
	 * This way several different simulations with different starting parameters are possible. */
	char infile[1024];
	getcwd(infile, sizeof(infile));

	// append "/" since getcwd doesn't include that, but it is needed for good user input
	strncat(infile, "/", 1);
	unsigned int length = sizeof(infile) - strlen(infile);

	if(argcount == 2) {
		strncat(infile, argvektor[1], length);
	} else {
		strncat(infile,"Config_Files/config.txt", length);		// default
	}

	// get current working directory and pass arguments to file reader
	getcwd(outfile, sizeof(outfile));
	int success_check = read_struct(infile);

	// check whether file could be properly read
	if(success_check != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

	// Initiate all variables, here a memory problem is most likely the issue at hand
	success_check = init();
	if(success_check != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

	// run simulation
	simulation();

	return EXIT_SUCCESS;
}
