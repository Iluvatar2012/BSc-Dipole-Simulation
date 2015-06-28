#define _GNU_SOURCE

// define constants of the simulation
#define		N 					300
#define		thread_number		0
#define 	L_y					2
#define 	X_A					0.0

#define 	kT					1.0
#define		tau_B				1.0
#define		D_Brown_A			1.0

// standard header files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <pthread.h>

#include <unistd.h>

// own header files
#include "structs.h"
#include "simulation.h"
#include "read_config.h"
#include "extendedmath.h"
#include "hdf5_output.h"



/*-------------------------------------------------------------------------------------------------------*/
// Basic values of Simulation
static double cutoff;
static double cutoff_squared;

static int 	  tau;
static double delta_t;
static double timestep;
static int	  max_timesteps;

static int 	  write_step;
static int 	  sim_number;

// Particle properties
static double Gamma_A;
static double m;

// Arrays for mechanical movement
static double* force;
static double* position;
static double* displacement;

// System properties
static double L_x;
static double gamma_shear;
static double kappa;
static double v_A;
static double v_B;

static double D_Brown_B;
static double D_rat;
static double weigh_brown_A;
static double weigh_brown_B;

// Variables required for Verlet list creation
static int*    verlet;
static double* verlet_max;
static double* verlet_distance;
static double  verlet_max_1;
static double  verlet_max_2;
static double  d_cutoff_verlet;
static double  force_cutoff;

// Other miscellaneous system variables, especially for thread handling
static char 	outfile[1024];
static int* 	numbers;
static short 	cont;

static pthread_t* 			threads;
static pthread_barrier_t 	barrier_main_one;
static pthread_barrier_t 	barrier_main_two;
static pthread_barrier_t 	barrier_internal;



/*-------------------------------------------------------------------------------------------------------*/
void update_verlet (void) {
	// define necessary variables: i, j and temporary position variables,
	double xi, yi, xj, yj;
	double dx, dy;
	double r_squared;

	// count how many values there are in the list
	int k;

	// iterate over all particles, read positions
	for (int i=0; i<N; i++) {
		xi = position[2*i];
		yi = position[2*i+1];

		// save the amount of neighbors to particle i in k
		k = 0;

		for (int j=0; j<N; j++) {

			// ignore entry where particle i itself is inspected
			if (i == j)
				continue;

			// get the position of particle j
			xj = position[2*j];
			yj = position[2*j+1];

			// get the distance between both particles
			dx = xi - xj;
			dy = yi - yj;

 			// find images through altering dx
			dx -= dround(dx/L_x)*L_x;

			// find out squared distance and check against the verlet cutoff
			r_squared = dx*dx + dy*dy;

			// squaring is a strict monotonous function, thus we can check with the squares of the values (saves N*N sqrt()-calls)
			if(cutoff_squared >= r_squared) {

				// add neighbor and signums into verlet and sign list
				verlet[N*i+k] 	= j;
				k++;
			}
		}

		// edit last entry of row, here we store how many neighbors particle i has, and reset distance counting vector
		verlet[N*(i+1)-1] 		= k;
		verlet_distance[2*i]	= 0;
		verlet_distance[2*i+1] 	= 0;
	}

	// reset thread specific verlet counter
	for (int i=0; i<thread_number; i++) {
		verlet_max[2*i] 	= 0;
		verlet_max[2*i+1]	= 0;
	}
}



/*-------------------------------------------------------------------------------------------------------*/
int init(struct parameters *param, double* init_positions) {

	// copy values from incoming struct
	Gamma_A			= param->Gamma_A;
	m 				= param->m;

	gamma_shear 	= param->gamma_shear;
	D_rat			= param->D_rat;

	timestep 		= param->timestep;
	tau				= param->tau;
	write_step		= param->write_step;
	sim_number		= param->sim_number;

	// check whether particle species orientation is correct
	if (D_rat < 1) {
		fprintf(stderr, "Value D_rat is smaller than 1.0, consider switching particle species properties. D_rat: %lf1.3\n", D_rat);
		return EXIT_FAILURE;
	}

	// copy filename to write to
	strncpy(outfile, param->outfile, 1024);

	// release valuable memory space
	free(param);

	// Initiate random number generator (0,1), use system time as seed
    time_t t;
    time(&t);
	srand((unsigned int)t);

	// compute Boxlength, narrow a as side of a box so that: a = sqrt(Lx*Ly/N) = 1
	L_x 	= N/L_y;

	// set the interaction potential for the walls
	kappa = 10.0;

	// compute diffusion value of particle B, compute box speeds for particles A
	D_Brown_B 			= D_Brown_A * D_rat;
	v_A					= gamma_shear*D_Brown_A/kT;

	// set values of remaining static variables
	delta_t				= tau_B * timestep;
	max_timesteps		= (int)(tau/timestep);

	weigh_brown_A 		= sqrt(2.0 * D_Brown_A * delta_t);
	weigh_brown_B 		= sqrt(2.0 * D_Brown_B * delta_t);

	cutoff 				= (fmax(L_x, L_y)/2.0);
	cutoff_squared 		= cutoff*cutoff;
	d_cutoff_verlet 	= 0.16 * cutoff;	// equals 2.9/2.5-1, estimate for best runtime

	// compute the force at cutoff value, this force will be deducted from the system
	force_cutoff	= 3*Gamma_A/(cutoff_squared*cutoff_squared);

	// check whether an old configuration of data can be used or whether new memory has to be allocated
	if (init_positions != NULL) {
		position = 		init_positions;
	} else {
		position = 		malloc(2*N*sizeof(double));
	}

	// allocate memory for fundamentally important arrays
	force = 			malloc(2*N*sizeof(double));
	displacement = 		malloc(2*N*sizeof(double));
	verlet = 			malloc(N*N*sizeof(int));
	verlet_distance = 	malloc(2*N*sizeof(double));
	threads = 			malloc(thread_number*sizeof(pthread_t));
	numbers = 			malloc(thread_number*sizeof(int));
	verlet_max = 		malloc(2*thread_number*sizeof(double));

	// check if successful
	if((position == NULL) || (displacement == NULL) || (force == NULL) || (verlet == NULL) || (verlet_distance == NULL) || (threads == NULL) || (numbers == NULL) || (verlet_max == NULL)) {
		fprintf(stderr, "Memory space for fundamentally important arrays could not be allocated.\n");
		return EXIT_FAILURE;
	}

	// initiate the numbers array
	for (int i=0; i<thread_number; i++) {
		numbers[i] = i;
	}

	// initiate all displacement values to 0
	for (int i=0; i<N; i++) {
		displacement[2*i]	= 0;
		displacement[2*i+1] = 0;
	}

	// check whether we still need to set the position values, initiate Verlet list
	if (init_positions != NULL) {
		update_verlet();
		return EXIT_SUCCESS;
	}

	// compute initial particle positions, the particles shall be placed randomly with minimum distance criteria
	short done;
	double dx;
	double dy;
	double minDist = 0.2;

	// iterate over all particles
	for (int i=0; i<(2*N); i+=2) {

		// repeat until all particles meet the minimum distance criterium
		do {
			//reset checking variable (1 = TRUE, 0 = FALSE)
			done = 1;

			// try using random positions for all particles
			position[i]		= (rand()/(double)RAND_MAX)*L_x;
			position[i+1]	= (rand()/(double)RAND_MAX)*L_y;

			// check all particles already placed
			for (int j=0; j<i; j+=2) {
				// compute distance of the two particles
				dx = position[j] 	- position[i];
				dy = position[j+1] 	- position[i+1];

				// check whether the minimum distance of particles is met, break loop if not
				if (dx*dx+dy*dy < minDist*minDist || position[i+1] < minDist || position[i+1] > L_y-minDist) {
					done = 0;
					break;
				}
			}
		} while (done == 0);
	}

	// initiate Verlet list
	update_verlet();

	return EXIT_SUCCESS;
}



/*-------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------*/
void simulation (void) {

	char 	job[1024] = "Jobs/";
	char 	buf[32];

	double 	temp;
	double 	perc;

	int 	timesteps;
	int 	ret_thread;

	time_t 	init_time;
	time_t 	current_time;
	time_t 	rem_time;

	char*	time_string;

	// struct for passing all arguments to the hdf5 file
	struct attributes *attr = malloc(sizeof(struct attributes));

	// create a file to which we will write user output
	sprintf(buf, "%d", sim_number);
	strncat(job, buf, 32);

	FILE *file = fopen(job, "w+");
	if (file == NULL) {
		fprintf(stderr, "File \"%s\" could not be opened, program will terminate. \n", job);
		return;
	}

	// Initiate timestep-counter and continuation value
	timesteps = 0;
	cont = 1;

	// fill the attribute struct with values
	attr->Num 			= N;
	attr->Gamma_A		= Gamma_A;
	attr->m 			= m;
	attr->X				= X_A;
	attr->gamma_shear 	= gamma_shear;
	attr->D_rat			= D_rat;
	attr->L_x			= L_x;
	attr->L_y_attr		= L_y;
	attr->max_writeouts = max_timesteps/write_step;
	attr->tau 			= tau;
	attr->write_step	= write_step;
	attr->timestep 		= timestep;


	// initialize file and write first time setup of the system
	if (create_file(outfile, attr) == EXIT_FAILURE) {
		
		// variable for storing the amount of data already written
		int written;
		if (reopen_file(position, displacement, max_timesteps/write_step, &written) == EXIT_FAILURE)
			exit(EXIT_FAILURE);

		// iterate what timestep we are at
		timesteps = write_step * written;
		fprintf(file, "Timestep: %d\n", timesteps);
	}
	else {
		write_data(timesteps, position, displacement);
	}

	// increase timestep counter
	timesteps ++;

	// initiate the threads own numbers
	for (int i=0; i<thread_number; i++) {
		numbers[i] = i;
	}

	// Initiate Threads and barrier, catch problems, number of threads is given in config-file
	pthread_barrier_init(&barrier_main_one, NULL, thread_number+1);
	pthread_barrier_init(&barrier_main_two, NULL, thread_number+1);
	pthread_barrier_init(&barrier_internal, NULL, thread_number);


	// save time at iteration start and output simulation start to user
	time(&init_time);
	time_string = ctime(&init_time);

	fprintf(file, "Starting simulation ID: %d at %s"
					"Parameters N: %d, m: %.2lf, Gamma: %.0lf, Shear: %.0lf, Steps: %d\n\n", sim_number, time_string, N, m, Gamma_A, gamma_shear, max_timesteps);
	fflush(file);

	// iterate over all timesteps
	while (timesteps <= max_timesteps) {
		// Synchronize all threads
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

		// check if verlet list has to be updated
		if ((verlet_max_1+verlet_max_2) > d_cutoff_verlet) {
			update_verlet();
			verlet_max_1 = 0;
			verlet_max_2 = 0;
		}

		// check whether parameters should be written to declared external file
		if ((timesteps%write_step) == 0) {
			write_data(timesteps, position, displacement);

			// compute percentage of program already completed
			perc = (100.*timesteps)/max_timesteps;

			// compute elapsed program time
			time(&current_time);
			rem_time = current_time - init_time;

			// give user a coherent overview
			fprintf(file, "ID: %d, Progress: \t\t[", sim_number);

			for(int i=0; i<floor(perc); i++) {
				fprintf(file, "=");
			}

			if (perc != 100)
				fprintf(file, ">");

			for (int i=floor(perc)+1; i<100; i++) {
				fprintf(file, ".");
			}

			fprintf(file, "] \t%.1lf%%\t\telapsed time: %d s\r", perc, (int)(rem_time));

			fflush(file);
		}

		// increase timesteps and continue waiting threads
		timesteps++;

		pthread_barrier_wait(&barrier_main_two);
	}

	// terminate all threads, wait for every one to finish
	for (int i=0; i<thread_number; i++) {
		pthread_cancel(threads[i]);
	}

	sleep(1);

	// end of simulation, free all allocated memory space
	free(position);
	free(force);
	free(displacement);
	free(verlet);
	free(verlet_distance);
	free(verlet_max);
	free(threads);
	free(numbers);

	// compute runtime, give information to user and leave program
	time(&current_time);
	time_string = ctime(&current_time);

	fprintf(file, "\n\nFinished simulation ID: %d at %s"
					"Parameters N: %d, m: %.2lf, Gamma: %.0lf, Shear: %.0lf, Steps: %d\n"
					"Elapsed time: %d seconds\n\n", sim_number, time_string, N, m, Gamma_A, gamma_shear, max_timesteps, (int)(current_time - init_time));
	fflush(file);
	fclose(file);
}
