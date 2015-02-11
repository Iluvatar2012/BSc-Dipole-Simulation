#define _GNU_SOURCE

// define constants of the simulation
#define 	N 					1000
#define		thread_number		16

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
static double L;
static double Li;
static double v_s;
static double v_A;
static double v_B;
static double kappa;

static double D_Brown_B;
static double D_rat;
static double weigh_brown_A;
static double weigh_brown_B;

// Variables required for Verlet list creation
static int 	   N_Verlet;
static int*    verlet;
static double* verlet_max;
static double* verlet_distance;
static double  verlet_max_1;
static double  verlet_max_2;
static double  d_cutoff_verlet;
static double  force_cutoff;

// Other miscellaneous system variables, especially for thread handling
static char 	outfile[1024];
static int* 	borders;
static int* 	numbers;
static short 	cont;

static pthread_t* 			threads;
static pthread_barrier_t 	barrier_main_one;
static pthread_barrier_t 	barrier_main_two;
static pthread_barrier_t 	barrier_internal;



/*-------------------------------------------------------------------------------------------------------*/
int init(struct sim_struct *param, double* init_positions) {

	// copy values from incoming struct
	Gamma_A			= param->Gamma_A;
	m 				= param->m;
	v_s 			= param->shear;
	D_rat			= param->D_rat;

	timestep 		= param->timestep;
	max_timesteps	= param->max_timesteps;
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

	// compute Boxlength, narrow a as side of a box so that: a = sqrt(L²/N_A) = 1
	L 	= sqrt(N/2.);	// N_A = N/2.;
	Li 	= 1.0/L;

	// set the interaction potential for the walls
	kappa = 100;

	// compute diffusion value of particle B, compute box speeds for particles A and B
	// D_Brown_B 			= D_rat*D_Brown_A;
	D_Brown_B 			= D_Brown_A;
	// v_A					= v_s/D_Brown_A*kT;
	v_A					= v_s/D_Brown_A*kT - v_s/D_Brown_B*kT;
	// v_B					= v_s/D_Brown_B*kT;
	v_B					= 0;

	// set values of remaining static variables
	delta_t				= tau_B * timestep;

	weigh_brown_A 		= sqrt(2.0 * D_Brown_A * delta_t);
	weigh_brown_B 		= sqrt(2.0 * D_Brown_B * delta_t);

	cutoff 				= (L/2.0);
	cutoff_squared 		= cutoff*cutoff;
	d_cutoff_verlet 	= 0.16 * cutoff;	// equals 2.9/2.5-1, estimate for best runtime

	// compute the force at cutoff value, this force will be deducted from the system
	force_cutoff	= 3*Gamma_A/(cutoff_squared*cutoff_squared);

	// compute the size of the Verlet list, add 10% as safety margin
	N_Verlet = N*PI*cutoff_squared/(L*L);
	N_Verlet *= 1.1;

	// check whether an old configuration of data can be used or whether new memory has to be allocated
	if (init_positions != NULL) {
		position = 		init_positions;
	} else {
		position = 		malloc(2*N*sizeof(double));
	}

	// allocate memory for fundamentally important arrays
	force = 			malloc(2*N*sizeof(double));
	displacement = 		malloc(2*N*sizeof(double));
	verlet = 			malloc(N*N_Verlet*sizeof(int));
	verlet_distance = 	malloc(2*N*sizeof(double));
	threads = 			malloc(thread_number*sizeof(pthread_t));
	borders = 			malloc((thread_number+1)*sizeof(int));
	numbers = 			malloc(thread_number*sizeof(int));
	verlet_max = 		malloc(2*thread_number*sizeof(double));

	// check if successful
	if((position == NULL) || (displacement == NULL) || (force == NULL) || (verlet == NULL) || (verlet_distance == NULL) || \
		(threads == NULL) || (borders == NULL) || (numbers == NULL) || (verlet_max == NULL)) {
		fprintf(stderr, "Memory space for fundamentally important arrays could not be allocated.\n");
		return EXIT_FAILURE;
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
			position[i]		= (rand()/(double)RAND_MAX)*L;
			position[i+1]	= (rand()/(double)RAND_MAX)*L;

			// check all particles already placed
			for (int j=0; j<i; j+=2) {
				// compute distance of the two particles
				dx = position[j] 	- position[i];
				dy = position[j+1] 	- position[i+1];

				// check whether the minimum distance of particles is met, break loop if not
				if (dx*dx+dy*dy < minDist*minDist) {
					done = 0;
					break;
				}
			}
		} while (done == 0);
	}

	// initiate Verlet list
	update_verlet();

	fprintf(stderr, "D_A: %lf\n", D_Brown_A);
	fprintf(stderr, "D_B: %lf\n", D_Brown_B);
	fprintf(stderr, "v_A: %lf\n", v_A);
	fprintf(stderr, "v_B: %lf\n", v_B);



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

	double m_i_one;
	double m_i_two;

	double weigh_brown_one;
	double weigh_brown_two;

	double D_kT_one;
	double D_kT_two;

	double shear_one;
	double shear_two;

	double m_j;

	if (min%2 == 0){
		// get the relation of all even values, this is basically the interaction relation this particle will have with other particles
		m_i_one = 1.0;
		m_i_two = m;

		// set the appropriate diffusion parameters according to index of particle
		weigh_brown_one = weigh_brown_A;
		weigh_brown_two = weigh_brown_B;

		D_kT_one = D_Brown_A/kT;
		D_kT_two = D_Brown_B/kT;

		// shear_one = v_s / (D_kT_one * L) * delta_t;
		// shear_two = v_s / (D_kT_two * L) * delta_t;
		shear_one = v_A / (L) * delta_t;
		shear_two = v_B / (L) * delta_t;
	}
	else {
		// get the relation of all even values, this is basically the interaction relation this particle will have with other particles
		m_i_one = m;
		m_i_two = 1.0;

		// set the appropriate diffusion parameters according to index of particle
		weigh_brown_one = weigh_brown_B;
		weigh_brown_two = weigh_brown_A;

		D_kT_one = D_Brown_B/kT;
		D_kT_two = D_Brown_A/kT;

		// shear_one = v_s / (D_kT_one * L) * delta_t;
		// shear_two = v_s / (D_kT_two * L) * delta_t;
		shear_one = v_B / (L) * delta_t;
		shear_two = v_A / (L) * delta_t;
	}

	// let the simulation run until the thread is terminated
	while(cont == 1) {
		// iterate 2D-forces over all particles A in Verlet-List
		for (int i=min; i<max; i+=2) {
			// get the position of the current particle
			xi = position[2*i];
			yi = position[2*i+1];

			// initiate forces for this round
			force[2*i]	 = 0;
			force[2*i+1] = 0;

			// add the potential for the walls
			force[2*i+1] += Gamma_A/kappa *(exp(-kappa*yi)*(1/yi + 1/(yi*yi)) + exp(-kappa*(yi-L))*(1/(yi-L) + 1/((yi-L)*(yi-L))));

			// get the amount of particles we have to iterate
			iterate = verlet[N_Verlet*(i+1)-1];

			// iterate to number given in last entry of i in Verlet-List
			for (int k=0; k<iterate; k+=1) {
				// get the next particle and it's position within the root box
				j  = verlet[N_Verlet*i+k];
				xj = position[2*j];
				yj = position[2*j+1];

				// assign the appropriate interaction relation to the particle, depending whether it has even or uneven index
				m_j = (j%2)*m + (j+1)%2;

				// Get the distance between both particles
				dx = xi - xj;
				dy = yi - yj;

 				// find images through altering dx and dy
				dx 		-= dround(dx*Li)*L;

				// get square of distance and compute force in x and y direction
				r_squared = dx*dx + dy*dy;
				temp_force = (m_i_one*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff); // WARNING: this is F/r

				force[2*i] 		+= temp_force*dx; // equals F*x/r = F*cos(phi) (x-component)
				force[2*i+1]	+= temp_force*dy; // equals F*y/r = F*sin(phi) (y-component)
			}
		}

		// iterate 2D-forces over all particles B in Verlet-List
		for (int i=min+1; i<max; i+=2) {
			// get the position of the current particle
			xi = position[2*i];
			yi = position[2*i+1];

			// initiate forces for this round
			force[2*i]	 = 0;
			force[2*i+1] = 0;

			// add the potential for the walls
			force[2*i+1] += Gamma_A/kappa *(exp(-kappa*yi)*(1/yi + 1/(yi*yi)) + exp(-kappa*(yi-L))*(1/(yi-L) + 1/((yi-L)*(yi-L))));

			// get the amount of particles we have to iterate
			iterate = verlet[N_Verlet*(i+1)-1];

			// iterate to number given in last entry of i in Verlet-List
			for (int k=0; k<iterate; k+=1) {
				// get the next particle and it's position within the root box
				j  = verlet[N_Verlet*i+k];
				xj = position[2*j];
				yj = position[2*j+1];

				// assign the appropriate interaction relation to the particle, depending whether it has even or uneven number
				m_j = (j%2)*m + (j+1)%2;

				// Get the distance between both particles
				dx = xi - xj;
				dy = yi - yj;

 				// find images through altering dx and dy
				dx		-= dround(dx*Li)*L;

				// get square of distance and compute force in x and y direction
				r_squared = dx*dx + dy*dy;
				temp_force = (m_i_two*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff); // WARNING: this is F/r

				force[2*i] 		+= temp_force*dx; // equals F*x/r = F*cos(phi) (x-component)
				force[2*i+1]	+= temp_force*dy; // equals F*y/r = F*sin(phi) (y-component)
			}
		}

		// wait for all threads to finish iteration of forces, then continue with writing position data
		pthread_barrier_wait(&barrier_internal);

		// compute new positions for particles A from forces, remember periodic boundary conditions
		for (int i=min; i<max; i+=2) {

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

			// Calculate next positions and total displacement for all particles
			dx = D_kT_one*delta_t*force[2*i]   + weigh_brown_one * g1 + position[2*i+1]*shear_one;
			dy = D_kT_one*delta_t*force[2*i+1] + weigh_brown_one * g2;

			position[2*i]		+= dx;
			position[2*i+1] 	+= dy;

			displacement[2*i]	+= dx;
			displacement[2*i+1] += dy;

			// update list of total displacement for each particle after last verlet-list update
			verlet_distance[2*i] 	+= (xi - position[2*i]);
			verlet_distance[2*i+1] 	+= (yi - position[2*i+1]);

			if ((temp = sqrt(verlet_distance[2*i]*verlet_distance[2*i] + verlet_distance[2*i+1]*verlet_distance[2*i+1])) > verlet_max[2*(*no)]) {
				verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
				verlet_max[2*(*no)] = temp;
			} else if (temp > verlet_max[2*(*no)+1]) {
				verlet_max[2*(*no)+1] = temp;
			}

			// Calculate x positions with periodic boundary conditions, check whether the particle moved from one row to another
			position[2*i] 	-= floor(position[2*i]*Li)*L;
		}

		// compute new positions for particles B from forces, remember periodic boundary conditions
		for (int i=min+1; i<max; i+=2) {

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

			// Calculate next positions and total displacement for all particles
			dx = D_kT_two*delta_t*force[2*i]   + weigh_brown_two * g1 + position[2*i+1]*shear_two;
			dy = D_kT_two*delta_t*force[2*i+1] + weigh_brown_two * g2;

			position[2*i]		+= dx;
			position[2*i+1] 	+= dy;

			displacement[2*i]	+= dx;
			displacement[2*i+1] += dy;

			// update list of total displacement for each particle after last verlet-list update
			verlet_distance[2*i] 	+= (xi - position[2*i]);
			verlet_distance[2*i+1] 	+= (yi - position[2*i+1]);

			if ((temp = sqrt(verlet_distance[2*i]*verlet_distance[2*i] + verlet_distance[2*i+1]*verlet_distance[2*i+1])) > verlet_max[2*(*no)]) {
				verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
				verlet_max[2*(*no)] = temp;
			} else if (temp > verlet_max[2*(*no)+1]) {
				verlet_max[2*(*no)+1] = temp;
			}

			// Calculate x positions with periodic boundary conditions, check whether the particle moved from one row to another
			position[2*i] 	-= floor(position[2*i]*Li)*L;
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

	// initialize file and write first time setup of the system
	if (create_file(outfile, N, max_timesteps/write_step) == EXIT_FAILURE) {
		
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

	// compute borders of all threads
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
			exit(EXIT_FAILURE);
		}
	}

	// save time at iteration start and output simulation start to user
	time(&init_time);
	time_string = ctime(&init_time);

	fprintf(file, "Starting simulation ID: %d at %s"
					"Parameters N: %d, m: %.2lf, Gamma: %.0lf, Shear: %.0lf, Steps: %d\n\n", sim_number, time_string, N, m, Gamma_A, v_s, max_timesteps);
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
	free(borders);
	free(numbers);

	// compute runtime, give information to user and leave program
	time(&current_time);
	time_string = ctime(&current_time);

	fprintf(file, "\n\nFinished simulation ID: %d at %s"
					"Parameters N: %d, m: %.2lf, Gamma: %.0lf, Shear: %.0lf, Steps: %d\n"
					"Elapsed time: %d seconds\n\n", sim_number, time_string, N, m, Gamma_A, v_s, max_timesteps, (int)(current_time - init_time));
	fflush(file);
	fclose(file);
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

 			// find images through altering dx and dy
			dx -= dround(dx*Li)*L;

			// find out squared distance and check against the verlet cutoff
			r_squared = dx*dx + dy*dy;

			// squaring is a strict monotonous function, thus we can check with the squares of the values (saves N*N sqrt()-calls)
			if(cutoff_squared >= r_squared) {
				// check whether there is enough space in the verlet list
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
}