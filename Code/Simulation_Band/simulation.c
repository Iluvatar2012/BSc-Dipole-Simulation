#define _GNU_SOURCE

// define constants of the simulation
#define		N 					300
#define		thread_number		8
#define 	L_y					2
#define 	X_A					0.4

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
static void *iteration_0 (int *no) {
	// Define necessary variables
	double xi, xj, yi, yj;
	double dx, dy;
	double r_squared;

	double temp_force;
	double temp;

	double u1, u2, g1, g2;

	int iterate;
	int j;

	double m_i_A = 1.0;
	double m_i_B = m;

	double D_kT_A = D_Brown_A/kT;
	double D_kT_B = D_Brown_B/kT;

	double shear_A = v_A * delta_t;

	double m_j;

	double dist_bottom;
	double dist_top;

	// cutoff for the walls force, this will lead to a smoother transition
	double dist_cutoff 			= 1.0;
	double prefactor			= 3.0;
	double force_wall_cutoff 	= 1.0*prefactor *(kappa/dist_cutoff + 1.0/(dist_cutoff*dist_cutoff))*exp(-kappa*dist_cutoff);

	while(cont == 1) {
			xi = position[2*0];
			yi = position[2*0+1];

			force[2*0]	 = 0;
			force[2*0+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*0+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*0+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*0+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*0+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(0+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*0+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*0] 		+= temp_force*dx;
			force[2*0+1]	+= temp_force*dy;
		}

			xi = position[2*1];
			yi = position[2*1+1];

			force[2*1]	 = 0;
			force[2*1+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*1+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*1+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*1+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*1+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(1+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*1+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*1] 		+= temp_force*dx;
			force[2*1+1]	+= temp_force*dy;
		}

			xi = position[2*2];
			yi = position[2*2+1];

			force[2*2]	 = 0;
			force[2*2+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*2+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*2+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*2+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*2+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(2+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*2+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*2] 		+= temp_force*dx;
			force[2*2+1]	+= temp_force*dy;
		}

			xi = position[2*3];
			yi = position[2*3+1];

			force[2*3]	 = 0;
			force[2*3+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*3+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*3+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*3+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*3+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(3+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*3+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*3] 		+= temp_force*dx;
			force[2*3+1]	+= temp_force*dy;
		}

			xi = position[2*4];
			yi = position[2*4+1];

			force[2*4]	 = 0;
			force[2*4+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*4+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*4+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*4+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*4+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(4+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*4+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*4] 		+= temp_force*dx;
			force[2*4+1]	+= temp_force*dy;
		}

			xi = position[2*5];
			yi = position[2*5+1];

			force[2*5]	 = 0;
			force[2*5+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*5+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*5+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*5+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*5+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(5+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*5+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*5] 		+= temp_force*dx;
			force[2*5+1]	+= temp_force*dy;
		}

			xi = position[2*6];
			yi = position[2*6+1];

			force[2*6]	 = 0;
			force[2*6+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*6+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*6+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*6+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*6+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(6+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*6+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*6] 		+= temp_force*dx;
			force[2*6+1]	+= temp_force*dy;
		}

			xi = position[2*7];
			yi = position[2*7+1];

			force[2*7]	 = 0;
			force[2*7+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*7+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*7+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*7+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*7+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(7+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*7+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*7] 		+= temp_force*dx;
			force[2*7+1]	+= temp_force*dy;
		}

			xi = position[2*8];
			yi = position[2*8+1];

			force[2*8]	 = 0;
			force[2*8+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*8+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*8+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*8+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*8+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(8+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*8+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*8] 		+= temp_force*dx;
			force[2*8+1]	+= temp_force*dy;
		}

			xi = position[2*9];
			yi = position[2*9+1];

			force[2*9]	 = 0;
			force[2*9+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*9+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*9+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*9+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*9+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(9+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*9+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*9] 		+= temp_force*dx;
			force[2*9+1]	+= temp_force*dy;
		}

			xi = position[2*10];
			yi = position[2*10+1];

			force[2*10]	 = 0;
			force[2*10+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*10+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*10+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*10+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*10+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(10+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*10+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*10] 		+= temp_force*dx;
			force[2*10+1]	+= temp_force*dy;
		}

			xi = position[2*11];
			yi = position[2*11+1];

			force[2*11]	 = 0;
			force[2*11+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*11+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*11+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*11+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*11+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(11+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*11+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*11] 		+= temp_force*dx;
			force[2*11+1]	+= temp_force*dy;
		}

			xi = position[2*12];
			yi = position[2*12+1];

			force[2*12]	 = 0;
			force[2*12+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*12+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*12+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*12+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*12+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(12+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*12+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*12] 		+= temp_force*dx;
			force[2*12+1]	+= temp_force*dy;
		}

			xi = position[2*13];
			yi = position[2*13+1];

			force[2*13]	 = 0;
			force[2*13+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*13+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*13+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*13+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*13+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(13+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*13+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*13] 		+= temp_force*dx;
			force[2*13+1]	+= temp_force*dy;
		}

			xi = position[2*14];
			yi = position[2*14+1];

			force[2*14]	 = 0;
			force[2*14+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*14+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*14+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*14+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*14+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(14+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*14+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*14] 		+= temp_force*dx;
			force[2*14+1]	+= temp_force*dy;
		}

			xi = position[2*15];
			yi = position[2*15+1];

			force[2*15]	 = 0;
			force[2*15+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*15+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*15+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*15+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*15+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(15+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*15+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*15] 		+= temp_force*dx;
			force[2*15+1]	+= temp_force*dy;
		}

			xi = position[2*16];
			yi = position[2*16+1];

			force[2*16]	 = 0;
			force[2*16+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*16+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*16+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*16+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*16+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(16+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*16+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*16] 		+= temp_force*dx;
			force[2*16+1]	+= temp_force*dy;
		}

			xi = position[2*17];
			yi = position[2*17+1];

			force[2*17]	 = 0;
			force[2*17+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*17+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*17+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*17+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*17+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(17+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*17+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*17] 		+= temp_force*dx;
			force[2*17+1]	+= temp_force*dy;
		}

			xi = position[2*18];
			yi = position[2*18+1];

			force[2*18]	 = 0;
			force[2*18+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*18+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*18+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*18+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*18+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(18+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*18+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*18] 		+= temp_force*dx;
			force[2*18+1]	+= temp_force*dy;
		}

			xi = position[2*19];
			yi = position[2*19+1];

			force[2*19]	 = 0;
			force[2*19+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*19+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*19+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*19+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*19+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(19+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*19+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*19] 		+= temp_force*dx;
			force[2*19+1]	+= temp_force*dy;
		}

			xi = position[2*20];
			yi = position[2*20+1];

			force[2*20]	 = 0;
			force[2*20+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*20+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*20+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*20+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*20+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(20+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*20+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*20] 		+= temp_force*dx;
			force[2*20+1]	+= temp_force*dy;
		}

			xi = position[2*21];
			yi = position[2*21+1];

			force[2*21]	 = 0;
			force[2*21+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*21+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*21+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*21+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*21+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(21+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*21+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*21] 		+= temp_force*dx;
			force[2*21+1]	+= temp_force*dy;
		}

			xi = position[2*22];
			yi = position[2*22+1];

			force[2*22]	 = 0;
			force[2*22+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*22+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*22+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*22+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*22+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(22+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*22+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*22] 		+= temp_force*dx;
			force[2*22+1]	+= temp_force*dy;
		}

			xi = position[2*23];
			yi = position[2*23+1];

			force[2*23]	 = 0;
			force[2*23+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*23+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*23+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*23+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*23+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(23+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*23+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*23] 		+= temp_force*dx;
			force[2*23+1]	+= temp_force*dy;
		}

			xi = position[2*24];
			yi = position[2*24+1];

			force[2*24]	 = 0;
			force[2*24+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*24+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*24+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*24+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*24+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(24+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*24+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*24] 		+= temp_force*dx;
			force[2*24+1]	+= temp_force*dy;
		}

			xi = position[2*25];
			yi = position[2*25+1];

			force[2*25]	 = 0;
			force[2*25+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*25+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*25+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*25+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*25+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(25+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*25+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*25] 		+= temp_force*dx;
			force[2*25+1]	+= temp_force*dy;
		}

			xi = position[2*26];
			yi = position[2*26+1];

			force[2*26]	 = 0;
			force[2*26+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*26+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*26+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*26+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*26+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(26+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*26+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*26] 		+= temp_force*dx;
			force[2*26+1]	+= temp_force*dy;
		}

			xi = position[2*27];
			yi = position[2*27+1];

			force[2*27]	 = 0;
			force[2*27+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*27+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*27+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*27+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*27+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(27+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*27+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*27] 		+= temp_force*dx;
			force[2*27+1]	+= temp_force*dy;
		}

			xi = position[2*28];
			yi = position[2*28+1];

			force[2*28]	 = 0;
			force[2*28+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*28+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*28+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*28+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*28+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(28+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*28+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*28] 		+= temp_force*dx;
			force[2*28+1]	+= temp_force*dy;
		}

			xi = position[2*29];
			yi = position[2*29+1];

			force[2*29]	 = 0;
			force[2*29+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*29+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*29+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*29+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*29+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(29+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*29+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*29] 		+= temp_force*dx;
			force[2*29+1]	+= temp_force*dy;
		}

			xi = position[2*30];
			yi = position[2*30+1];

			force[2*30]	 = 0;
			force[2*30+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*30+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*30+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*30+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*30+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(30+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*30+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*30] 		+= temp_force*dx;
			force[2*30+1]	+= temp_force*dy;
		}

			xi = position[2*31];
			yi = position[2*31+1];

			force[2*31]	 = 0;
			force[2*31+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*31+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*31+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*31+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*31+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(31+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*31+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*31] 		+= temp_force*dx;
			force[2*31+1]	+= temp_force*dy;
		}

			xi = position[2*32];
			yi = position[2*32+1];

			force[2*32]	 = 0;
			force[2*32+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*32+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*32+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*32+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*32+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(32+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*32+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*32] 		+= temp_force*dx;
			force[2*32+1]	+= temp_force*dy;
		}

			xi = position[2*33];
			yi = position[2*33+1];

			force[2*33]	 = 0;
			force[2*33+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*33+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*33+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*33+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*33+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(33+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*33+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*33] 		+= temp_force*dx;
			force[2*33+1]	+= temp_force*dy;
		}

			xi = position[2*34];
			yi = position[2*34+1];

			force[2*34]	 = 0;
			force[2*34+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*34+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*34+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*34+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*34+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(34+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*34+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*34] 		+= temp_force*dx;
			force[2*34+1]	+= temp_force*dy;
		}

			xi = position[2*35];
			yi = position[2*35+1];

			force[2*35]	 = 0;
			force[2*35+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*35+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*35+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*35+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*35+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(35+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*35+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*35] 		+= temp_force*dx;
			force[2*35+1]	+= temp_force*dy;
		}

			xi = position[2*36];
			yi = position[2*36+1];

			force[2*36]	 = 0;
			force[2*36+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*36+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*36+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*36+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*36+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(36+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*36+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*36] 		+= temp_force*dx;
			force[2*36+1]	+= temp_force*dy;
		}

			xi = position[2*37];
			yi = position[2*37+1];

			force[2*37]	 = 0;
			force[2*37+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*37+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*37+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*37+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*37+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(37+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*37+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*37] 		+= temp_force*dx;
			force[2*37+1]	+= temp_force*dy;
		}

		pthread_barrier_wait(&barrier_internal);

		xi = position[2*0];
		yi = position[2*0+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*0]   + weigh_brown_A * g1 + (position[2*0+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*0+1] + weigh_brown_A * g2;

		position[2*0]		+= dx;
		position[2*0+1] 	+= dy;

		displacement[2*0]	+= dx;
		displacement[2*0+1] += dy;

		verlet_distance[2*0] 	+= (xi - position[2*0]);
		verlet_distance[2*0+1] 	+= (yi - position[2*0+1]);

		if ((temp = sqrt(verlet_distance[2*0]*verlet_distance[2*0] + verlet_distance[2*0+1]*verlet_distance[2*0+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*0] 	-= floor(position[2*0]/L_x)*L_x;

		xi = position[2*1];
		yi = position[2*1+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*1]   + weigh_brown_A * g1 + (position[2*1+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*1+1] + weigh_brown_A * g2;

		position[2*1]		+= dx;
		position[2*1+1] 	+= dy;

		displacement[2*1]	+= dx;
		displacement[2*1+1] += dy;

		verlet_distance[2*1] 	+= (xi - position[2*1]);
		verlet_distance[2*1+1] 	+= (yi - position[2*1+1]);

		if ((temp = sqrt(verlet_distance[2*1]*verlet_distance[2*1] + verlet_distance[2*1+1]*verlet_distance[2*1+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*1] 	-= floor(position[2*1]/L_x)*L_x;

		xi = position[2*2];
		yi = position[2*2+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*2]   + weigh_brown_A * g1 + (position[2*2+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*2+1] + weigh_brown_A * g2;

		position[2*2]		+= dx;
		position[2*2+1] 	+= dy;

		displacement[2*2]	+= dx;
		displacement[2*2+1] += dy;

		verlet_distance[2*2] 	+= (xi - position[2*2]);
		verlet_distance[2*2+1] 	+= (yi - position[2*2+1]);

		if ((temp = sqrt(verlet_distance[2*2]*verlet_distance[2*2] + verlet_distance[2*2+1]*verlet_distance[2*2+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*2] 	-= floor(position[2*2]/L_x)*L_x;

		xi = position[2*3];
		yi = position[2*3+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*3]   + weigh_brown_A * g1 + (position[2*3+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*3+1] + weigh_brown_A * g2;

		position[2*3]		+= dx;
		position[2*3+1] 	+= dy;

		displacement[2*3]	+= dx;
		displacement[2*3+1] += dy;

		verlet_distance[2*3] 	+= (xi - position[2*3]);
		verlet_distance[2*3+1] 	+= (yi - position[2*3+1]);

		if ((temp = sqrt(verlet_distance[2*3]*verlet_distance[2*3] + verlet_distance[2*3+1]*verlet_distance[2*3+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*3] 	-= floor(position[2*3]/L_x)*L_x;

		xi = position[2*4];
		yi = position[2*4+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*4]   + weigh_brown_A * g1 + (position[2*4+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*4+1] + weigh_brown_A * g2;

		position[2*4]		+= dx;
		position[2*4+1] 	+= dy;

		displacement[2*4]	+= dx;
		displacement[2*4+1] += dy;

		verlet_distance[2*4] 	+= (xi - position[2*4]);
		verlet_distance[2*4+1] 	+= (yi - position[2*4+1]);

		if ((temp = sqrt(verlet_distance[2*4]*verlet_distance[2*4] + verlet_distance[2*4+1]*verlet_distance[2*4+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*4] 	-= floor(position[2*4]/L_x)*L_x;

		xi = position[2*5];
		yi = position[2*5+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*5]   + weigh_brown_A * g1 + (position[2*5+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*5+1] + weigh_brown_A * g2;

		position[2*5]		+= dx;
		position[2*5+1] 	+= dy;

		displacement[2*5]	+= dx;
		displacement[2*5+1] += dy;

		verlet_distance[2*5] 	+= (xi - position[2*5]);
		verlet_distance[2*5+1] 	+= (yi - position[2*5+1]);

		if ((temp = sqrt(verlet_distance[2*5]*verlet_distance[2*5] + verlet_distance[2*5+1]*verlet_distance[2*5+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*5] 	-= floor(position[2*5]/L_x)*L_x;

		xi = position[2*6];
		yi = position[2*6+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*6]   + weigh_brown_A * g1 + (position[2*6+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*6+1] + weigh_brown_A * g2;

		position[2*6]		+= dx;
		position[2*6+1] 	+= dy;

		displacement[2*6]	+= dx;
		displacement[2*6+1] += dy;

		verlet_distance[2*6] 	+= (xi - position[2*6]);
		verlet_distance[2*6+1] 	+= (yi - position[2*6+1]);

		if ((temp = sqrt(verlet_distance[2*6]*verlet_distance[2*6] + verlet_distance[2*6+1]*verlet_distance[2*6+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*6] 	-= floor(position[2*6]/L_x)*L_x;

		xi = position[2*7];
		yi = position[2*7+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*7]   + weigh_brown_A * g1 + (position[2*7+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*7+1] + weigh_brown_A * g2;

		position[2*7]		+= dx;
		position[2*7+1] 	+= dy;

		displacement[2*7]	+= dx;
		displacement[2*7+1] += dy;

		verlet_distance[2*7] 	+= (xi - position[2*7]);
		verlet_distance[2*7+1] 	+= (yi - position[2*7+1]);

		if ((temp = sqrt(verlet_distance[2*7]*verlet_distance[2*7] + verlet_distance[2*7+1]*verlet_distance[2*7+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*7] 	-= floor(position[2*7]/L_x)*L_x;

		xi = position[2*8];
		yi = position[2*8+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*8]   + weigh_brown_A * g1 + (position[2*8+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*8+1] + weigh_brown_A * g2;

		position[2*8]		+= dx;
		position[2*8+1] 	+= dy;

		displacement[2*8]	+= dx;
		displacement[2*8+1] += dy;

		verlet_distance[2*8] 	+= (xi - position[2*8]);
		verlet_distance[2*8+1] 	+= (yi - position[2*8+1]);

		if ((temp = sqrt(verlet_distance[2*8]*verlet_distance[2*8] + verlet_distance[2*8+1]*verlet_distance[2*8+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*8] 	-= floor(position[2*8]/L_x)*L_x;

		xi = position[2*9];
		yi = position[2*9+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*9]   + weigh_brown_A * g1 + (position[2*9+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*9+1] + weigh_brown_A * g2;

		position[2*9]		+= dx;
		position[2*9+1] 	+= dy;

		displacement[2*9]	+= dx;
		displacement[2*9+1] += dy;

		verlet_distance[2*9] 	+= (xi - position[2*9]);
		verlet_distance[2*9+1] 	+= (yi - position[2*9+1]);

		if ((temp = sqrt(verlet_distance[2*9]*verlet_distance[2*9] + verlet_distance[2*9+1]*verlet_distance[2*9+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*9] 	-= floor(position[2*9]/L_x)*L_x;

		xi = position[2*10];
		yi = position[2*10+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*10]   + weigh_brown_A * g1 + (position[2*10+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*10+1] + weigh_brown_A * g2;

		position[2*10]		+= dx;
		position[2*10+1] 	+= dy;

		displacement[2*10]	+= dx;
		displacement[2*10+1] += dy;

		verlet_distance[2*10] 	+= (xi - position[2*10]);
		verlet_distance[2*10+1] 	+= (yi - position[2*10+1]);

		if ((temp = sqrt(verlet_distance[2*10]*verlet_distance[2*10] + verlet_distance[2*10+1]*verlet_distance[2*10+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*10] 	-= floor(position[2*10]/L_x)*L_x;

		xi = position[2*11];
		yi = position[2*11+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*11]   + weigh_brown_A * g1 + (position[2*11+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*11+1] + weigh_brown_A * g2;

		position[2*11]		+= dx;
		position[2*11+1] 	+= dy;

		displacement[2*11]	+= dx;
		displacement[2*11+1] += dy;

		verlet_distance[2*11] 	+= (xi - position[2*11]);
		verlet_distance[2*11+1] 	+= (yi - position[2*11+1]);

		if ((temp = sqrt(verlet_distance[2*11]*verlet_distance[2*11] + verlet_distance[2*11+1]*verlet_distance[2*11+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*11] 	-= floor(position[2*11]/L_x)*L_x;

		xi = position[2*12];
		yi = position[2*12+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*12]   + weigh_brown_A * g1 + (position[2*12+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*12+1] + weigh_brown_A * g2;

		position[2*12]		+= dx;
		position[2*12+1] 	+= dy;

		displacement[2*12]	+= dx;
		displacement[2*12+1] += dy;

		verlet_distance[2*12] 	+= (xi - position[2*12]);
		verlet_distance[2*12+1] 	+= (yi - position[2*12+1]);

		if ((temp = sqrt(verlet_distance[2*12]*verlet_distance[2*12] + verlet_distance[2*12+1]*verlet_distance[2*12+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*12] 	-= floor(position[2*12]/L_x)*L_x;

		xi = position[2*13];
		yi = position[2*13+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*13]   + weigh_brown_A * g1 + (position[2*13+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*13+1] + weigh_brown_A * g2;

		position[2*13]		+= dx;
		position[2*13+1] 	+= dy;

		displacement[2*13]	+= dx;
		displacement[2*13+1] += dy;

		verlet_distance[2*13] 	+= (xi - position[2*13]);
		verlet_distance[2*13+1] 	+= (yi - position[2*13+1]);

		if ((temp = sqrt(verlet_distance[2*13]*verlet_distance[2*13] + verlet_distance[2*13+1]*verlet_distance[2*13+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*13] 	-= floor(position[2*13]/L_x)*L_x;

		xi = position[2*14];
		yi = position[2*14+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*14]   + weigh_brown_A * g1 + (position[2*14+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*14+1] + weigh_brown_A * g2;

		position[2*14]		+= dx;
		position[2*14+1] 	+= dy;

		displacement[2*14]	+= dx;
		displacement[2*14+1] += dy;

		verlet_distance[2*14] 	+= (xi - position[2*14]);
		verlet_distance[2*14+1] 	+= (yi - position[2*14+1]);

		if ((temp = sqrt(verlet_distance[2*14]*verlet_distance[2*14] + verlet_distance[2*14+1]*verlet_distance[2*14+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*14] 	-= floor(position[2*14]/L_x)*L_x;

		xi = position[2*15];
		yi = position[2*15+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*15]   + weigh_brown_A * g1 + (position[2*15+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*15+1] + weigh_brown_A * g2;

		position[2*15]		+= dx;
		position[2*15+1] 	+= dy;

		displacement[2*15]	+= dx;
		displacement[2*15+1] += dy;

		verlet_distance[2*15] 	+= (xi - position[2*15]);
		verlet_distance[2*15+1] 	+= (yi - position[2*15+1]);

		if ((temp = sqrt(verlet_distance[2*15]*verlet_distance[2*15] + verlet_distance[2*15+1]*verlet_distance[2*15+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*15] 	-= floor(position[2*15]/L_x)*L_x;

		xi = position[2*16];
		yi = position[2*16+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*16]   + weigh_brown_A * g1 + (position[2*16+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*16+1] + weigh_brown_A * g2;

		position[2*16]		+= dx;
		position[2*16+1] 	+= dy;

		displacement[2*16]	+= dx;
		displacement[2*16+1] += dy;

		verlet_distance[2*16] 	+= (xi - position[2*16]);
		verlet_distance[2*16+1] 	+= (yi - position[2*16+1]);

		if ((temp = sqrt(verlet_distance[2*16]*verlet_distance[2*16] + verlet_distance[2*16+1]*verlet_distance[2*16+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*16] 	-= floor(position[2*16]/L_x)*L_x;

		xi = position[2*17];
		yi = position[2*17+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*17]   + weigh_brown_A * g1 + (position[2*17+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*17+1] + weigh_brown_A * g2;

		position[2*17]		+= dx;
		position[2*17+1] 	+= dy;

		displacement[2*17]	+= dx;
		displacement[2*17+1] += dy;

		verlet_distance[2*17] 	+= (xi - position[2*17]);
		verlet_distance[2*17+1] 	+= (yi - position[2*17+1]);

		if ((temp = sqrt(verlet_distance[2*17]*verlet_distance[2*17] + verlet_distance[2*17+1]*verlet_distance[2*17+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*17] 	-= floor(position[2*17]/L_x)*L_x;

		xi = position[2*18];
		yi = position[2*18+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*18]   + weigh_brown_A * g1 + (position[2*18+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*18+1] + weigh_brown_A * g2;

		position[2*18]		+= dx;
		position[2*18+1] 	+= dy;

		displacement[2*18]	+= dx;
		displacement[2*18+1] += dy;

		verlet_distance[2*18] 	+= (xi - position[2*18]);
		verlet_distance[2*18+1] 	+= (yi - position[2*18+1]);

		if ((temp = sqrt(verlet_distance[2*18]*verlet_distance[2*18] + verlet_distance[2*18+1]*verlet_distance[2*18+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*18] 	-= floor(position[2*18]/L_x)*L_x;

		xi = position[2*19];
		yi = position[2*19+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*19]   + weigh_brown_A * g1 + (position[2*19+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*19+1] + weigh_brown_A * g2;

		position[2*19]		+= dx;
		position[2*19+1] 	+= dy;

		displacement[2*19]	+= dx;
		displacement[2*19+1] += dy;

		verlet_distance[2*19] 	+= (xi - position[2*19]);
		verlet_distance[2*19+1] 	+= (yi - position[2*19+1]);

		if ((temp = sqrt(verlet_distance[2*19]*verlet_distance[2*19] + verlet_distance[2*19+1]*verlet_distance[2*19+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*19] 	-= floor(position[2*19]/L_x)*L_x;

		xi = position[2*20];
		yi = position[2*20+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*20]   + weigh_brown_A * g1 + (position[2*20+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*20+1] + weigh_brown_A * g2;

		position[2*20]		+= dx;
		position[2*20+1] 	+= dy;

		displacement[2*20]	+= dx;
		displacement[2*20+1] += dy;

		verlet_distance[2*20] 	+= (xi - position[2*20]);
		verlet_distance[2*20+1] 	+= (yi - position[2*20+1]);

		if ((temp = sqrt(verlet_distance[2*20]*verlet_distance[2*20] + verlet_distance[2*20+1]*verlet_distance[2*20+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*20] 	-= floor(position[2*20]/L_x)*L_x;

		xi = position[2*21];
		yi = position[2*21+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*21]   + weigh_brown_A * g1 + (position[2*21+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*21+1] + weigh_brown_A * g2;

		position[2*21]		+= dx;
		position[2*21+1] 	+= dy;

		displacement[2*21]	+= dx;
		displacement[2*21+1] += dy;

		verlet_distance[2*21] 	+= (xi - position[2*21]);
		verlet_distance[2*21+1] 	+= (yi - position[2*21+1]);

		if ((temp = sqrt(verlet_distance[2*21]*verlet_distance[2*21] + verlet_distance[2*21+1]*verlet_distance[2*21+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*21] 	-= floor(position[2*21]/L_x)*L_x;

		xi = position[2*22];
		yi = position[2*22+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*22]   + weigh_brown_A * g1 + (position[2*22+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*22+1] + weigh_brown_A * g2;

		position[2*22]		+= dx;
		position[2*22+1] 	+= dy;

		displacement[2*22]	+= dx;
		displacement[2*22+1] += dy;

		verlet_distance[2*22] 	+= (xi - position[2*22]);
		verlet_distance[2*22+1] 	+= (yi - position[2*22+1]);

		if ((temp = sqrt(verlet_distance[2*22]*verlet_distance[2*22] + verlet_distance[2*22+1]*verlet_distance[2*22+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*22] 	-= floor(position[2*22]/L_x)*L_x;

		xi = position[2*23];
		yi = position[2*23+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*23]   + weigh_brown_A * g1 + (position[2*23+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*23+1] + weigh_brown_A * g2;

		position[2*23]		+= dx;
		position[2*23+1] 	+= dy;

		displacement[2*23]	+= dx;
		displacement[2*23+1] += dy;

		verlet_distance[2*23] 	+= (xi - position[2*23]);
		verlet_distance[2*23+1] 	+= (yi - position[2*23+1]);

		if ((temp = sqrt(verlet_distance[2*23]*verlet_distance[2*23] + verlet_distance[2*23+1]*verlet_distance[2*23+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*23] 	-= floor(position[2*23]/L_x)*L_x;

		xi = position[2*24];
		yi = position[2*24+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*24]   + weigh_brown_A * g1 + (position[2*24+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*24+1] + weigh_brown_A * g2;

		position[2*24]		+= dx;
		position[2*24+1] 	+= dy;

		displacement[2*24]	+= dx;
		displacement[2*24+1] += dy;

		verlet_distance[2*24] 	+= (xi - position[2*24]);
		verlet_distance[2*24+1] 	+= (yi - position[2*24+1]);

		if ((temp = sqrt(verlet_distance[2*24]*verlet_distance[2*24] + verlet_distance[2*24+1]*verlet_distance[2*24+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*24] 	-= floor(position[2*24]/L_x)*L_x;

		xi = position[2*25];
		yi = position[2*25+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*25]   + weigh_brown_A * g1 + (position[2*25+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*25+1] + weigh_brown_A * g2;

		position[2*25]		+= dx;
		position[2*25+1] 	+= dy;

		displacement[2*25]	+= dx;
		displacement[2*25+1] += dy;

		verlet_distance[2*25] 	+= (xi - position[2*25]);
		verlet_distance[2*25+1] 	+= (yi - position[2*25+1]);

		if ((temp = sqrt(verlet_distance[2*25]*verlet_distance[2*25] + verlet_distance[2*25+1]*verlet_distance[2*25+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*25] 	-= floor(position[2*25]/L_x)*L_x;

		xi = position[2*26];
		yi = position[2*26+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*26]   + weigh_brown_A * g1 + (position[2*26+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*26+1] + weigh_brown_A * g2;

		position[2*26]		+= dx;
		position[2*26+1] 	+= dy;

		displacement[2*26]	+= dx;
		displacement[2*26+1] += dy;

		verlet_distance[2*26] 	+= (xi - position[2*26]);
		verlet_distance[2*26+1] 	+= (yi - position[2*26+1]);

		if ((temp = sqrt(verlet_distance[2*26]*verlet_distance[2*26] + verlet_distance[2*26+1]*verlet_distance[2*26+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*26] 	-= floor(position[2*26]/L_x)*L_x;

		xi = position[2*27];
		yi = position[2*27+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*27]   + weigh_brown_A * g1 + (position[2*27+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*27+1] + weigh_brown_A * g2;

		position[2*27]		+= dx;
		position[2*27+1] 	+= dy;

		displacement[2*27]	+= dx;
		displacement[2*27+1] += dy;

		verlet_distance[2*27] 	+= (xi - position[2*27]);
		verlet_distance[2*27+1] 	+= (yi - position[2*27+1]);

		if ((temp = sqrt(verlet_distance[2*27]*verlet_distance[2*27] + verlet_distance[2*27+1]*verlet_distance[2*27+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*27] 	-= floor(position[2*27]/L_x)*L_x;

		xi = position[2*28];
		yi = position[2*28+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*28]   + weigh_brown_A * g1 + (position[2*28+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*28+1] + weigh_brown_A * g2;

		position[2*28]		+= dx;
		position[2*28+1] 	+= dy;

		displacement[2*28]	+= dx;
		displacement[2*28+1] += dy;

		verlet_distance[2*28] 	+= (xi - position[2*28]);
		verlet_distance[2*28+1] 	+= (yi - position[2*28+1]);

		if ((temp = sqrt(verlet_distance[2*28]*verlet_distance[2*28] + verlet_distance[2*28+1]*verlet_distance[2*28+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*28] 	-= floor(position[2*28]/L_x)*L_x;

		xi = position[2*29];
		yi = position[2*29+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*29]   + weigh_brown_A * g1 + (position[2*29+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*29+1] + weigh_brown_A * g2;

		position[2*29]		+= dx;
		position[2*29+1] 	+= dy;

		displacement[2*29]	+= dx;
		displacement[2*29+1] += dy;

		verlet_distance[2*29] 	+= (xi - position[2*29]);
		verlet_distance[2*29+1] 	+= (yi - position[2*29+1]);

		if ((temp = sqrt(verlet_distance[2*29]*verlet_distance[2*29] + verlet_distance[2*29+1]*verlet_distance[2*29+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*29] 	-= floor(position[2*29]/L_x)*L_x;

		xi = position[2*30];
		yi = position[2*30+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*30]   + weigh_brown_A * g1 + (position[2*30+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*30+1] + weigh_brown_A * g2;

		position[2*30]		+= dx;
		position[2*30+1] 	+= dy;

		displacement[2*30]	+= dx;
		displacement[2*30+1] += dy;

		verlet_distance[2*30] 	+= (xi - position[2*30]);
		verlet_distance[2*30+1] 	+= (yi - position[2*30+1]);

		if ((temp = sqrt(verlet_distance[2*30]*verlet_distance[2*30] + verlet_distance[2*30+1]*verlet_distance[2*30+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*30] 	-= floor(position[2*30]/L_x)*L_x;

		xi = position[2*31];
		yi = position[2*31+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*31]   + weigh_brown_A * g1 + (position[2*31+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*31+1] + weigh_brown_A * g2;

		position[2*31]		+= dx;
		position[2*31+1] 	+= dy;

		displacement[2*31]	+= dx;
		displacement[2*31+1] += dy;

		verlet_distance[2*31] 	+= (xi - position[2*31]);
		verlet_distance[2*31+1] 	+= (yi - position[2*31+1]);

		if ((temp = sqrt(verlet_distance[2*31]*verlet_distance[2*31] + verlet_distance[2*31+1]*verlet_distance[2*31+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*31] 	-= floor(position[2*31]/L_x)*L_x;

		xi = position[2*32];
		yi = position[2*32+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*32]   + weigh_brown_A * g1 + (position[2*32+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*32+1] + weigh_brown_A * g2;

		position[2*32]		+= dx;
		position[2*32+1] 	+= dy;

		displacement[2*32]	+= dx;
		displacement[2*32+1] += dy;

		verlet_distance[2*32] 	+= (xi - position[2*32]);
		verlet_distance[2*32+1] 	+= (yi - position[2*32+1]);

		if ((temp = sqrt(verlet_distance[2*32]*verlet_distance[2*32] + verlet_distance[2*32+1]*verlet_distance[2*32+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*32] 	-= floor(position[2*32]/L_x)*L_x;

		xi = position[2*33];
		yi = position[2*33+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*33]   + weigh_brown_A * g1 + (position[2*33+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*33+1] + weigh_brown_A * g2;

		position[2*33]		+= dx;
		position[2*33+1] 	+= dy;

		displacement[2*33]	+= dx;
		displacement[2*33+1] += dy;

		verlet_distance[2*33] 	+= (xi - position[2*33]);
		verlet_distance[2*33+1] 	+= (yi - position[2*33+1]);

		if ((temp = sqrt(verlet_distance[2*33]*verlet_distance[2*33] + verlet_distance[2*33+1]*verlet_distance[2*33+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*33] 	-= floor(position[2*33]/L_x)*L_x;

		xi = position[2*34];
		yi = position[2*34+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*34]   + weigh_brown_A * g1 + (position[2*34+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*34+1] + weigh_brown_A * g2;

		position[2*34]		+= dx;
		position[2*34+1] 	+= dy;

		displacement[2*34]	+= dx;
		displacement[2*34+1] += dy;

		verlet_distance[2*34] 	+= (xi - position[2*34]);
		verlet_distance[2*34+1] 	+= (yi - position[2*34+1]);

		if ((temp = sqrt(verlet_distance[2*34]*verlet_distance[2*34] + verlet_distance[2*34+1]*verlet_distance[2*34+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*34] 	-= floor(position[2*34]/L_x)*L_x;

		xi = position[2*35];
		yi = position[2*35+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*35]   + weigh_brown_A * g1 + (position[2*35+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*35+1] + weigh_brown_A * g2;

		position[2*35]		+= dx;
		position[2*35+1] 	+= dy;

		displacement[2*35]	+= dx;
		displacement[2*35+1] += dy;

		verlet_distance[2*35] 	+= (xi - position[2*35]);
		verlet_distance[2*35+1] 	+= (yi - position[2*35+1]);

		if ((temp = sqrt(verlet_distance[2*35]*verlet_distance[2*35] + verlet_distance[2*35+1]*verlet_distance[2*35+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*35] 	-= floor(position[2*35]/L_x)*L_x;

		xi = position[2*36];
		yi = position[2*36+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*36]   + weigh_brown_A * g1 + (position[2*36+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*36+1] + weigh_brown_A * g2;

		position[2*36]		+= dx;
		position[2*36+1] 	+= dy;

		displacement[2*36]	+= dx;
		displacement[2*36+1] += dy;

		verlet_distance[2*36] 	+= (xi - position[2*36]);
		verlet_distance[2*36+1] 	+= (yi - position[2*36+1]);

		if ((temp = sqrt(verlet_distance[2*36]*verlet_distance[2*36] + verlet_distance[2*36+1]*verlet_distance[2*36+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*36] 	-= floor(position[2*36]/L_x)*L_x;

		xi = position[2*37];
		yi = position[2*37+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*37]   + weigh_brown_A * g1 + (position[2*37+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*37+1] + weigh_brown_A * g2;

		position[2*37]		+= dx;
		position[2*37+1] 	+= dy;

		displacement[2*37]	+= dx;
		displacement[2*37+1] += dy;

		verlet_distance[2*37] 	+= (xi - position[2*37]);
		verlet_distance[2*37+1] 	+= (yi - position[2*37+1]);

		if ((temp = sqrt(verlet_distance[2*37]*verlet_distance[2*37] + verlet_distance[2*37+1]*verlet_distance[2*37+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*37] 	-= floor(position[2*37]/L_x)*L_x;

		// Signal to main thread, that all threads have finished their iteration
		pthread_barrier_wait(&barrier_main_one);

		// thread is done with one iteration, it will now wait for a signal from the main thread to continue
		pthread_barrier_wait(&barrier_main_two);
	}

	return NULL;
}



static void *iteration_1 (int *no) {
	// Define necessary variables
	double xi, xj, yi, yj;
	double dx, dy;
	double r_squared;

	double temp_force;
	double temp;

	double u1, u2, g1, g2;

	int iterate;
	int j;

	double m_i_A = 1.0;
	double m_i_B = m;

	double D_kT_A = D_Brown_A/kT;
	double D_kT_B = D_Brown_B/kT;

	double shear_A = v_A * delta_t;

	double m_j;

	double dist_bottom;
	double dist_top;

	// cutoff for the walls force, this will lead to a smoother transition
	double dist_cutoff 			= 1.0;
	double prefactor			= 3.0;
	double force_wall_cutoff 	= 1.0*prefactor *(kappa/dist_cutoff + 1.0/(dist_cutoff*dist_cutoff))*exp(-kappa*dist_cutoff);

	while(cont == 1) {
			xi = position[2*38];
			yi = position[2*38+1];

			force[2*38]	 = 0;
			force[2*38+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*38+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*38+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*38+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*38+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(38+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*38+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*38] 		+= temp_force*dx;
			force[2*38+1]	+= temp_force*dy;
		}

			xi = position[2*39];
			yi = position[2*39+1];

			force[2*39]	 = 0;
			force[2*39+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*39+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*39+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*39+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*39+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(39+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*39+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*39] 		+= temp_force*dx;
			force[2*39+1]	+= temp_force*dy;
		}

			xi = position[2*40];
			yi = position[2*40+1];

			force[2*40]	 = 0;
			force[2*40+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*40+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*40+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*40+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*40+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(40+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*40+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*40] 		+= temp_force*dx;
			force[2*40+1]	+= temp_force*dy;
		}

			xi = position[2*41];
			yi = position[2*41+1];

			force[2*41]	 = 0;
			force[2*41+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*41+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*41+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*41+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*41+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(41+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*41+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*41] 		+= temp_force*dx;
			force[2*41+1]	+= temp_force*dy;
		}

			xi = position[2*42];
			yi = position[2*42+1];

			force[2*42]	 = 0;
			force[2*42+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*42+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*42+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*42+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*42+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(42+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*42+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*42] 		+= temp_force*dx;
			force[2*42+1]	+= temp_force*dy;
		}

			xi = position[2*43];
			yi = position[2*43+1];

			force[2*43]	 = 0;
			force[2*43+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*43+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*43+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*43+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*43+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(43+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*43+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*43] 		+= temp_force*dx;
			force[2*43+1]	+= temp_force*dy;
		}

			xi = position[2*44];
			yi = position[2*44+1];

			force[2*44]	 = 0;
			force[2*44+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*44+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*44+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*44+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*44+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(44+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*44+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*44] 		+= temp_force*dx;
			force[2*44+1]	+= temp_force*dy;
		}

			xi = position[2*45];
			yi = position[2*45+1];

			force[2*45]	 = 0;
			force[2*45+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*45+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*45+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*45+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*45+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(45+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*45+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*45] 		+= temp_force*dx;
			force[2*45+1]	+= temp_force*dy;
		}

			xi = position[2*46];
			yi = position[2*46+1];

			force[2*46]	 = 0;
			force[2*46+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*46+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*46+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*46+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*46+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(46+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*46+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*46] 		+= temp_force*dx;
			force[2*46+1]	+= temp_force*dy;
		}

			xi = position[2*47];
			yi = position[2*47+1];

			force[2*47]	 = 0;
			force[2*47+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*47+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*47+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*47+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*47+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(47+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*47+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*47] 		+= temp_force*dx;
			force[2*47+1]	+= temp_force*dy;
		}

			xi = position[2*48];
			yi = position[2*48+1];

			force[2*48]	 = 0;
			force[2*48+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*48+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*48+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*48+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*48+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(48+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*48+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*48] 		+= temp_force*dx;
			force[2*48+1]	+= temp_force*dy;
		}

			xi = position[2*49];
			yi = position[2*49+1];

			force[2*49]	 = 0;
			force[2*49+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*49+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*49+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*49+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*49+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(49+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*49+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*49] 		+= temp_force*dx;
			force[2*49+1]	+= temp_force*dy;
		}

			xi = position[2*50];
			yi = position[2*50+1];

			force[2*50]	 = 0;
			force[2*50+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*50+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*50+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*50+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*50+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(50+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*50+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*50] 		+= temp_force*dx;
			force[2*50+1]	+= temp_force*dy;
		}

			xi = position[2*51];
			yi = position[2*51+1];

			force[2*51]	 = 0;
			force[2*51+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*51+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*51+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*51+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*51+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(51+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*51+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*51] 		+= temp_force*dx;
			force[2*51+1]	+= temp_force*dy;
		}

			xi = position[2*52];
			yi = position[2*52+1];

			force[2*52]	 = 0;
			force[2*52+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*52+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*52+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*52+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*52+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(52+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*52+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*52] 		+= temp_force*dx;
			force[2*52+1]	+= temp_force*dy;
		}

			xi = position[2*53];
			yi = position[2*53+1];

			force[2*53]	 = 0;
			force[2*53+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*53+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*53+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*53+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*53+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(53+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*53+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*53] 		+= temp_force*dx;
			force[2*53+1]	+= temp_force*dy;
		}

			xi = position[2*54];
			yi = position[2*54+1];

			force[2*54]	 = 0;
			force[2*54+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*54+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*54+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*54+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*54+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(54+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*54+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*54] 		+= temp_force*dx;
			force[2*54+1]	+= temp_force*dy;
		}

			xi = position[2*55];
			yi = position[2*55+1];

			force[2*55]	 = 0;
			force[2*55+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*55+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*55+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*55+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*55+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(55+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*55+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*55] 		+= temp_force*dx;
			force[2*55+1]	+= temp_force*dy;
		}

			xi = position[2*56];
			yi = position[2*56+1];

			force[2*56]	 = 0;
			force[2*56+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*56+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*56+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*56+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*56+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(56+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*56+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*56] 		+= temp_force*dx;
			force[2*56+1]	+= temp_force*dy;
		}

			xi = position[2*57];
			yi = position[2*57+1];

			force[2*57]	 = 0;
			force[2*57+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*57+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*57+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*57+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*57+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(57+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*57+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*57] 		+= temp_force*dx;
			force[2*57+1]	+= temp_force*dy;
		}

			xi = position[2*58];
			yi = position[2*58+1];

			force[2*58]	 = 0;
			force[2*58+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*58+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*58+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*58+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*58+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(58+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*58+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*58] 		+= temp_force*dx;
			force[2*58+1]	+= temp_force*dy;
		}

			xi = position[2*59];
			yi = position[2*59+1];

			force[2*59]	 = 0;
			force[2*59+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*59+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*59+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*59+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*59+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(59+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*59+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*59] 		+= temp_force*dx;
			force[2*59+1]	+= temp_force*dy;
		}

			xi = position[2*60];
			yi = position[2*60+1];

			force[2*60]	 = 0;
			force[2*60+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*60+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*60+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*60+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*60+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(60+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*60+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*60] 		+= temp_force*dx;
			force[2*60+1]	+= temp_force*dy;
		}

			xi = position[2*61];
			yi = position[2*61+1];

			force[2*61]	 = 0;
			force[2*61+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*61+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*61+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*61+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*61+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(61+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*61+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*61] 		+= temp_force*dx;
			force[2*61+1]	+= temp_force*dy;
		}

			xi = position[2*62];
			yi = position[2*62+1];

			force[2*62]	 = 0;
			force[2*62+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*62+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*62+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*62+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*62+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(62+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*62+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*62] 		+= temp_force*dx;
			force[2*62+1]	+= temp_force*dy;
		}

			xi = position[2*63];
			yi = position[2*63+1];

			force[2*63]	 = 0;
			force[2*63+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*63+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*63+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*63+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*63+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(63+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*63+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*63] 		+= temp_force*dx;
			force[2*63+1]	+= temp_force*dy;
		}

			xi = position[2*64];
			yi = position[2*64+1];

			force[2*64]	 = 0;
			force[2*64+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*64+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*64+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*64+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*64+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(64+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*64+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*64] 		+= temp_force*dx;
			force[2*64+1]	+= temp_force*dy;
		}

			xi = position[2*65];
			yi = position[2*65+1];

			force[2*65]	 = 0;
			force[2*65+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*65+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*65+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*65+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*65+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(65+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*65+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*65] 		+= temp_force*dx;
			force[2*65+1]	+= temp_force*dy;
		}

			xi = position[2*66];
			yi = position[2*66+1];

			force[2*66]	 = 0;
			force[2*66+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*66+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*66+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*66+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*66+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(66+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*66+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*66] 		+= temp_force*dx;
			force[2*66+1]	+= temp_force*dy;
		}

			xi = position[2*67];
			yi = position[2*67+1];

			force[2*67]	 = 0;
			force[2*67+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*67+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*67+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*67+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*67+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(67+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*67+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*67] 		+= temp_force*dx;
			force[2*67+1]	+= temp_force*dy;
		}

			xi = position[2*68];
			yi = position[2*68+1];

			force[2*68]	 = 0;
			force[2*68+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*68+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*68+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*68+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*68+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(68+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*68+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*68] 		+= temp_force*dx;
			force[2*68+1]	+= temp_force*dy;
		}

			xi = position[2*69];
			yi = position[2*69+1];

			force[2*69]	 = 0;
			force[2*69+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*69+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*69+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*69+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*69+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(69+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*69+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*69] 		+= temp_force*dx;
			force[2*69+1]	+= temp_force*dy;
		}

			xi = position[2*70];
			yi = position[2*70+1];

			force[2*70]	 = 0;
			force[2*70+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*70+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*70+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*70+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*70+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(70+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*70+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*70] 		+= temp_force*dx;
			force[2*70+1]	+= temp_force*dy;
		}

			xi = position[2*71];
			yi = position[2*71+1];

			force[2*71]	 = 0;
			force[2*71+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*71+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*71+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*71+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*71+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(71+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*71+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*71] 		+= temp_force*dx;
			force[2*71+1]	+= temp_force*dy;
		}

			xi = position[2*72];
			yi = position[2*72+1];

			force[2*72]	 = 0;
			force[2*72+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*72+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*72+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*72+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*72+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(72+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*72+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*72] 		+= temp_force*dx;
			force[2*72+1]	+= temp_force*dy;
		}

			xi = position[2*73];
			yi = position[2*73+1];

			force[2*73]	 = 0;
			force[2*73+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*73+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*73+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*73+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*73+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(73+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*73+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*73] 		+= temp_force*dx;
			force[2*73+1]	+= temp_force*dy;
		}

			xi = position[2*74];
			yi = position[2*74+1];

			force[2*74]	 = 0;
			force[2*74+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*74+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*74+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*74+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*74+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(74+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*74+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*74] 		+= temp_force*dx;
			force[2*74+1]	+= temp_force*dy;
		}

		pthread_barrier_wait(&barrier_internal);

		xi = position[2*38];
		yi = position[2*38+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*38]   + weigh_brown_A * g1 + (position[2*38+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*38+1] + weigh_brown_A * g2;

		position[2*38]		+= dx;
		position[2*38+1] 	+= dy;

		displacement[2*38]	+= dx;
		displacement[2*38+1] += dy;

		verlet_distance[2*38] 	+= (xi - position[2*38]);
		verlet_distance[2*38+1] 	+= (yi - position[2*38+1]);

		if ((temp = sqrt(verlet_distance[2*38]*verlet_distance[2*38] + verlet_distance[2*38+1]*verlet_distance[2*38+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*38] 	-= floor(position[2*38]/L_x)*L_x;

		xi = position[2*39];
		yi = position[2*39+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*39]   + weigh_brown_A * g1 + (position[2*39+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*39+1] + weigh_brown_A * g2;

		position[2*39]		+= dx;
		position[2*39+1] 	+= dy;

		displacement[2*39]	+= dx;
		displacement[2*39+1] += dy;

		verlet_distance[2*39] 	+= (xi - position[2*39]);
		verlet_distance[2*39+1] 	+= (yi - position[2*39+1]);

		if ((temp = sqrt(verlet_distance[2*39]*verlet_distance[2*39] + verlet_distance[2*39+1]*verlet_distance[2*39+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*39] 	-= floor(position[2*39]/L_x)*L_x;

		xi = position[2*40];
		yi = position[2*40+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*40]   + weigh_brown_A * g1 + (position[2*40+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*40+1] + weigh_brown_A * g2;

		position[2*40]		+= dx;
		position[2*40+1] 	+= dy;

		displacement[2*40]	+= dx;
		displacement[2*40+1] += dy;

		verlet_distance[2*40] 	+= (xi - position[2*40]);
		verlet_distance[2*40+1] 	+= (yi - position[2*40+1]);

		if ((temp = sqrt(verlet_distance[2*40]*verlet_distance[2*40] + verlet_distance[2*40+1]*verlet_distance[2*40+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*40] 	-= floor(position[2*40]/L_x)*L_x;

		xi = position[2*41];
		yi = position[2*41+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*41]   + weigh_brown_A * g1 + (position[2*41+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*41+1] + weigh_brown_A * g2;

		position[2*41]		+= dx;
		position[2*41+1] 	+= dy;

		displacement[2*41]	+= dx;
		displacement[2*41+1] += dy;

		verlet_distance[2*41] 	+= (xi - position[2*41]);
		verlet_distance[2*41+1] 	+= (yi - position[2*41+1]);

		if ((temp = sqrt(verlet_distance[2*41]*verlet_distance[2*41] + verlet_distance[2*41+1]*verlet_distance[2*41+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*41] 	-= floor(position[2*41]/L_x)*L_x;

		xi = position[2*42];
		yi = position[2*42+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*42]   + weigh_brown_A * g1 + (position[2*42+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*42+1] + weigh_brown_A * g2;

		position[2*42]		+= dx;
		position[2*42+1] 	+= dy;

		displacement[2*42]	+= dx;
		displacement[2*42+1] += dy;

		verlet_distance[2*42] 	+= (xi - position[2*42]);
		verlet_distance[2*42+1] 	+= (yi - position[2*42+1]);

		if ((temp = sqrt(verlet_distance[2*42]*verlet_distance[2*42] + verlet_distance[2*42+1]*verlet_distance[2*42+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*42] 	-= floor(position[2*42]/L_x)*L_x;

		xi = position[2*43];
		yi = position[2*43+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*43]   + weigh_brown_A * g1 + (position[2*43+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*43+1] + weigh_brown_A * g2;

		position[2*43]		+= dx;
		position[2*43+1] 	+= dy;

		displacement[2*43]	+= dx;
		displacement[2*43+1] += dy;

		verlet_distance[2*43] 	+= (xi - position[2*43]);
		verlet_distance[2*43+1] 	+= (yi - position[2*43+1]);

		if ((temp = sqrt(verlet_distance[2*43]*verlet_distance[2*43] + verlet_distance[2*43+1]*verlet_distance[2*43+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*43] 	-= floor(position[2*43]/L_x)*L_x;

		xi = position[2*44];
		yi = position[2*44+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*44]   + weigh_brown_A * g1 + (position[2*44+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*44+1] + weigh_brown_A * g2;

		position[2*44]		+= dx;
		position[2*44+1] 	+= dy;

		displacement[2*44]	+= dx;
		displacement[2*44+1] += dy;

		verlet_distance[2*44] 	+= (xi - position[2*44]);
		verlet_distance[2*44+1] 	+= (yi - position[2*44+1]);

		if ((temp = sqrt(verlet_distance[2*44]*verlet_distance[2*44] + verlet_distance[2*44+1]*verlet_distance[2*44+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*44] 	-= floor(position[2*44]/L_x)*L_x;

		xi = position[2*45];
		yi = position[2*45+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*45]   + weigh_brown_A * g1 + (position[2*45+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*45+1] + weigh_brown_A * g2;

		position[2*45]		+= dx;
		position[2*45+1] 	+= dy;

		displacement[2*45]	+= dx;
		displacement[2*45+1] += dy;

		verlet_distance[2*45] 	+= (xi - position[2*45]);
		verlet_distance[2*45+1] 	+= (yi - position[2*45+1]);

		if ((temp = sqrt(verlet_distance[2*45]*verlet_distance[2*45] + verlet_distance[2*45+1]*verlet_distance[2*45+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*45] 	-= floor(position[2*45]/L_x)*L_x;

		xi = position[2*46];
		yi = position[2*46+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*46]   + weigh_brown_A * g1 + (position[2*46+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*46+1] + weigh_brown_A * g2;

		position[2*46]		+= dx;
		position[2*46+1] 	+= dy;

		displacement[2*46]	+= dx;
		displacement[2*46+1] += dy;

		verlet_distance[2*46] 	+= (xi - position[2*46]);
		verlet_distance[2*46+1] 	+= (yi - position[2*46+1]);

		if ((temp = sqrt(verlet_distance[2*46]*verlet_distance[2*46] + verlet_distance[2*46+1]*verlet_distance[2*46+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*46] 	-= floor(position[2*46]/L_x)*L_x;

		xi = position[2*47];
		yi = position[2*47+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*47]   + weigh_brown_A * g1 + (position[2*47+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*47+1] + weigh_brown_A * g2;

		position[2*47]		+= dx;
		position[2*47+1] 	+= dy;

		displacement[2*47]	+= dx;
		displacement[2*47+1] += dy;

		verlet_distance[2*47] 	+= (xi - position[2*47]);
		verlet_distance[2*47+1] 	+= (yi - position[2*47+1]);

		if ((temp = sqrt(verlet_distance[2*47]*verlet_distance[2*47] + verlet_distance[2*47+1]*verlet_distance[2*47+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*47] 	-= floor(position[2*47]/L_x)*L_x;

		xi = position[2*48];
		yi = position[2*48+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*48]   + weigh_brown_A * g1 + (position[2*48+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*48+1] + weigh_brown_A * g2;

		position[2*48]		+= dx;
		position[2*48+1] 	+= dy;

		displacement[2*48]	+= dx;
		displacement[2*48+1] += dy;

		verlet_distance[2*48] 	+= (xi - position[2*48]);
		verlet_distance[2*48+1] 	+= (yi - position[2*48+1]);

		if ((temp = sqrt(verlet_distance[2*48]*verlet_distance[2*48] + verlet_distance[2*48+1]*verlet_distance[2*48+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*48] 	-= floor(position[2*48]/L_x)*L_x;

		xi = position[2*49];
		yi = position[2*49+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*49]   + weigh_brown_A * g1 + (position[2*49+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*49+1] + weigh_brown_A * g2;

		position[2*49]		+= dx;
		position[2*49+1] 	+= dy;

		displacement[2*49]	+= dx;
		displacement[2*49+1] += dy;

		verlet_distance[2*49] 	+= (xi - position[2*49]);
		verlet_distance[2*49+1] 	+= (yi - position[2*49+1]);

		if ((temp = sqrt(verlet_distance[2*49]*verlet_distance[2*49] + verlet_distance[2*49+1]*verlet_distance[2*49+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*49] 	-= floor(position[2*49]/L_x)*L_x;

		xi = position[2*50];
		yi = position[2*50+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*50]   + weigh_brown_A * g1 + (position[2*50+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*50+1] + weigh_brown_A * g2;

		position[2*50]		+= dx;
		position[2*50+1] 	+= dy;

		displacement[2*50]	+= dx;
		displacement[2*50+1] += dy;

		verlet_distance[2*50] 	+= (xi - position[2*50]);
		verlet_distance[2*50+1] 	+= (yi - position[2*50+1]);

		if ((temp = sqrt(verlet_distance[2*50]*verlet_distance[2*50] + verlet_distance[2*50+1]*verlet_distance[2*50+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*50] 	-= floor(position[2*50]/L_x)*L_x;

		xi = position[2*51];
		yi = position[2*51+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*51]   + weigh_brown_A * g1 + (position[2*51+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*51+1] + weigh_brown_A * g2;

		position[2*51]		+= dx;
		position[2*51+1] 	+= dy;

		displacement[2*51]	+= dx;
		displacement[2*51+1] += dy;

		verlet_distance[2*51] 	+= (xi - position[2*51]);
		verlet_distance[2*51+1] 	+= (yi - position[2*51+1]);

		if ((temp = sqrt(verlet_distance[2*51]*verlet_distance[2*51] + verlet_distance[2*51+1]*verlet_distance[2*51+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*51] 	-= floor(position[2*51]/L_x)*L_x;

		xi = position[2*52];
		yi = position[2*52+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*52]   + weigh_brown_A * g1 + (position[2*52+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*52+1] + weigh_brown_A * g2;

		position[2*52]		+= dx;
		position[2*52+1] 	+= dy;

		displacement[2*52]	+= dx;
		displacement[2*52+1] += dy;

		verlet_distance[2*52] 	+= (xi - position[2*52]);
		verlet_distance[2*52+1] 	+= (yi - position[2*52+1]);

		if ((temp = sqrt(verlet_distance[2*52]*verlet_distance[2*52] + verlet_distance[2*52+1]*verlet_distance[2*52+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*52] 	-= floor(position[2*52]/L_x)*L_x;

		xi = position[2*53];
		yi = position[2*53+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*53]   + weigh_brown_A * g1 + (position[2*53+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*53+1] + weigh_brown_A * g2;

		position[2*53]		+= dx;
		position[2*53+1] 	+= dy;

		displacement[2*53]	+= dx;
		displacement[2*53+1] += dy;

		verlet_distance[2*53] 	+= (xi - position[2*53]);
		verlet_distance[2*53+1] 	+= (yi - position[2*53+1]);

		if ((temp = sqrt(verlet_distance[2*53]*verlet_distance[2*53] + verlet_distance[2*53+1]*verlet_distance[2*53+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*53] 	-= floor(position[2*53]/L_x)*L_x;

		xi = position[2*54];
		yi = position[2*54+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*54]   + weigh_brown_A * g1 + (position[2*54+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*54+1] + weigh_brown_A * g2;

		position[2*54]		+= dx;
		position[2*54+1] 	+= dy;

		displacement[2*54]	+= dx;
		displacement[2*54+1] += dy;

		verlet_distance[2*54] 	+= (xi - position[2*54]);
		verlet_distance[2*54+1] 	+= (yi - position[2*54+1]);

		if ((temp = sqrt(verlet_distance[2*54]*verlet_distance[2*54] + verlet_distance[2*54+1]*verlet_distance[2*54+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*54] 	-= floor(position[2*54]/L_x)*L_x;

		xi = position[2*55];
		yi = position[2*55+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*55]   + weigh_brown_A * g1 + (position[2*55+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*55+1] + weigh_brown_A * g2;

		position[2*55]		+= dx;
		position[2*55+1] 	+= dy;

		displacement[2*55]	+= dx;
		displacement[2*55+1] += dy;

		verlet_distance[2*55] 	+= (xi - position[2*55]);
		verlet_distance[2*55+1] 	+= (yi - position[2*55+1]);

		if ((temp = sqrt(verlet_distance[2*55]*verlet_distance[2*55] + verlet_distance[2*55+1]*verlet_distance[2*55+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*55] 	-= floor(position[2*55]/L_x)*L_x;

		xi = position[2*56];
		yi = position[2*56+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*56]   + weigh_brown_A * g1 + (position[2*56+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*56+1] + weigh_brown_A * g2;

		position[2*56]		+= dx;
		position[2*56+1] 	+= dy;

		displacement[2*56]	+= dx;
		displacement[2*56+1] += dy;

		verlet_distance[2*56] 	+= (xi - position[2*56]);
		verlet_distance[2*56+1] 	+= (yi - position[2*56+1]);

		if ((temp = sqrt(verlet_distance[2*56]*verlet_distance[2*56] + verlet_distance[2*56+1]*verlet_distance[2*56+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*56] 	-= floor(position[2*56]/L_x)*L_x;

		xi = position[2*57];
		yi = position[2*57+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*57]   + weigh_brown_A * g1 + (position[2*57+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*57+1] + weigh_brown_A * g2;

		position[2*57]		+= dx;
		position[2*57+1] 	+= dy;

		displacement[2*57]	+= dx;
		displacement[2*57+1] += dy;

		verlet_distance[2*57] 	+= (xi - position[2*57]);
		verlet_distance[2*57+1] 	+= (yi - position[2*57+1]);

		if ((temp = sqrt(verlet_distance[2*57]*verlet_distance[2*57] + verlet_distance[2*57+1]*verlet_distance[2*57+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*57] 	-= floor(position[2*57]/L_x)*L_x;

		xi = position[2*58];
		yi = position[2*58+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*58]   + weigh_brown_A * g1 + (position[2*58+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*58+1] + weigh_brown_A * g2;

		position[2*58]		+= dx;
		position[2*58+1] 	+= dy;

		displacement[2*58]	+= dx;
		displacement[2*58+1] += dy;

		verlet_distance[2*58] 	+= (xi - position[2*58]);
		verlet_distance[2*58+1] 	+= (yi - position[2*58+1]);

		if ((temp = sqrt(verlet_distance[2*58]*verlet_distance[2*58] + verlet_distance[2*58+1]*verlet_distance[2*58+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*58] 	-= floor(position[2*58]/L_x)*L_x;

		xi = position[2*59];
		yi = position[2*59+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*59]   + weigh_brown_A * g1 + (position[2*59+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*59+1] + weigh_brown_A * g2;

		position[2*59]		+= dx;
		position[2*59+1] 	+= dy;

		displacement[2*59]	+= dx;
		displacement[2*59+1] += dy;

		verlet_distance[2*59] 	+= (xi - position[2*59]);
		verlet_distance[2*59+1] 	+= (yi - position[2*59+1]);

		if ((temp = sqrt(verlet_distance[2*59]*verlet_distance[2*59] + verlet_distance[2*59+1]*verlet_distance[2*59+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*59] 	-= floor(position[2*59]/L_x)*L_x;

		xi = position[2*60];
		yi = position[2*60+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*60]   + weigh_brown_A * g1 + (position[2*60+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*60+1] + weigh_brown_A * g2;

		position[2*60]		+= dx;
		position[2*60+1] 	+= dy;

		displacement[2*60]	+= dx;
		displacement[2*60+1] += dy;

		verlet_distance[2*60] 	+= (xi - position[2*60]);
		verlet_distance[2*60+1] 	+= (yi - position[2*60+1]);

		if ((temp = sqrt(verlet_distance[2*60]*verlet_distance[2*60] + verlet_distance[2*60+1]*verlet_distance[2*60+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*60] 	-= floor(position[2*60]/L_x)*L_x;

		xi = position[2*61];
		yi = position[2*61+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*61]   + weigh_brown_A * g1 + (position[2*61+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*61+1] + weigh_brown_A * g2;

		position[2*61]		+= dx;
		position[2*61+1] 	+= dy;

		displacement[2*61]	+= dx;
		displacement[2*61+1] += dy;

		verlet_distance[2*61] 	+= (xi - position[2*61]);
		verlet_distance[2*61+1] 	+= (yi - position[2*61+1]);

		if ((temp = sqrt(verlet_distance[2*61]*verlet_distance[2*61] + verlet_distance[2*61+1]*verlet_distance[2*61+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*61] 	-= floor(position[2*61]/L_x)*L_x;

		xi = position[2*62];
		yi = position[2*62+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*62]   + weigh_brown_A * g1 + (position[2*62+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*62+1] + weigh_brown_A * g2;

		position[2*62]		+= dx;
		position[2*62+1] 	+= dy;

		displacement[2*62]	+= dx;
		displacement[2*62+1] += dy;

		verlet_distance[2*62] 	+= (xi - position[2*62]);
		verlet_distance[2*62+1] 	+= (yi - position[2*62+1]);

		if ((temp = sqrt(verlet_distance[2*62]*verlet_distance[2*62] + verlet_distance[2*62+1]*verlet_distance[2*62+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*62] 	-= floor(position[2*62]/L_x)*L_x;

		xi = position[2*63];
		yi = position[2*63+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*63]   + weigh_brown_A * g1 + (position[2*63+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*63+1] + weigh_brown_A * g2;

		position[2*63]		+= dx;
		position[2*63+1] 	+= dy;

		displacement[2*63]	+= dx;
		displacement[2*63+1] += dy;

		verlet_distance[2*63] 	+= (xi - position[2*63]);
		verlet_distance[2*63+1] 	+= (yi - position[2*63+1]);

		if ((temp = sqrt(verlet_distance[2*63]*verlet_distance[2*63] + verlet_distance[2*63+1]*verlet_distance[2*63+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*63] 	-= floor(position[2*63]/L_x)*L_x;

		xi = position[2*64];
		yi = position[2*64+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*64]   + weigh_brown_A * g1 + (position[2*64+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*64+1] + weigh_brown_A * g2;

		position[2*64]		+= dx;
		position[2*64+1] 	+= dy;

		displacement[2*64]	+= dx;
		displacement[2*64+1] += dy;

		verlet_distance[2*64] 	+= (xi - position[2*64]);
		verlet_distance[2*64+1] 	+= (yi - position[2*64+1]);

		if ((temp = sqrt(verlet_distance[2*64]*verlet_distance[2*64] + verlet_distance[2*64+1]*verlet_distance[2*64+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*64] 	-= floor(position[2*64]/L_x)*L_x;

		xi = position[2*65];
		yi = position[2*65+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*65]   + weigh_brown_A * g1 + (position[2*65+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*65+1] + weigh_brown_A * g2;

		position[2*65]		+= dx;
		position[2*65+1] 	+= dy;

		displacement[2*65]	+= dx;
		displacement[2*65+1] += dy;

		verlet_distance[2*65] 	+= (xi - position[2*65]);
		verlet_distance[2*65+1] 	+= (yi - position[2*65+1]);

		if ((temp = sqrt(verlet_distance[2*65]*verlet_distance[2*65] + verlet_distance[2*65+1]*verlet_distance[2*65+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*65] 	-= floor(position[2*65]/L_x)*L_x;

		xi = position[2*66];
		yi = position[2*66+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*66]   + weigh_brown_A * g1 + (position[2*66+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*66+1] + weigh_brown_A * g2;

		position[2*66]		+= dx;
		position[2*66+1] 	+= dy;

		displacement[2*66]	+= dx;
		displacement[2*66+1] += dy;

		verlet_distance[2*66] 	+= (xi - position[2*66]);
		verlet_distance[2*66+1] 	+= (yi - position[2*66+1]);

		if ((temp = sqrt(verlet_distance[2*66]*verlet_distance[2*66] + verlet_distance[2*66+1]*verlet_distance[2*66+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*66] 	-= floor(position[2*66]/L_x)*L_x;

		xi = position[2*67];
		yi = position[2*67+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*67]   + weigh_brown_A * g1 + (position[2*67+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*67+1] + weigh_brown_A * g2;

		position[2*67]		+= dx;
		position[2*67+1] 	+= dy;

		displacement[2*67]	+= dx;
		displacement[2*67+1] += dy;

		verlet_distance[2*67] 	+= (xi - position[2*67]);
		verlet_distance[2*67+1] 	+= (yi - position[2*67+1]);

		if ((temp = sqrt(verlet_distance[2*67]*verlet_distance[2*67] + verlet_distance[2*67+1]*verlet_distance[2*67+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*67] 	-= floor(position[2*67]/L_x)*L_x;

		xi = position[2*68];
		yi = position[2*68+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*68]   + weigh_brown_A * g1 + (position[2*68+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*68+1] + weigh_brown_A * g2;

		position[2*68]		+= dx;
		position[2*68+1] 	+= dy;

		displacement[2*68]	+= dx;
		displacement[2*68+1] += dy;

		verlet_distance[2*68] 	+= (xi - position[2*68]);
		verlet_distance[2*68+1] 	+= (yi - position[2*68+1]);

		if ((temp = sqrt(verlet_distance[2*68]*verlet_distance[2*68] + verlet_distance[2*68+1]*verlet_distance[2*68+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*68] 	-= floor(position[2*68]/L_x)*L_x;

		xi = position[2*69];
		yi = position[2*69+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*69]   + weigh_brown_A * g1 + (position[2*69+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*69+1] + weigh_brown_A * g2;

		position[2*69]		+= dx;
		position[2*69+1] 	+= dy;

		displacement[2*69]	+= dx;
		displacement[2*69+1] += dy;

		verlet_distance[2*69] 	+= (xi - position[2*69]);
		verlet_distance[2*69+1] 	+= (yi - position[2*69+1]);

		if ((temp = sqrt(verlet_distance[2*69]*verlet_distance[2*69] + verlet_distance[2*69+1]*verlet_distance[2*69+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*69] 	-= floor(position[2*69]/L_x)*L_x;

		xi = position[2*70];
		yi = position[2*70+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*70]   + weigh_brown_A * g1 + (position[2*70+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*70+1] + weigh_brown_A * g2;

		position[2*70]		+= dx;
		position[2*70+1] 	+= dy;

		displacement[2*70]	+= dx;
		displacement[2*70+1] += dy;

		verlet_distance[2*70] 	+= (xi - position[2*70]);
		verlet_distance[2*70+1] 	+= (yi - position[2*70+1]);

		if ((temp = sqrt(verlet_distance[2*70]*verlet_distance[2*70] + verlet_distance[2*70+1]*verlet_distance[2*70+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*70] 	-= floor(position[2*70]/L_x)*L_x;

		xi = position[2*71];
		yi = position[2*71+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*71]   + weigh_brown_A * g1 + (position[2*71+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*71+1] + weigh_brown_A * g2;

		position[2*71]		+= dx;
		position[2*71+1] 	+= dy;

		displacement[2*71]	+= dx;
		displacement[2*71+1] += dy;

		verlet_distance[2*71] 	+= (xi - position[2*71]);
		verlet_distance[2*71+1] 	+= (yi - position[2*71+1]);

		if ((temp = sqrt(verlet_distance[2*71]*verlet_distance[2*71] + verlet_distance[2*71+1]*verlet_distance[2*71+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*71] 	-= floor(position[2*71]/L_x)*L_x;

		xi = position[2*72];
		yi = position[2*72+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*72]   + weigh_brown_A * g1 + (position[2*72+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*72+1] + weigh_brown_A * g2;

		position[2*72]		+= dx;
		position[2*72+1] 	+= dy;

		displacement[2*72]	+= dx;
		displacement[2*72+1] += dy;

		verlet_distance[2*72] 	+= (xi - position[2*72]);
		verlet_distance[2*72+1] 	+= (yi - position[2*72+1]);

		if ((temp = sqrt(verlet_distance[2*72]*verlet_distance[2*72] + verlet_distance[2*72+1]*verlet_distance[2*72+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*72] 	-= floor(position[2*72]/L_x)*L_x;

		xi = position[2*73];
		yi = position[2*73+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*73]   + weigh_brown_A * g1 + (position[2*73+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*73+1] + weigh_brown_A * g2;

		position[2*73]		+= dx;
		position[2*73+1] 	+= dy;

		displacement[2*73]	+= dx;
		displacement[2*73+1] += dy;

		verlet_distance[2*73] 	+= (xi - position[2*73]);
		verlet_distance[2*73+1] 	+= (yi - position[2*73+1]);

		if ((temp = sqrt(verlet_distance[2*73]*verlet_distance[2*73] + verlet_distance[2*73+1]*verlet_distance[2*73+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*73] 	-= floor(position[2*73]/L_x)*L_x;

		xi = position[2*74];
		yi = position[2*74+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*74]   + weigh_brown_A * g1 + (position[2*74+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*74+1] + weigh_brown_A * g2;

		position[2*74]		+= dx;
		position[2*74+1] 	+= dy;

		displacement[2*74]	+= dx;
		displacement[2*74+1] += dy;

		verlet_distance[2*74] 	+= (xi - position[2*74]);
		verlet_distance[2*74+1] 	+= (yi - position[2*74+1]);

		if ((temp = sqrt(verlet_distance[2*74]*verlet_distance[2*74] + verlet_distance[2*74+1]*verlet_distance[2*74+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*74] 	-= floor(position[2*74]/L_x)*L_x;

		// Signal to main thread, that all threads have finished their iteration
		pthread_barrier_wait(&barrier_main_one);

		// thread is done with one iteration, it will now wait for a signal from the main thread to continue
		pthread_barrier_wait(&barrier_main_two);
	}

	return NULL;
}



static void *iteration_2 (int *no) {
	// Define necessary variables
	double xi, xj, yi, yj;
	double dx, dy;
	double r_squared;

	double temp_force;
	double temp;

	double u1, u2, g1, g2;

	int iterate;
	int j;

	double m_i_A = 1.0;
	double m_i_B = m;

	double D_kT_A = D_Brown_A/kT;
	double D_kT_B = D_Brown_B/kT;

	double shear_A = v_A * delta_t;

	double m_j;

	double dist_bottom;
	double dist_top;

	// cutoff for the walls force, this will lead to a smoother transition
	double dist_cutoff 			= 1.0;
	double prefactor			= 3.0;
	double force_wall_cutoff 	= 1.0*prefactor *(kappa/dist_cutoff + 1.0/(dist_cutoff*dist_cutoff))*exp(-kappa*dist_cutoff);

	while(cont == 1) {
			xi = position[2*75];
			yi = position[2*75+1];

			force[2*75]	 = 0;
			force[2*75+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*75+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*75+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*75+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*75+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(75+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*75+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*75] 		+= temp_force*dx;
			force[2*75+1]	+= temp_force*dy;
		}

			xi = position[2*76];
			yi = position[2*76+1];

			force[2*76]	 = 0;
			force[2*76+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*76+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*76+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*76+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*76+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(76+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*76+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*76] 		+= temp_force*dx;
			force[2*76+1]	+= temp_force*dy;
		}

			xi = position[2*77];
			yi = position[2*77+1];

			force[2*77]	 = 0;
			force[2*77+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*77+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*77+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*77+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*77+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(77+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*77+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*77] 		+= temp_force*dx;
			force[2*77+1]	+= temp_force*dy;
		}

			xi = position[2*78];
			yi = position[2*78+1];

			force[2*78]	 = 0;
			force[2*78+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*78+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*78+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*78+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*78+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(78+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*78+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*78] 		+= temp_force*dx;
			force[2*78+1]	+= temp_force*dy;
		}

			xi = position[2*79];
			yi = position[2*79+1];

			force[2*79]	 = 0;
			force[2*79+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*79+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*79+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*79+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*79+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(79+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*79+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*79] 		+= temp_force*dx;
			force[2*79+1]	+= temp_force*dy;
		}

			xi = position[2*80];
			yi = position[2*80+1];

			force[2*80]	 = 0;
			force[2*80+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*80+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*80+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*80+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*80+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(80+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*80+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*80] 		+= temp_force*dx;
			force[2*80+1]	+= temp_force*dy;
		}

			xi = position[2*81];
			yi = position[2*81+1];

			force[2*81]	 = 0;
			force[2*81+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*81+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*81+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*81+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*81+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(81+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*81+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*81] 		+= temp_force*dx;
			force[2*81+1]	+= temp_force*dy;
		}

			xi = position[2*82];
			yi = position[2*82+1];

			force[2*82]	 = 0;
			force[2*82+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*82+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*82+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*82+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*82+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(82+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*82+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*82] 		+= temp_force*dx;
			force[2*82+1]	+= temp_force*dy;
		}

			xi = position[2*83];
			yi = position[2*83+1];

			force[2*83]	 = 0;
			force[2*83+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*83+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*83+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*83+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*83+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(83+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*83+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*83] 		+= temp_force*dx;
			force[2*83+1]	+= temp_force*dy;
		}

			xi = position[2*84];
			yi = position[2*84+1];

			force[2*84]	 = 0;
			force[2*84+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*84+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*84+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*84+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*84+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(84+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*84+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*84] 		+= temp_force*dx;
			force[2*84+1]	+= temp_force*dy;
		}

			xi = position[2*85];
			yi = position[2*85+1];

			force[2*85]	 = 0;
			force[2*85+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*85+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*85+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*85+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*85+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(85+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*85+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*85] 		+= temp_force*dx;
			force[2*85+1]	+= temp_force*dy;
		}

			xi = position[2*86];
			yi = position[2*86+1];

			force[2*86]	 = 0;
			force[2*86+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*86+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*86+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*86+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*86+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(86+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*86+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*86] 		+= temp_force*dx;
			force[2*86+1]	+= temp_force*dy;
		}

			xi = position[2*87];
			yi = position[2*87+1];

			force[2*87]	 = 0;
			force[2*87+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*87+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*87+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*87+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*87+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(87+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*87+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*87] 		+= temp_force*dx;
			force[2*87+1]	+= temp_force*dy;
		}

			xi = position[2*88];
			yi = position[2*88+1];

			force[2*88]	 = 0;
			force[2*88+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*88+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*88+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*88+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*88+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(88+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*88+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*88] 		+= temp_force*dx;
			force[2*88+1]	+= temp_force*dy;
		}

			xi = position[2*89];
			yi = position[2*89+1];

			force[2*89]	 = 0;
			force[2*89+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*89+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*89+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*89+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*89+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(89+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*89+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*89] 		+= temp_force*dx;
			force[2*89+1]	+= temp_force*dy;
		}

			xi = position[2*90];
			yi = position[2*90+1];

			force[2*90]	 = 0;
			force[2*90+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*90+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*90+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*90+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*90+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(90+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*90+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*90] 		+= temp_force*dx;
			force[2*90+1]	+= temp_force*dy;
		}

			xi = position[2*91];
			yi = position[2*91+1];

			force[2*91]	 = 0;
			force[2*91+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*91+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*91+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*91+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*91+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(91+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*91+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*91] 		+= temp_force*dx;
			force[2*91+1]	+= temp_force*dy;
		}

			xi = position[2*92];
			yi = position[2*92+1];

			force[2*92]	 = 0;
			force[2*92+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*92+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*92+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*92+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*92+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(92+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*92+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*92] 		+= temp_force*dx;
			force[2*92+1]	+= temp_force*dy;
		}

			xi = position[2*93];
			yi = position[2*93+1];

			force[2*93]	 = 0;
			force[2*93+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*93+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*93+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*93+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*93+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(93+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*93+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*93] 		+= temp_force*dx;
			force[2*93+1]	+= temp_force*dy;
		}

			xi = position[2*94];
			yi = position[2*94+1];

			force[2*94]	 = 0;
			force[2*94+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*94+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*94+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*94+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*94+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(94+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*94+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*94] 		+= temp_force*dx;
			force[2*94+1]	+= temp_force*dy;
		}

			xi = position[2*95];
			yi = position[2*95+1];

			force[2*95]	 = 0;
			force[2*95+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*95+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*95+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*95+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*95+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(95+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*95+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*95] 		+= temp_force*dx;
			force[2*95+1]	+= temp_force*dy;
		}

			xi = position[2*96];
			yi = position[2*96+1];

			force[2*96]	 = 0;
			force[2*96+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*96+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*96+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*96+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*96+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(96+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*96+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*96] 		+= temp_force*dx;
			force[2*96+1]	+= temp_force*dy;
		}

			xi = position[2*97];
			yi = position[2*97+1];

			force[2*97]	 = 0;
			force[2*97+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*97+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*97+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*97+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*97+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(97+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*97+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*97] 		+= temp_force*dx;
			force[2*97+1]	+= temp_force*dy;
		}

			xi = position[2*98];
			yi = position[2*98+1];

			force[2*98]	 = 0;
			force[2*98+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*98+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*98+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*98+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*98+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(98+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*98+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*98] 		+= temp_force*dx;
			force[2*98+1]	+= temp_force*dy;
		}

			xi = position[2*99];
			yi = position[2*99+1];

			force[2*99]	 = 0;
			force[2*99+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*99+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*99+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*99+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*99+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(99+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*99+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*99] 		+= temp_force*dx;
			force[2*99+1]	+= temp_force*dy;
		}

			xi = position[2*100];
			yi = position[2*100+1];

			force[2*100]	 = 0;
			force[2*100+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*100+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*100+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*100+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*100+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(100+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*100+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*100] 		+= temp_force*dx;
			force[2*100+1]	+= temp_force*dy;
		}

			xi = position[2*101];
			yi = position[2*101+1];

			force[2*101]	 = 0;
			force[2*101+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*101+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*101+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*101+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*101+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(101+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*101+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*101] 		+= temp_force*dx;
			force[2*101+1]	+= temp_force*dy;
		}

			xi = position[2*102];
			yi = position[2*102+1];

			force[2*102]	 = 0;
			force[2*102+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*102+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*102+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*102+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*102+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(102+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*102+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*102] 		+= temp_force*dx;
			force[2*102+1]	+= temp_force*dy;
		}

			xi = position[2*103];
			yi = position[2*103+1];

			force[2*103]	 = 0;
			force[2*103+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*103+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*103+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*103+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*103+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(103+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*103+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*103] 		+= temp_force*dx;
			force[2*103+1]	+= temp_force*dy;
		}

			xi = position[2*104];
			yi = position[2*104+1];

			force[2*104]	 = 0;
			force[2*104+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*104+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*104+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*104+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*104+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(104+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*104+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*104] 		+= temp_force*dx;
			force[2*104+1]	+= temp_force*dy;
		}

			xi = position[2*105];
			yi = position[2*105+1];

			force[2*105]	 = 0;
			force[2*105+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*105+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*105+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*105+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*105+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(105+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*105+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*105] 		+= temp_force*dx;
			force[2*105+1]	+= temp_force*dy;
		}

			xi = position[2*106];
			yi = position[2*106+1];

			force[2*106]	 = 0;
			force[2*106+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*106+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*106+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*106+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*106+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(106+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*106+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*106] 		+= temp_force*dx;
			force[2*106+1]	+= temp_force*dy;
		}

			xi = position[2*107];
			yi = position[2*107+1];

			force[2*107]	 = 0;
			force[2*107+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*107+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*107+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*107+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*107+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(107+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*107+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*107] 		+= temp_force*dx;
			force[2*107+1]	+= temp_force*dy;
		}

			xi = position[2*108];
			yi = position[2*108+1];

			force[2*108]	 = 0;
			force[2*108+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*108+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*108+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*108+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*108+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(108+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*108+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*108] 		+= temp_force*dx;
			force[2*108+1]	+= temp_force*dy;
		}

			xi = position[2*109];
			yi = position[2*109+1];

			force[2*109]	 = 0;
			force[2*109+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*109+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*109+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*109+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*109+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(109+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*109+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*109] 		+= temp_force*dx;
			force[2*109+1]	+= temp_force*dy;
		}

			xi = position[2*110];
			yi = position[2*110+1];

			force[2*110]	 = 0;
			force[2*110+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*110+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*110+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*110+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*110+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(110+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*110+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*110] 		+= temp_force*dx;
			force[2*110+1]	+= temp_force*dy;
		}

			xi = position[2*111];
			yi = position[2*111+1];

			force[2*111]	 = 0;
			force[2*111+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*111+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*111+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*111+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*111+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(111+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*111+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*111] 		+= temp_force*dx;
			force[2*111+1]	+= temp_force*dy;
		}

			xi = position[2*112];
			yi = position[2*112+1];

			force[2*112]	 = 0;
			force[2*112+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*112+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*112+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*112+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*112+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(112+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*112+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*112] 		+= temp_force*dx;
			force[2*112+1]	+= temp_force*dy;
		}

		pthread_barrier_wait(&barrier_internal);

		xi = position[2*75];
		yi = position[2*75+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*75]   + weigh_brown_A * g1 + (position[2*75+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*75+1] + weigh_brown_A * g2;

		position[2*75]		+= dx;
		position[2*75+1] 	+= dy;

		displacement[2*75]	+= dx;
		displacement[2*75+1] += dy;

		verlet_distance[2*75] 	+= (xi - position[2*75]);
		verlet_distance[2*75+1] 	+= (yi - position[2*75+1]);

		if ((temp = sqrt(verlet_distance[2*75]*verlet_distance[2*75] + verlet_distance[2*75+1]*verlet_distance[2*75+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*75] 	-= floor(position[2*75]/L_x)*L_x;

		xi = position[2*76];
		yi = position[2*76+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*76]   + weigh_brown_A * g1 + (position[2*76+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*76+1] + weigh_brown_A * g2;

		position[2*76]		+= dx;
		position[2*76+1] 	+= dy;

		displacement[2*76]	+= dx;
		displacement[2*76+1] += dy;

		verlet_distance[2*76] 	+= (xi - position[2*76]);
		verlet_distance[2*76+1] 	+= (yi - position[2*76+1]);

		if ((temp = sqrt(verlet_distance[2*76]*verlet_distance[2*76] + verlet_distance[2*76+1]*verlet_distance[2*76+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*76] 	-= floor(position[2*76]/L_x)*L_x;

		xi = position[2*77];
		yi = position[2*77+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*77]   + weigh_brown_A * g1 + (position[2*77+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*77+1] + weigh_brown_A * g2;

		position[2*77]		+= dx;
		position[2*77+1] 	+= dy;

		displacement[2*77]	+= dx;
		displacement[2*77+1] += dy;

		verlet_distance[2*77] 	+= (xi - position[2*77]);
		verlet_distance[2*77+1] 	+= (yi - position[2*77+1]);

		if ((temp = sqrt(verlet_distance[2*77]*verlet_distance[2*77] + verlet_distance[2*77+1]*verlet_distance[2*77+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*77] 	-= floor(position[2*77]/L_x)*L_x;

		xi = position[2*78];
		yi = position[2*78+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*78]   + weigh_brown_A * g1 + (position[2*78+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*78+1] + weigh_brown_A * g2;

		position[2*78]		+= dx;
		position[2*78+1] 	+= dy;

		displacement[2*78]	+= dx;
		displacement[2*78+1] += dy;

		verlet_distance[2*78] 	+= (xi - position[2*78]);
		verlet_distance[2*78+1] 	+= (yi - position[2*78+1]);

		if ((temp = sqrt(verlet_distance[2*78]*verlet_distance[2*78] + verlet_distance[2*78+1]*verlet_distance[2*78+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*78] 	-= floor(position[2*78]/L_x)*L_x;

		xi = position[2*79];
		yi = position[2*79+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*79]   + weigh_brown_A * g1 + (position[2*79+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*79+1] + weigh_brown_A * g2;

		position[2*79]		+= dx;
		position[2*79+1] 	+= dy;

		displacement[2*79]	+= dx;
		displacement[2*79+1] += dy;

		verlet_distance[2*79] 	+= (xi - position[2*79]);
		verlet_distance[2*79+1] 	+= (yi - position[2*79+1]);

		if ((temp = sqrt(verlet_distance[2*79]*verlet_distance[2*79] + verlet_distance[2*79+1]*verlet_distance[2*79+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*79] 	-= floor(position[2*79]/L_x)*L_x;

		xi = position[2*80];
		yi = position[2*80+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*80]   + weigh_brown_A * g1 + (position[2*80+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*80+1] + weigh_brown_A * g2;

		position[2*80]		+= dx;
		position[2*80+1] 	+= dy;

		displacement[2*80]	+= dx;
		displacement[2*80+1] += dy;

		verlet_distance[2*80] 	+= (xi - position[2*80]);
		verlet_distance[2*80+1] 	+= (yi - position[2*80+1]);

		if ((temp = sqrt(verlet_distance[2*80]*verlet_distance[2*80] + verlet_distance[2*80+1]*verlet_distance[2*80+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*80] 	-= floor(position[2*80]/L_x)*L_x;

		xi = position[2*81];
		yi = position[2*81+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*81]   + weigh_brown_A * g1 + (position[2*81+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*81+1] + weigh_brown_A * g2;

		position[2*81]		+= dx;
		position[2*81+1] 	+= dy;

		displacement[2*81]	+= dx;
		displacement[2*81+1] += dy;

		verlet_distance[2*81] 	+= (xi - position[2*81]);
		verlet_distance[2*81+1] 	+= (yi - position[2*81+1]);

		if ((temp = sqrt(verlet_distance[2*81]*verlet_distance[2*81] + verlet_distance[2*81+1]*verlet_distance[2*81+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*81] 	-= floor(position[2*81]/L_x)*L_x;

		xi = position[2*82];
		yi = position[2*82+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*82]   + weigh_brown_A * g1 + (position[2*82+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*82+1] + weigh_brown_A * g2;

		position[2*82]		+= dx;
		position[2*82+1] 	+= dy;

		displacement[2*82]	+= dx;
		displacement[2*82+1] += dy;

		verlet_distance[2*82] 	+= (xi - position[2*82]);
		verlet_distance[2*82+1] 	+= (yi - position[2*82+1]);

		if ((temp = sqrt(verlet_distance[2*82]*verlet_distance[2*82] + verlet_distance[2*82+1]*verlet_distance[2*82+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*82] 	-= floor(position[2*82]/L_x)*L_x;

		xi = position[2*83];
		yi = position[2*83+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*83]   + weigh_brown_A * g1 + (position[2*83+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*83+1] + weigh_brown_A * g2;

		position[2*83]		+= dx;
		position[2*83+1] 	+= dy;

		displacement[2*83]	+= dx;
		displacement[2*83+1] += dy;

		verlet_distance[2*83] 	+= (xi - position[2*83]);
		verlet_distance[2*83+1] 	+= (yi - position[2*83+1]);

		if ((temp = sqrt(verlet_distance[2*83]*verlet_distance[2*83] + verlet_distance[2*83+1]*verlet_distance[2*83+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*83] 	-= floor(position[2*83]/L_x)*L_x;

		xi = position[2*84];
		yi = position[2*84+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*84]   + weigh_brown_A * g1 + (position[2*84+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*84+1] + weigh_brown_A * g2;

		position[2*84]		+= dx;
		position[2*84+1] 	+= dy;

		displacement[2*84]	+= dx;
		displacement[2*84+1] += dy;

		verlet_distance[2*84] 	+= (xi - position[2*84]);
		verlet_distance[2*84+1] 	+= (yi - position[2*84+1]);

		if ((temp = sqrt(verlet_distance[2*84]*verlet_distance[2*84] + verlet_distance[2*84+1]*verlet_distance[2*84+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*84] 	-= floor(position[2*84]/L_x)*L_x;

		xi = position[2*85];
		yi = position[2*85+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*85]   + weigh_brown_A * g1 + (position[2*85+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*85+1] + weigh_brown_A * g2;

		position[2*85]		+= dx;
		position[2*85+1] 	+= dy;

		displacement[2*85]	+= dx;
		displacement[2*85+1] += dy;

		verlet_distance[2*85] 	+= (xi - position[2*85]);
		verlet_distance[2*85+1] 	+= (yi - position[2*85+1]);

		if ((temp = sqrt(verlet_distance[2*85]*verlet_distance[2*85] + verlet_distance[2*85+1]*verlet_distance[2*85+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*85] 	-= floor(position[2*85]/L_x)*L_x;

		xi = position[2*86];
		yi = position[2*86+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*86]   + weigh_brown_A * g1 + (position[2*86+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*86+1] + weigh_brown_A * g2;

		position[2*86]		+= dx;
		position[2*86+1] 	+= dy;

		displacement[2*86]	+= dx;
		displacement[2*86+1] += dy;

		verlet_distance[2*86] 	+= (xi - position[2*86]);
		verlet_distance[2*86+1] 	+= (yi - position[2*86+1]);

		if ((temp = sqrt(verlet_distance[2*86]*verlet_distance[2*86] + verlet_distance[2*86+1]*verlet_distance[2*86+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*86] 	-= floor(position[2*86]/L_x)*L_x;

		xi = position[2*87];
		yi = position[2*87+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*87]   + weigh_brown_A * g1 + (position[2*87+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*87+1] + weigh_brown_A * g2;

		position[2*87]		+= dx;
		position[2*87+1] 	+= dy;

		displacement[2*87]	+= dx;
		displacement[2*87+1] += dy;

		verlet_distance[2*87] 	+= (xi - position[2*87]);
		verlet_distance[2*87+1] 	+= (yi - position[2*87+1]);

		if ((temp = sqrt(verlet_distance[2*87]*verlet_distance[2*87] + verlet_distance[2*87+1]*verlet_distance[2*87+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*87] 	-= floor(position[2*87]/L_x)*L_x;

		xi = position[2*88];
		yi = position[2*88+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*88]   + weigh_brown_A * g1 + (position[2*88+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*88+1] + weigh_brown_A * g2;

		position[2*88]		+= dx;
		position[2*88+1] 	+= dy;

		displacement[2*88]	+= dx;
		displacement[2*88+1] += dy;

		verlet_distance[2*88] 	+= (xi - position[2*88]);
		verlet_distance[2*88+1] 	+= (yi - position[2*88+1]);

		if ((temp = sqrt(verlet_distance[2*88]*verlet_distance[2*88] + verlet_distance[2*88+1]*verlet_distance[2*88+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*88] 	-= floor(position[2*88]/L_x)*L_x;

		xi = position[2*89];
		yi = position[2*89+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*89]   + weigh_brown_A * g1 + (position[2*89+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*89+1] + weigh_brown_A * g2;

		position[2*89]		+= dx;
		position[2*89+1] 	+= dy;

		displacement[2*89]	+= dx;
		displacement[2*89+1] += dy;

		verlet_distance[2*89] 	+= (xi - position[2*89]);
		verlet_distance[2*89+1] 	+= (yi - position[2*89+1]);

		if ((temp = sqrt(verlet_distance[2*89]*verlet_distance[2*89] + verlet_distance[2*89+1]*verlet_distance[2*89+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*89] 	-= floor(position[2*89]/L_x)*L_x;

		xi = position[2*90];
		yi = position[2*90+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*90]   + weigh_brown_A * g1 + (position[2*90+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*90+1] + weigh_brown_A * g2;

		position[2*90]		+= dx;
		position[2*90+1] 	+= dy;

		displacement[2*90]	+= dx;
		displacement[2*90+1] += dy;

		verlet_distance[2*90] 	+= (xi - position[2*90]);
		verlet_distance[2*90+1] 	+= (yi - position[2*90+1]);

		if ((temp = sqrt(verlet_distance[2*90]*verlet_distance[2*90] + verlet_distance[2*90+1]*verlet_distance[2*90+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*90] 	-= floor(position[2*90]/L_x)*L_x;

		xi = position[2*91];
		yi = position[2*91+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*91]   + weigh_brown_A * g1 + (position[2*91+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*91+1] + weigh_brown_A * g2;

		position[2*91]		+= dx;
		position[2*91+1] 	+= dy;

		displacement[2*91]	+= dx;
		displacement[2*91+1] += dy;

		verlet_distance[2*91] 	+= (xi - position[2*91]);
		verlet_distance[2*91+1] 	+= (yi - position[2*91+1]);

		if ((temp = sqrt(verlet_distance[2*91]*verlet_distance[2*91] + verlet_distance[2*91+1]*verlet_distance[2*91+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*91] 	-= floor(position[2*91]/L_x)*L_x;

		xi = position[2*92];
		yi = position[2*92+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*92]   + weigh_brown_A * g1 + (position[2*92+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*92+1] + weigh_brown_A * g2;

		position[2*92]		+= dx;
		position[2*92+1] 	+= dy;

		displacement[2*92]	+= dx;
		displacement[2*92+1] += dy;

		verlet_distance[2*92] 	+= (xi - position[2*92]);
		verlet_distance[2*92+1] 	+= (yi - position[2*92+1]);

		if ((temp = sqrt(verlet_distance[2*92]*verlet_distance[2*92] + verlet_distance[2*92+1]*verlet_distance[2*92+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*92] 	-= floor(position[2*92]/L_x)*L_x;

		xi = position[2*93];
		yi = position[2*93+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*93]   + weigh_brown_A * g1 + (position[2*93+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*93+1] + weigh_brown_A * g2;

		position[2*93]		+= dx;
		position[2*93+1] 	+= dy;

		displacement[2*93]	+= dx;
		displacement[2*93+1] += dy;

		verlet_distance[2*93] 	+= (xi - position[2*93]);
		verlet_distance[2*93+1] 	+= (yi - position[2*93+1]);

		if ((temp = sqrt(verlet_distance[2*93]*verlet_distance[2*93] + verlet_distance[2*93+1]*verlet_distance[2*93+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*93] 	-= floor(position[2*93]/L_x)*L_x;

		xi = position[2*94];
		yi = position[2*94+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*94]   + weigh_brown_A * g1 + (position[2*94+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*94+1] + weigh_brown_A * g2;

		position[2*94]		+= dx;
		position[2*94+1] 	+= dy;

		displacement[2*94]	+= dx;
		displacement[2*94+1] += dy;

		verlet_distance[2*94] 	+= (xi - position[2*94]);
		verlet_distance[2*94+1] 	+= (yi - position[2*94+1]);

		if ((temp = sqrt(verlet_distance[2*94]*verlet_distance[2*94] + verlet_distance[2*94+1]*verlet_distance[2*94+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*94] 	-= floor(position[2*94]/L_x)*L_x;

		xi = position[2*95];
		yi = position[2*95+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*95]   + weigh_brown_A * g1 + (position[2*95+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*95+1] + weigh_brown_A * g2;

		position[2*95]		+= dx;
		position[2*95+1] 	+= dy;

		displacement[2*95]	+= dx;
		displacement[2*95+1] += dy;

		verlet_distance[2*95] 	+= (xi - position[2*95]);
		verlet_distance[2*95+1] 	+= (yi - position[2*95+1]);

		if ((temp = sqrt(verlet_distance[2*95]*verlet_distance[2*95] + verlet_distance[2*95+1]*verlet_distance[2*95+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*95] 	-= floor(position[2*95]/L_x)*L_x;

		xi = position[2*96];
		yi = position[2*96+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*96]   + weigh_brown_A * g1 + (position[2*96+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*96+1] + weigh_brown_A * g2;

		position[2*96]		+= dx;
		position[2*96+1] 	+= dy;

		displacement[2*96]	+= dx;
		displacement[2*96+1] += dy;

		verlet_distance[2*96] 	+= (xi - position[2*96]);
		verlet_distance[2*96+1] 	+= (yi - position[2*96+1]);

		if ((temp = sqrt(verlet_distance[2*96]*verlet_distance[2*96] + verlet_distance[2*96+1]*verlet_distance[2*96+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*96] 	-= floor(position[2*96]/L_x)*L_x;

		xi = position[2*97];
		yi = position[2*97+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*97]   + weigh_brown_A * g1 + (position[2*97+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*97+1] + weigh_brown_A * g2;

		position[2*97]		+= dx;
		position[2*97+1] 	+= dy;

		displacement[2*97]	+= dx;
		displacement[2*97+1] += dy;

		verlet_distance[2*97] 	+= (xi - position[2*97]);
		verlet_distance[2*97+1] 	+= (yi - position[2*97+1]);

		if ((temp = sqrt(verlet_distance[2*97]*verlet_distance[2*97] + verlet_distance[2*97+1]*verlet_distance[2*97+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*97] 	-= floor(position[2*97]/L_x)*L_x;

		xi = position[2*98];
		yi = position[2*98+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*98]   + weigh_brown_A * g1 + (position[2*98+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*98+1] + weigh_brown_A * g2;

		position[2*98]		+= dx;
		position[2*98+1] 	+= dy;

		displacement[2*98]	+= dx;
		displacement[2*98+1] += dy;

		verlet_distance[2*98] 	+= (xi - position[2*98]);
		verlet_distance[2*98+1] 	+= (yi - position[2*98+1]);

		if ((temp = sqrt(verlet_distance[2*98]*verlet_distance[2*98] + verlet_distance[2*98+1]*verlet_distance[2*98+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*98] 	-= floor(position[2*98]/L_x)*L_x;

		xi = position[2*99];
		yi = position[2*99+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*99]   + weigh_brown_A * g1 + (position[2*99+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*99+1] + weigh_brown_A * g2;

		position[2*99]		+= dx;
		position[2*99+1] 	+= dy;

		displacement[2*99]	+= dx;
		displacement[2*99+1] += dy;

		verlet_distance[2*99] 	+= (xi - position[2*99]);
		verlet_distance[2*99+1] 	+= (yi - position[2*99+1]);

		if ((temp = sqrt(verlet_distance[2*99]*verlet_distance[2*99] + verlet_distance[2*99+1]*verlet_distance[2*99+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*99] 	-= floor(position[2*99]/L_x)*L_x;

		xi = position[2*100];
		yi = position[2*100+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*100]   + weigh_brown_A * g1 + (position[2*100+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*100+1] + weigh_brown_A * g2;

		position[2*100]		+= dx;
		position[2*100+1] 	+= dy;

		displacement[2*100]	+= dx;
		displacement[2*100+1] += dy;

		verlet_distance[2*100] 	+= (xi - position[2*100]);
		verlet_distance[2*100+1] 	+= (yi - position[2*100+1]);

		if ((temp = sqrt(verlet_distance[2*100]*verlet_distance[2*100] + verlet_distance[2*100+1]*verlet_distance[2*100+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*100] 	-= floor(position[2*100]/L_x)*L_x;

		xi = position[2*101];
		yi = position[2*101+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*101]   + weigh_brown_A * g1 + (position[2*101+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*101+1] + weigh_brown_A * g2;

		position[2*101]		+= dx;
		position[2*101+1] 	+= dy;

		displacement[2*101]	+= dx;
		displacement[2*101+1] += dy;

		verlet_distance[2*101] 	+= (xi - position[2*101]);
		verlet_distance[2*101+1] 	+= (yi - position[2*101+1]);

		if ((temp = sqrt(verlet_distance[2*101]*verlet_distance[2*101] + verlet_distance[2*101+1]*verlet_distance[2*101+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*101] 	-= floor(position[2*101]/L_x)*L_x;

		xi = position[2*102];
		yi = position[2*102+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*102]   + weigh_brown_A * g1 + (position[2*102+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*102+1] + weigh_brown_A * g2;

		position[2*102]		+= dx;
		position[2*102+1] 	+= dy;

		displacement[2*102]	+= dx;
		displacement[2*102+1] += dy;

		verlet_distance[2*102] 	+= (xi - position[2*102]);
		verlet_distance[2*102+1] 	+= (yi - position[2*102+1]);

		if ((temp = sqrt(verlet_distance[2*102]*verlet_distance[2*102] + verlet_distance[2*102+1]*verlet_distance[2*102+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*102] 	-= floor(position[2*102]/L_x)*L_x;

		xi = position[2*103];
		yi = position[2*103+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*103]   + weigh_brown_A * g1 + (position[2*103+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*103+1] + weigh_brown_A * g2;

		position[2*103]		+= dx;
		position[2*103+1] 	+= dy;

		displacement[2*103]	+= dx;
		displacement[2*103+1] += dy;

		verlet_distance[2*103] 	+= (xi - position[2*103]);
		verlet_distance[2*103+1] 	+= (yi - position[2*103+1]);

		if ((temp = sqrt(verlet_distance[2*103]*verlet_distance[2*103] + verlet_distance[2*103+1]*verlet_distance[2*103+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*103] 	-= floor(position[2*103]/L_x)*L_x;

		xi = position[2*104];
		yi = position[2*104+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*104]   + weigh_brown_A * g1 + (position[2*104+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*104+1] + weigh_brown_A * g2;

		position[2*104]		+= dx;
		position[2*104+1] 	+= dy;

		displacement[2*104]	+= dx;
		displacement[2*104+1] += dy;

		verlet_distance[2*104] 	+= (xi - position[2*104]);
		verlet_distance[2*104+1] 	+= (yi - position[2*104+1]);

		if ((temp = sqrt(verlet_distance[2*104]*verlet_distance[2*104] + verlet_distance[2*104+1]*verlet_distance[2*104+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*104] 	-= floor(position[2*104]/L_x)*L_x;

		xi = position[2*105];
		yi = position[2*105+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*105]   + weigh_brown_A * g1 + (position[2*105+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*105+1] + weigh_brown_A * g2;

		position[2*105]		+= dx;
		position[2*105+1] 	+= dy;

		displacement[2*105]	+= dx;
		displacement[2*105+1] += dy;

		verlet_distance[2*105] 	+= (xi - position[2*105]);
		verlet_distance[2*105+1] 	+= (yi - position[2*105+1]);

		if ((temp = sqrt(verlet_distance[2*105]*verlet_distance[2*105] + verlet_distance[2*105+1]*verlet_distance[2*105+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*105] 	-= floor(position[2*105]/L_x)*L_x;

		xi = position[2*106];
		yi = position[2*106+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*106]   + weigh_brown_A * g1 + (position[2*106+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*106+1] + weigh_brown_A * g2;

		position[2*106]		+= dx;
		position[2*106+1] 	+= dy;

		displacement[2*106]	+= dx;
		displacement[2*106+1] += dy;

		verlet_distance[2*106] 	+= (xi - position[2*106]);
		verlet_distance[2*106+1] 	+= (yi - position[2*106+1]);

		if ((temp = sqrt(verlet_distance[2*106]*verlet_distance[2*106] + verlet_distance[2*106+1]*verlet_distance[2*106+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*106] 	-= floor(position[2*106]/L_x)*L_x;

		xi = position[2*107];
		yi = position[2*107+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*107]   + weigh_brown_A * g1 + (position[2*107+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*107+1] + weigh_brown_A * g2;

		position[2*107]		+= dx;
		position[2*107+1] 	+= dy;

		displacement[2*107]	+= dx;
		displacement[2*107+1] += dy;

		verlet_distance[2*107] 	+= (xi - position[2*107]);
		verlet_distance[2*107+1] 	+= (yi - position[2*107+1]);

		if ((temp = sqrt(verlet_distance[2*107]*verlet_distance[2*107] + verlet_distance[2*107+1]*verlet_distance[2*107+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*107] 	-= floor(position[2*107]/L_x)*L_x;

		xi = position[2*108];
		yi = position[2*108+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*108]   + weigh_brown_A * g1 + (position[2*108+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*108+1] + weigh_brown_A * g2;

		position[2*108]		+= dx;
		position[2*108+1] 	+= dy;

		displacement[2*108]	+= dx;
		displacement[2*108+1] += dy;

		verlet_distance[2*108] 	+= (xi - position[2*108]);
		verlet_distance[2*108+1] 	+= (yi - position[2*108+1]);

		if ((temp = sqrt(verlet_distance[2*108]*verlet_distance[2*108] + verlet_distance[2*108+1]*verlet_distance[2*108+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*108] 	-= floor(position[2*108]/L_x)*L_x;

		xi = position[2*109];
		yi = position[2*109+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*109]   + weigh_brown_A * g1 + (position[2*109+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*109+1] + weigh_brown_A * g2;

		position[2*109]		+= dx;
		position[2*109+1] 	+= dy;

		displacement[2*109]	+= dx;
		displacement[2*109+1] += dy;

		verlet_distance[2*109] 	+= (xi - position[2*109]);
		verlet_distance[2*109+1] 	+= (yi - position[2*109+1]);

		if ((temp = sqrt(verlet_distance[2*109]*verlet_distance[2*109] + verlet_distance[2*109+1]*verlet_distance[2*109+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*109] 	-= floor(position[2*109]/L_x)*L_x;

		xi = position[2*110];
		yi = position[2*110+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*110]   + weigh_brown_A * g1 + (position[2*110+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*110+1] + weigh_brown_A * g2;

		position[2*110]		+= dx;
		position[2*110+1] 	+= dy;

		displacement[2*110]	+= dx;
		displacement[2*110+1] += dy;

		verlet_distance[2*110] 	+= (xi - position[2*110]);
		verlet_distance[2*110+1] 	+= (yi - position[2*110+1]);

		if ((temp = sqrt(verlet_distance[2*110]*verlet_distance[2*110] + verlet_distance[2*110+1]*verlet_distance[2*110+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*110] 	-= floor(position[2*110]/L_x)*L_x;

		xi = position[2*111];
		yi = position[2*111+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*111]   + weigh_brown_A * g1 + (position[2*111+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*111+1] + weigh_brown_A * g2;

		position[2*111]		+= dx;
		position[2*111+1] 	+= dy;

		displacement[2*111]	+= dx;
		displacement[2*111+1] += dy;

		verlet_distance[2*111] 	+= (xi - position[2*111]);
		verlet_distance[2*111+1] 	+= (yi - position[2*111+1]);

		if ((temp = sqrt(verlet_distance[2*111]*verlet_distance[2*111] + verlet_distance[2*111+1]*verlet_distance[2*111+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*111] 	-= floor(position[2*111]/L_x)*L_x;

		xi = position[2*112];
		yi = position[2*112+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*112]   + weigh_brown_A * g1 + (position[2*112+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*112+1] + weigh_brown_A * g2;

		position[2*112]		+= dx;
		position[2*112+1] 	+= dy;

		displacement[2*112]	+= dx;
		displacement[2*112+1] += dy;

		verlet_distance[2*112] 	+= (xi - position[2*112]);
		verlet_distance[2*112+1] 	+= (yi - position[2*112+1]);

		if ((temp = sqrt(verlet_distance[2*112]*verlet_distance[2*112] + verlet_distance[2*112+1]*verlet_distance[2*112+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*112] 	-= floor(position[2*112]/L_x)*L_x;

		// Signal to main thread, that all threads have finished their iteration
		pthread_barrier_wait(&barrier_main_one);

		// thread is done with one iteration, it will now wait for a signal from the main thread to continue
		pthread_barrier_wait(&barrier_main_two);
	}

	return NULL;
}



static void *iteration_3 (int *no) {
	// Define necessary variables
	double xi, xj, yi, yj;
	double dx, dy;
	double r_squared;

	double temp_force;
	double temp;

	double u1, u2, g1, g2;

	int iterate;
	int j;

	double m_i_A = 1.0;
	double m_i_B = m;

	double D_kT_A = D_Brown_A/kT;
	double D_kT_B = D_Brown_B/kT;

	double shear_A = v_A * delta_t;

	double m_j;

	double dist_bottom;
	double dist_top;

	// cutoff for the walls force, this will lead to a smoother transition
	double dist_cutoff 			= 1.0;
	double prefactor			= 3.0;
	double force_wall_cutoff 	= 1.0*prefactor *(kappa/dist_cutoff + 1.0/(dist_cutoff*dist_cutoff))*exp(-kappa*dist_cutoff);

	while(cont == 1) {
			xi = position[2*113];
			yi = position[2*113+1];

			force[2*113]	 = 0;
			force[2*113+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*113+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*113+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*113+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*113+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(113+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*113+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*113] 		+= temp_force*dx;
			force[2*113+1]	+= temp_force*dy;
		}

			xi = position[2*114];
			yi = position[2*114+1];

			force[2*114]	 = 0;
			force[2*114+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*114+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*114+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*114+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*114+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(114+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*114+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*114] 		+= temp_force*dx;
			force[2*114+1]	+= temp_force*dy;
		}

			xi = position[2*115];
			yi = position[2*115+1];

			force[2*115]	 = 0;
			force[2*115+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*115+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*115+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*115+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*115+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(115+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*115+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*115] 		+= temp_force*dx;
			force[2*115+1]	+= temp_force*dy;
		}

			xi = position[2*116];
			yi = position[2*116+1];

			force[2*116]	 = 0;
			force[2*116+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*116+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*116+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*116+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*116+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(116+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*116+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*116] 		+= temp_force*dx;
			force[2*116+1]	+= temp_force*dy;
		}

			xi = position[2*117];
			yi = position[2*117+1];

			force[2*117]	 = 0;
			force[2*117+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*117+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*117+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*117+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*117+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(117+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*117+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*117] 		+= temp_force*dx;
			force[2*117+1]	+= temp_force*dy;
		}

			xi = position[2*118];
			yi = position[2*118+1];

			force[2*118]	 = 0;
			force[2*118+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*118+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*118+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*118+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*118+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(118+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*118+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*118] 		+= temp_force*dx;
			force[2*118+1]	+= temp_force*dy;
		}

			xi = position[2*119];
			yi = position[2*119+1];

			force[2*119]	 = 0;
			force[2*119+1] = 0;

			if (yi <= 1e-12) {
				dist_bottom = 1e-12;
				force[2*119+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}
			else if (yi <= dist_cutoff) {
				dist_bottom = yi;
				force[2*119+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
			}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*119+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*119+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

			iterate = verlet[N*(119+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*119+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_A*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*119] 		+= temp_force*dx;
			force[2*119+1]	+= temp_force*dy;
		}

		xi = position[2*120];
		yi = position[2*120+1];

		force[2*120]	 = 0;
		force[2*120+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*120+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*120+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*120+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*120+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(120+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*120+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*120] 		+= temp_force*dx;
			force[2*120+1]	+= temp_force*dy;
		}

		xi = position[2*121];
		yi = position[2*121+1];

		force[2*121]	 = 0;
		force[2*121+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*121+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*121+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*121+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*121+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(121+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*121+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*121] 		+= temp_force*dx;
			force[2*121+1]	+= temp_force*dy;
		}

		xi = position[2*122];
		yi = position[2*122+1];

		force[2*122]	 = 0;
		force[2*122+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*122+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*122+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*122+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*122+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(122+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*122+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*122] 		+= temp_force*dx;
			force[2*122+1]	+= temp_force*dy;
		}

		xi = position[2*123];
		yi = position[2*123+1];

		force[2*123]	 = 0;
		force[2*123+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*123+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*123+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*123+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*123+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(123+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*123+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*123] 		+= temp_force*dx;
			force[2*123+1]	+= temp_force*dy;
		}

		xi = position[2*124];
		yi = position[2*124+1];

		force[2*124]	 = 0;
		force[2*124+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*124+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*124+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*124+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*124+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(124+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*124+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*124] 		+= temp_force*dx;
			force[2*124+1]	+= temp_force*dy;
		}

		xi = position[2*125];
		yi = position[2*125+1];

		force[2*125]	 = 0;
		force[2*125+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*125+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*125+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*125+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*125+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(125+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*125+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*125] 		+= temp_force*dx;
			force[2*125+1]	+= temp_force*dy;
		}

		xi = position[2*126];
		yi = position[2*126+1];

		force[2*126]	 = 0;
		force[2*126+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*126+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*126+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*126+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*126+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(126+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*126+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*126] 		+= temp_force*dx;
			force[2*126+1]	+= temp_force*dy;
		}

		xi = position[2*127];
		yi = position[2*127+1];

		force[2*127]	 = 0;
		force[2*127+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*127+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*127+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*127+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*127+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(127+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*127+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*127] 		+= temp_force*dx;
			force[2*127+1]	+= temp_force*dy;
		}

		xi = position[2*128];
		yi = position[2*128+1];

		force[2*128]	 = 0;
		force[2*128+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*128+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*128+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*128+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*128+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(128+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*128+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*128] 		+= temp_force*dx;
			force[2*128+1]	+= temp_force*dy;
		}

		xi = position[2*129];
		yi = position[2*129+1];

		force[2*129]	 = 0;
		force[2*129+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*129+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*129+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*129+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*129+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(129+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*129+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*129] 		+= temp_force*dx;
			force[2*129+1]	+= temp_force*dy;
		}

		xi = position[2*130];
		yi = position[2*130+1];

		force[2*130]	 = 0;
		force[2*130+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*130+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*130+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*130+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*130+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(130+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*130+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*130] 		+= temp_force*dx;
			force[2*130+1]	+= temp_force*dy;
		}

		xi = position[2*131];
		yi = position[2*131+1];

		force[2*131]	 = 0;
		force[2*131+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*131+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*131+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*131+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*131+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(131+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*131+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*131] 		+= temp_force*dx;
			force[2*131+1]	+= temp_force*dy;
		}

		xi = position[2*132];
		yi = position[2*132+1];

		force[2*132]	 = 0;
		force[2*132+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*132+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*132+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*132+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*132+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(132+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*132+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*132] 		+= temp_force*dx;
			force[2*132+1]	+= temp_force*dy;
		}

		xi = position[2*133];
		yi = position[2*133+1];

		force[2*133]	 = 0;
		force[2*133+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*133+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*133+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*133+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*133+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(133+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*133+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*133] 		+= temp_force*dx;
			force[2*133+1]	+= temp_force*dy;
		}

		xi = position[2*134];
		yi = position[2*134+1];

		force[2*134]	 = 0;
		force[2*134+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*134+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*134+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*134+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*134+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(134+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*134+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*134] 		+= temp_force*dx;
			force[2*134+1]	+= temp_force*dy;
		}

		xi = position[2*135];
		yi = position[2*135+1];

		force[2*135]	 = 0;
		force[2*135+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*135+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*135+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*135+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*135+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(135+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*135+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*135] 		+= temp_force*dx;
			force[2*135+1]	+= temp_force*dy;
		}

		xi = position[2*136];
		yi = position[2*136+1];

		force[2*136]	 = 0;
		force[2*136+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*136+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*136+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*136+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*136+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(136+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*136+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*136] 		+= temp_force*dx;
			force[2*136+1]	+= temp_force*dy;
		}

		xi = position[2*137];
		yi = position[2*137+1];

		force[2*137]	 = 0;
		force[2*137+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*137+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*137+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*137+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*137+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(137+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*137+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*137] 		+= temp_force*dx;
			force[2*137+1]	+= temp_force*dy;
		}

		xi = position[2*138];
		yi = position[2*138+1];

		force[2*138]	 = 0;
		force[2*138+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*138+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*138+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*138+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*138+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(138+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*138+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*138] 		+= temp_force*dx;
			force[2*138+1]	+= temp_force*dy;
		}

		xi = position[2*139];
		yi = position[2*139+1];

		force[2*139]	 = 0;
		force[2*139+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*139+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*139+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*139+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*139+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(139+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*139+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*139] 		+= temp_force*dx;
			force[2*139+1]	+= temp_force*dy;
		}

		xi = position[2*140];
		yi = position[2*140+1];

		force[2*140]	 = 0;
		force[2*140+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*140+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*140+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*140+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*140+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(140+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*140+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*140] 		+= temp_force*dx;
			force[2*140+1]	+= temp_force*dy;
		}

		xi = position[2*141];
		yi = position[2*141+1];

		force[2*141]	 = 0;
		force[2*141+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*141+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*141+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*141+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*141+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(141+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*141+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*141] 		+= temp_force*dx;
			force[2*141+1]	+= temp_force*dy;
		}

		xi = position[2*142];
		yi = position[2*142+1];

		force[2*142]	 = 0;
		force[2*142+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*142+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*142+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*142+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*142+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(142+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*142+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*142] 		+= temp_force*dx;
			force[2*142+1]	+= temp_force*dy;
		}

		xi = position[2*143];
		yi = position[2*143+1];

		force[2*143]	 = 0;
		force[2*143+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*143+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*143+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*143+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*143+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(143+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*143+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*143] 		+= temp_force*dx;
			force[2*143+1]	+= temp_force*dy;
		}

		xi = position[2*144];
		yi = position[2*144+1];

		force[2*144]	 = 0;
		force[2*144+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*144+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*144+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*144+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*144+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(144+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*144+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*144] 		+= temp_force*dx;
			force[2*144+1]	+= temp_force*dy;
		}

		xi = position[2*145];
		yi = position[2*145+1];

		force[2*145]	 = 0;
		force[2*145+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*145+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*145+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*145+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*145+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(145+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*145+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*145] 		+= temp_force*dx;
			force[2*145+1]	+= temp_force*dy;
		}

		xi = position[2*146];
		yi = position[2*146+1];

		force[2*146]	 = 0;
		force[2*146+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*146+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*146+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*146+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*146+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(146+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*146+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*146] 		+= temp_force*dx;
			force[2*146+1]	+= temp_force*dy;
		}

		xi = position[2*147];
		yi = position[2*147+1];

		force[2*147]	 = 0;
		force[2*147+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*147+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*147+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*147+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*147+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(147+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*147+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*147] 		+= temp_force*dx;
			force[2*147+1]	+= temp_force*dy;
		}

		xi = position[2*148];
		yi = position[2*148+1];

		force[2*148]	 = 0;
		force[2*148+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*148+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*148+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*148+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*148+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(148+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*148+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*148] 		+= temp_force*dx;
			force[2*148+1]	+= temp_force*dy;
		}

		xi = position[2*149];
		yi = position[2*149+1];

		force[2*149]	 = 0;
		force[2*149+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*149+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*149+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*149+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*149+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(149+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*149+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*149] 		+= temp_force*dx;
			force[2*149+1]	+= temp_force*dy;
		}

		pthread_barrier_wait(&barrier_internal);

		xi = position[2*113];
		yi = position[2*113+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*113]   + weigh_brown_A * g1 + (position[2*113+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*113+1] + weigh_brown_A * g2;

		position[2*113]		+= dx;
		position[2*113+1] 	+= dy;

		displacement[2*113]	+= dx;
		displacement[2*113+1] += dy;

		verlet_distance[2*113] 	+= (xi - position[2*113]);
		verlet_distance[2*113+1] 	+= (yi - position[2*113+1]);

		if ((temp = sqrt(verlet_distance[2*113]*verlet_distance[2*113] + verlet_distance[2*113+1]*verlet_distance[2*113+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*113] 	-= floor(position[2*113]/L_x)*L_x;

		xi = position[2*114];
		yi = position[2*114+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*114]   + weigh_brown_A * g1 + (position[2*114+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*114+1] + weigh_brown_A * g2;

		position[2*114]		+= dx;
		position[2*114+1] 	+= dy;

		displacement[2*114]	+= dx;
		displacement[2*114+1] += dy;

		verlet_distance[2*114] 	+= (xi - position[2*114]);
		verlet_distance[2*114+1] 	+= (yi - position[2*114+1]);

		if ((temp = sqrt(verlet_distance[2*114]*verlet_distance[2*114] + verlet_distance[2*114+1]*verlet_distance[2*114+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*114] 	-= floor(position[2*114]/L_x)*L_x;

		xi = position[2*115];
		yi = position[2*115+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*115]   + weigh_brown_A * g1 + (position[2*115+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*115+1] + weigh_brown_A * g2;

		position[2*115]		+= dx;
		position[2*115+1] 	+= dy;

		displacement[2*115]	+= dx;
		displacement[2*115+1] += dy;

		verlet_distance[2*115] 	+= (xi - position[2*115]);
		verlet_distance[2*115+1] 	+= (yi - position[2*115+1]);

		if ((temp = sqrt(verlet_distance[2*115]*verlet_distance[2*115] + verlet_distance[2*115+1]*verlet_distance[2*115+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*115] 	-= floor(position[2*115]/L_x)*L_x;

		xi = position[2*116];
		yi = position[2*116+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*116]   + weigh_brown_A * g1 + (position[2*116+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*116+1] + weigh_brown_A * g2;

		position[2*116]		+= dx;
		position[2*116+1] 	+= dy;

		displacement[2*116]	+= dx;
		displacement[2*116+1] += dy;

		verlet_distance[2*116] 	+= (xi - position[2*116]);
		verlet_distance[2*116+1] 	+= (yi - position[2*116+1]);

		if ((temp = sqrt(verlet_distance[2*116]*verlet_distance[2*116] + verlet_distance[2*116+1]*verlet_distance[2*116+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*116] 	-= floor(position[2*116]/L_x)*L_x;

		xi = position[2*117];
		yi = position[2*117+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*117]   + weigh_brown_A * g1 + (position[2*117+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*117+1] + weigh_brown_A * g2;

		position[2*117]		+= dx;
		position[2*117+1] 	+= dy;

		displacement[2*117]	+= dx;
		displacement[2*117+1] += dy;

		verlet_distance[2*117] 	+= (xi - position[2*117]);
		verlet_distance[2*117+1] 	+= (yi - position[2*117+1]);

		if ((temp = sqrt(verlet_distance[2*117]*verlet_distance[2*117] + verlet_distance[2*117+1]*verlet_distance[2*117+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*117] 	-= floor(position[2*117]/L_x)*L_x;

		xi = position[2*118];
		yi = position[2*118+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*118]   + weigh_brown_A * g1 + (position[2*118+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*118+1] + weigh_brown_A * g2;

		position[2*118]		+= dx;
		position[2*118+1] 	+= dy;

		displacement[2*118]	+= dx;
		displacement[2*118+1] += dy;

		verlet_distance[2*118] 	+= (xi - position[2*118]);
		verlet_distance[2*118+1] 	+= (yi - position[2*118+1]);

		if ((temp = sqrt(verlet_distance[2*118]*verlet_distance[2*118] + verlet_distance[2*118+1]*verlet_distance[2*118+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*118] 	-= floor(position[2*118]/L_x)*L_x;

		xi = position[2*119];
		yi = position[2*119+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_A*delta_t*force[2*119]   + weigh_brown_A * g1 + (position[2*119+1] - 0.5*L_y)*shear_A;
		dy = D_kT_A*delta_t*force[2*119+1] + weigh_brown_A * g2;

		position[2*119]		+= dx;
		position[2*119+1] 	+= dy;

		displacement[2*119]	+= dx;
		displacement[2*119+1] += dy;

		verlet_distance[2*119] 	+= (xi - position[2*119]);
		verlet_distance[2*119+1] 	+= (yi - position[2*119+1]);

		if ((temp = sqrt(verlet_distance[2*119]*verlet_distance[2*119] + verlet_distance[2*119+1]*verlet_distance[2*119+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*119] 	-= floor(position[2*119]/L_x)*L_x;

		xi = position[2*120];
		yi = position[2*120+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*120]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*120+1] + weigh_brown_B * g2;

		position[2*120]		+= dx;
		position[2*120+1] 	+= dy;

		displacement[2*120]	+= dx;
		displacement[2*120+1] += dy;

		verlet_distance[2*120] 	+= (xi - position[2*120]);
		verlet_distance[2*120+1] 	+= (yi - position[2*120+1]);

		if ((temp = sqrt(verlet_distance[2*120]*verlet_distance[2*120] + verlet_distance[2*120+1]*verlet_distance[2*120+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*120] 	-= floor(position[2*120]/L_x)*L_x;

		xi = position[2*121];
		yi = position[2*121+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*121]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*121+1] + weigh_brown_B * g2;

		position[2*121]		+= dx;
		position[2*121+1] 	+= dy;

		displacement[2*121]	+= dx;
		displacement[2*121+1] += dy;

		verlet_distance[2*121] 	+= (xi - position[2*121]);
		verlet_distance[2*121+1] 	+= (yi - position[2*121+1]);

		if ((temp = sqrt(verlet_distance[2*121]*verlet_distance[2*121] + verlet_distance[2*121+1]*verlet_distance[2*121+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*121] 	-= floor(position[2*121]/L_x)*L_x;

		xi = position[2*122];
		yi = position[2*122+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*122]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*122+1] + weigh_brown_B * g2;

		position[2*122]		+= dx;
		position[2*122+1] 	+= dy;

		displacement[2*122]	+= dx;
		displacement[2*122+1] += dy;

		verlet_distance[2*122] 	+= (xi - position[2*122]);
		verlet_distance[2*122+1] 	+= (yi - position[2*122+1]);

		if ((temp = sqrt(verlet_distance[2*122]*verlet_distance[2*122] + verlet_distance[2*122+1]*verlet_distance[2*122+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*122] 	-= floor(position[2*122]/L_x)*L_x;

		xi = position[2*123];
		yi = position[2*123+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*123]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*123+1] + weigh_brown_B * g2;

		position[2*123]		+= dx;
		position[2*123+1] 	+= dy;

		displacement[2*123]	+= dx;
		displacement[2*123+1] += dy;

		verlet_distance[2*123] 	+= (xi - position[2*123]);
		verlet_distance[2*123+1] 	+= (yi - position[2*123+1]);

		if ((temp = sqrt(verlet_distance[2*123]*verlet_distance[2*123] + verlet_distance[2*123+1]*verlet_distance[2*123+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*123] 	-= floor(position[2*123]/L_x)*L_x;

		xi = position[2*124];
		yi = position[2*124+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*124]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*124+1] + weigh_brown_B * g2;

		position[2*124]		+= dx;
		position[2*124+1] 	+= dy;

		displacement[2*124]	+= dx;
		displacement[2*124+1] += dy;

		verlet_distance[2*124] 	+= (xi - position[2*124]);
		verlet_distance[2*124+1] 	+= (yi - position[2*124+1]);

		if ((temp = sqrt(verlet_distance[2*124]*verlet_distance[2*124] + verlet_distance[2*124+1]*verlet_distance[2*124+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*124] 	-= floor(position[2*124]/L_x)*L_x;

		xi = position[2*125];
		yi = position[2*125+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*125]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*125+1] + weigh_brown_B * g2;

		position[2*125]		+= dx;
		position[2*125+1] 	+= dy;

		displacement[2*125]	+= dx;
		displacement[2*125+1] += dy;

		verlet_distance[2*125] 	+= (xi - position[2*125]);
		verlet_distance[2*125+1] 	+= (yi - position[2*125+1]);

		if ((temp = sqrt(verlet_distance[2*125]*verlet_distance[2*125] + verlet_distance[2*125+1]*verlet_distance[2*125+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*125] 	-= floor(position[2*125]/L_x)*L_x;

		xi = position[2*126];
		yi = position[2*126+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*126]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*126+1] + weigh_brown_B * g2;

		position[2*126]		+= dx;
		position[2*126+1] 	+= dy;

		displacement[2*126]	+= dx;
		displacement[2*126+1] += dy;

		verlet_distance[2*126] 	+= (xi - position[2*126]);
		verlet_distance[2*126+1] 	+= (yi - position[2*126+1]);

		if ((temp = sqrt(verlet_distance[2*126]*verlet_distance[2*126] + verlet_distance[2*126+1]*verlet_distance[2*126+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*126] 	-= floor(position[2*126]/L_x)*L_x;

		xi = position[2*127];
		yi = position[2*127+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*127]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*127+1] + weigh_brown_B * g2;

		position[2*127]		+= dx;
		position[2*127+1] 	+= dy;

		displacement[2*127]	+= dx;
		displacement[2*127+1] += dy;

		verlet_distance[2*127] 	+= (xi - position[2*127]);
		verlet_distance[2*127+1] 	+= (yi - position[2*127+1]);

		if ((temp = sqrt(verlet_distance[2*127]*verlet_distance[2*127] + verlet_distance[2*127+1]*verlet_distance[2*127+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*127] 	-= floor(position[2*127]/L_x)*L_x;

		xi = position[2*128];
		yi = position[2*128+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*128]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*128+1] + weigh_brown_B * g2;

		position[2*128]		+= dx;
		position[2*128+1] 	+= dy;

		displacement[2*128]	+= dx;
		displacement[2*128+1] += dy;

		verlet_distance[2*128] 	+= (xi - position[2*128]);
		verlet_distance[2*128+1] 	+= (yi - position[2*128+1]);

		if ((temp = sqrt(verlet_distance[2*128]*verlet_distance[2*128] + verlet_distance[2*128+1]*verlet_distance[2*128+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*128] 	-= floor(position[2*128]/L_x)*L_x;

		xi = position[2*129];
		yi = position[2*129+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*129]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*129+1] + weigh_brown_B * g2;

		position[2*129]		+= dx;
		position[2*129+1] 	+= dy;

		displacement[2*129]	+= dx;
		displacement[2*129+1] += dy;

		verlet_distance[2*129] 	+= (xi - position[2*129]);
		verlet_distance[2*129+1] 	+= (yi - position[2*129+1]);

		if ((temp = sqrt(verlet_distance[2*129]*verlet_distance[2*129] + verlet_distance[2*129+1]*verlet_distance[2*129+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*129] 	-= floor(position[2*129]/L_x)*L_x;

		xi = position[2*130];
		yi = position[2*130+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*130]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*130+1] + weigh_brown_B * g2;

		position[2*130]		+= dx;
		position[2*130+1] 	+= dy;

		displacement[2*130]	+= dx;
		displacement[2*130+1] += dy;

		verlet_distance[2*130] 	+= (xi - position[2*130]);
		verlet_distance[2*130+1] 	+= (yi - position[2*130+1]);

		if ((temp = sqrt(verlet_distance[2*130]*verlet_distance[2*130] + verlet_distance[2*130+1]*verlet_distance[2*130+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*130] 	-= floor(position[2*130]/L_x)*L_x;

		xi = position[2*131];
		yi = position[2*131+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*131]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*131+1] + weigh_brown_B * g2;

		position[2*131]		+= dx;
		position[2*131+1] 	+= dy;

		displacement[2*131]	+= dx;
		displacement[2*131+1] += dy;

		verlet_distance[2*131] 	+= (xi - position[2*131]);
		verlet_distance[2*131+1] 	+= (yi - position[2*131+1]);

		if ((temp = sqrt(verlet_distance[2*131]*verlet_distance[2*131] + verlet_distance[2*131+1]*verlet_distance[2*131+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*131] 	-= floor(position[2*131]/L_x)*L_x;

		xi = position[2*132];
		yi = position[2*132+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*132]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*132+1] + weigh_brown_B * g2;

		position[2*132]		+= dx;
		position[2*132+1] 	+= dy;

		displacement[2*132]	+= dx;
		displacement[2*132+1] += dy;

		verlet_distance[2*132] 	+= (xi - position[2*132]);
		verlet_distance[2*132+1] 	+= (yi - position[2*132+1]);

		if ((temp = sqrt(verlet_distance[2*132]*verlet_distance[2*132] + verlet_distance[2*132+1]*verlet_distance[2*132+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*132] 	-= floor(position[2*132]/L_x)*L_x;

		xi = position[2*133];
		yi = position[2*133+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*133]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*133+1] + weigh_brown_B * g2;

		position[2*133]		+= dx;
		position[2*133+1] 	+= dy;

		displacement[2*133]	+= dx;
		displacement[2*133+1] += dy;

		verlet_distance[2*133] 	+= (xi - position[2*133]);
		verlet_distance[2*133+1] 	+= (yi - position[2*133+1]);

		if ((temp = sqrt(verlet_distance[2*133]*verlet_distance[2*133] + verlet_distance[2*133+1]*verlet_distance[2*133+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*133] 	-= floor(position[2*133]/L_x)*L_x;

		xi = position[2*134];
		yi = position[2*134+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*134]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*134+1] + weigh_brown_B * g2;

		position[2*134]		+= dx;
		position[2*134+1] 	+= dy;

		displacement[2*134]	+= dx;
		displacement[2*134+1] += dy;

		verlet_distance[2*134] 	+= (xi - position[2*134]);
		verlet_distance[2*134+1] 	+= (yi - position[2*134+1]);

		if ((temp = sqrt(verlet_distance[2*134]*verlet_distance[2*134] + verlet_distance[2*134+1]*verlet_distance[2*134+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*134] 	-= floor(position[2*134]/L_x)*L_x;

		xi = position[2*135];
		yi = position[2*135+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*135]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*135+1] + weigh_brown_B * g2;

		position[2*135]		+= dx;
		position[2*135+1] 	+= dy;

		displacement[2*135]	+= dx;
		displacement[2*135+1] += dy;

		verlet_distance[2*135] 	+= (xi - position[2*135]);
		verlet_distance[2*135+1] 	+= (yi - position[2*135+1]);

		if ((temp = sqrt(verlet_distance[2*135]*verlet_distance[2*135] + verlet_distance[2*135+1]*verlet_distance[2*135+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*135] 	-= floor(position[2*135]/L_x)*L_x;

		xi = position[2*136];
		yi = position[2*136+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*136]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*136+1] + weigh_brown_B * g2;

		position[2*136]		+= dx;
		position[2*136+1] 	+= dy;

		displacement[2*136]	+= dx;
		displacement[2*136+1] += dy;

		verlet_distance[2*136] 	+= (xi - position[2*136]);
		verlet_distance[2*136+1] 	+= (yi - position[2*136+1]);

		if ((temp = sqrt(verlet_distance[2*136]*verlet_distance[2*136] + verlet_distance[2*136+1]*verlet_distance[2*136+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*136] 	-= floor(position[2*136]/L_x)*L_x;

		xi = position[2*137];
		yi = position[2*137+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*137]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*137+1] + weigh_brown_B * g2;

		position[2*137]		+= dx;
		position[2*137+1] 	+= dy;

		displacement[2*137]	+= dx;
		displacement[2*137+1] += dy;

		verlet_distance[2*137] 	+= (xi - position[2*137]);
		verlet_distance[2*137+1] 	+= (yi - position[2*137+1]);

		if ((temp = sqrt(verlet_distance[2*137]*verlet_distance[2*137] + verlet_distance[2*137+1]*verlet_distance[2*137+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*137] 	-= floor(position[2*137]/L_x)*L_x;

		xi = position[2*138];
		yi = position[2*138+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*138]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*138+1] + weigh_brown_B * g2;

		position[2*138]		+= dx;
		position[2*138+1] 	+= dy;

		displacement[2*138]	+= dx;
		displacement[2*138+1] += dy;

		verlet_distance[2*138] 	+= (xi - position[2*138]);
		verlet_distance[2*138+1] 	+= (yi - position[2*138+1]);

		if ((temp = sqrt(verlet_distance[2*138]*verlet_distance[2*138] + verlet_distance[2*138+1]*verlet_distance[2*138+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*138] 	-= floor(position[2*138]/L_x)*L_x;

		xi = position[2*139];
		yi = position[2*139+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*139]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*139+1] + weigh_brown_B * g2;

		position[2*139]		+= dx;
		position[2*139+1] 	+= dy;

		displacement[2*139]	+= dx;
		displacement[2*139+1] += dy;

		verlet_distance[2*139] 	+= (xi - position[2*139]);
		verlet_distance[2*139+1] 	+= (yi - position[2*139+1]);

		if ((temp = sqrt(verlet_distance[2*139]*verlet_distance[2*139] + verlet_distance[2*139+1]*verlet_distance[2*139+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*139] 	-= floor(position[2*139]/L_x)*L_x;

		xi = position[2*140];
		yi = position[2*140+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*140]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*140+1] + weigh_brown_B * g2;

		position[2*140]		+= dx;
		position[2*140+1] 	+= dy;

		displacement[2*140]	+= dx;
		displacement[2*140+1] += dy;

		verlet_distance[2*140] 	+= (xi - position[2*140]);
		verlet_distance[2*140+1] 	+= (yi - position[2*140+1]);

		if ((temp = sqrt(verlet_distance[2*140]*verlet_distance[2*140] + verlet_distance[2*140+1]*verlet_distance[2*140+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*140] 	-= floor(position[2*140]/L_x)*L_x;

		xi = position[2*141];
		yi = position[2*141+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*141]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*141+1] + weigh_brown_B * g2;

		position[2*141]		+= dx;
		position[2*141+1] 	+= dy;

		displacement[2*141]	+= dx;
		displacement[2*141+1] += dy;

		verlet_distance[2*141] 	+= (xi - position[2*141]);
		verlet_distance[2*141+1] 	+= (yi - position[2*141+1]);

		if ((temp = sqrt(verlet_distance[2*141]*verlet_distance[2*141] + verlet_distance[2*141+1]*verlet_distance[2*141+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*141] 	-= floor(position[2*141]/L_x)*L_x;

		xi = position[2*142];
		yi = position[2*142+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*142]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*142+1] + weigh_brown_B * g2;

		position[2*142]		+= dx;
		position[2*142+1] 	+= dy;

		displacement[2*142]	+= dx;
		displacement[2*142+1] += dy;

		verlet_distance[2*142] 	+= (xi - position[2*142]);
		verlet_distance[2*142+1] 	+= (yi - position[2*142+1]);

		if ((temp = sqrt(verlet_distance[2*142]*verlet_distance[2*142] + verlet_distance[2*142+1]*verlet_distance[2*142+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*142] 	-= floor(position[2*142]/L_x)*L_x;

		xi = position[2*143];
		yi = position[2*143+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*143]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*143+1] + weigh_brown_B * g2;

		position[2*143]		+= dx;
		position[2*143+1] 	+= dy;

		displacement[2*143]	+= dx;
		displacement[2*143+1] += dy;

		verlet_distance[2*143] 	+= (xi - position[2*143]);
		verlet_distance[2*143+1] 	+= (yi - position[2*143+1]);

		if ((temp = sqrt(verlet_distance[2*143]*verlet_distance[2*143] + verlet_distance[2*143+1]*verlet_distance[2*143+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*143] 	-= floor(position[2*143]/L_x)*L_x;

		xi = position[2*144];
		yi = position[2*144+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*144]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*144+1] + weigh_brown_B * g2;

		position[2*144]		+= dx;
		position[2*144+1] 	+= dy;

		displacement[2*144]	+= dx;
		displacement[2*144+1] += dy;

		verlet_distance[2*144] 	+= (xi - position[2*144]);
		verlet_distance[2*144+1] 	+= (yi - position[2*144+1]);

		if ((temp = sqrt(verlet_distance[2*144]*verlet_distance[2*144] + verlet_distance[2*144+1]*verlet_distance[2*144+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*144] 	-= floor(position[2*144]/L_x)*L_x;

		xi = position[2*145];
		yi = position[2*145+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*145]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*145+1] + weigh_brown_B * g2;

		position[2*145]		+= dx;
		position[2*145+1] 	+= dy;

		displacement[2*145]	+= dx;
		displacement[2*145+1] += dy;

		verlet_distance[2*145] 	+= (xi - position[2*145]);
		verlet_distance[2*145+1] 	+= (yi - position[2*145+1]);

		if ((temp = sqrt(verlet_distance[2*145]*verlet_distance[2*145] + verlet_distance[2*145+1]*verlet_distance[2*145+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*145] 	-= floor(position[2*145]/L_x)*L_x;

		xi = position[2*146];
		yi = position[2*146+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*146]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*146+1] + weigh_brown_B * g2;

		position[2*146]		+= dx;
		position[2*146+1] 	+= dy;

		displacement[2*146]	+= dx;
		displacement[2*146+1] += dy;

		verlet_distance[2*146] 	+= (xi - position[2*146]);
		verlet_distance[2*146+1] 	+= (yi - position[2*146+1]);

		if ((temp = sqrt(verlet_distance[2*146]*verlet_distance[2*146] + verlet_distance[2*146+1]*verlet_distance[2*146+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*146] 	-= floor(position[2*146]/L_x)*L_x;

		xi = position[2*147];
		yi = position[2*147+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*147]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*147+1] + weigh_brown_B * g2;

		position[2*147]		+= dx;
		position[2*147+1] 	+= dy;

		displacement[2*147]	+= dx;
		displacement[2*147+1] += dy;

		verlet_distance[2*147] 	+= (xi - position[2*147]);
		verlet_distance[2*147+1] 	+= (yi - position[2*147+1]);

		if ((temp = sqrt(verlet_distance[2*147]*verlet_distance[2*147] + verlet_distance[2*147+1]*verlet_distance[2*147+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*147] 	-= floor(position[2*147]/L_x)*L_x;

		xi = position[2*148];
		yi = position[2*148+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*148]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*148+1] + weigh_brown_B * g2;

		position[2*148]		+= dx;
		position[2*148+1] 	+= dy;

		displacement[2*148]	+= dx;
		displacement[2*148+1] += dy;

		verlet_distance[2*148] 	+= (xi - position[2*148]);
		verlet_distance[2*148+1] 	+= (yi - position[2*148+1]);

		if ((temp = sqrt(verlet_distance[2*148]*verlet_distance[2*148] + verlet_distance[2*148+1]*verlet_distance[2*148+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*148] 	-= floor(position[2*148]/L_x)*L_x;

		xi = position[2*149];
		yi = position[2*149+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*149]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*149+1] + weigh_brown_B * g2;

		position[2*149]		+= dx;
		position[2*149+1] 	+= dy;

		displacement[2*149]	+= dx;
		displacement[2*149+1] += dy;

		verlet_distance[2*149] 	+= (xi - position[2*149]);
		verlet_distance[2*149+1] 	+= (yi - position[2*149+1]);

		if ((temp = sqrt(verlet_distance[2*149]*verlet_distance[2*149] + verlet_distance[2*149+1]*verlet_distance[2*149+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*149] 	-= floor(position[2*149]/L_x)*L_x;

		// Signal to main thread, that all threads have finished their iteration
		pthread_barrier_wait(&barrier_main_one);

		// thread is done with one iteration, it will now wait for a signal from the main thread to continue
		pthread_barrier_wait(&barrier_main_two);
	}

	return NULL;
}



static void *iteration_4 (int *no) {
	// Define necessary variables
	double xi, xj, yi, yj;
	double dx, dy;
	double r_squared;

	double temp_force;
	double temp;

	double u1, u2, g1, g2;

	int iterate;
	int j;

	double m_i_A = 1.0;
	double m_i_B = m;

	double D_kT_A = D_Brown_A/kT;
	double D_kT_B = D_Brown_B/kT;

	double shear_A = v_A * delta_t;

	double m_j;

	double dist_bottom;
	double dist_top;

	// cutoff for the walls force, this will lead to a smoother transition
	double dist_cutoff 			= 1.0;
	double prefactor			= 3.0;
	double force_wall_cutoff 	= 1.0*prefactor *(kappa/dist_cutoff + 1.0/(dist_cutoff*dist_cutoff))*exp(-kappa*dist_cutoff);

	while(cont == 1) {
		xi = position[2*150];
		yi = position[2*150+1];

		force[2*150]	 = 0;
		force[2*150+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*150+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*150+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*150+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*150+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(150+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*150+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*150] 		+= temp_force*dx;
			force[2*150+1]	+= temp_force*dy;
		}

		xi = position[2*151];
		yi = position[2*151+1];

		force[2*151]	 = 0;
		force[2*151+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*151+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*151+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*151+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*151+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(151+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*151+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*151] 		+= temp_force*dx;
			force[2*151+1]	+= temp_force*dy;
		}

		xi = position[2*152];
		yi = position[2*152+1];

		force[2*152]	 = 0;
		force[2*152+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*152+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*152+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*152+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*152+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(152+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*152+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*152] 		+= temp_force*dx;
			force[2*152+1]	+= temp_force*dy;
		}

		xi = position[2*153];
		yi = position[2*153+1];

		force[2*153]	 = 0;
		force[2*153+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*153+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*153+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*153+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*153+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(153+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*153+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*153] 		+= temp_force*dx;
			force[2*153+1]	+= temp_force*dy;
		}

		xi = position[2*154];
		yi = position[2*154+1];

		force[2*154]	 = 0;
		force[2*154+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*154+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*154+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*154+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*154+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(154+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*154+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*154] 		+= temp_force*dx;
			force[2*154+1]	+= temp_force*dy;
		}

		xi = position[2*155];
		yi = position[2*155+1];

		force[2*155]	 = 0;
		force[2*155+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*155+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*155+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*155+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*155+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(155+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*155+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*155] 		+= temp_force*dx;
			force[2*155+1]	+= temp_force*dy;
		}

		xi = position[2*156];
		yi = position[2*156+1];

		force[2*156]	 = 0;
		force[2*156+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*156+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*156+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*156+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*156+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(156+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*156+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*156] 		+= temp_force*dx;
			force[2*156+1]	+= temp_force*dy;
		}

		xi = position[2*157];
		yi = position[2*157+1];

		force[2*157]	 = 0;
		force[2*157+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*157+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*157+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*157+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*157+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(157+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*157+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*157] 		+= temp_force*dx;
			force[2*157+1]	+= temp_force*dy;
		}

		xi = position[2*158];
		yi = position[2*158+1];

		force[2*158]	 = 0;
		force[2*158+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*158+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*158+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*158+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*158+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(158+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*158+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*158] 		+= temp_force*dx;
			force[2*158+1]	+= temp_force*dy;
		}

		xi = position[2*159];
		yi = position[2*159+1];

		force[2*159]	 = 0;
		force[2*159+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*159+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*159+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*159+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*159+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(159+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*159+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*159] 		+= temp_force*dx;
			force[2*159+1]	+= temp_force*dy;
		}

		xi = position[2*160];
		yi = position[2*160+1];

		force[2*160]	 = 0;
		force[2*160+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*160+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*160+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*160+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*160+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(160+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*160+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*160] 		+= temp_force*dx;
			force[2*160+1]	+= temp_force*dy;
		}

		xi = position[2*161];
		yi = position[2*161+1];

		force[2*161]	 = 0;
		force[2*161+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*161+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*161+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*161+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*161+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(161+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*161+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*161] 		+= temp_force*dx;
			force[2*161+1]	+= temp_force*dy;
		}

		xi = position[2*162];
		yi = position[2*162+1];

		force[2*162]	 = 0;
		force[2*162+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*162+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*162+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*162+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*162+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(162+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*162+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*162] 		+= temp_force*dx;
			force[2*162+1]	+= temp_force*dy;
		}

		xi = position[2*163];
		yi = position[2*163+1];

		force[2*163]	 = 0;
		force[2*163+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*163+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*163+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*163+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*163+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(163+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*163+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*163] 		+= temp_force*dx;
			force[2*163+1]	+= temp_force*dy;
		}

		xi = position[2*164];
		yi = position[2*164+1];

		force[2*164]	 = 0;
		force[2*164+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*164+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*164+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*164+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*164+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(164+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*164+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*164] 		+= temp_force*dx;
			force[2*164+1]	+= temp_force*dy;
		}

		xi = position[2*165];
		yi = position[2*165+1];

		force[2*165]	 = 0;
		force[2*165+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*165+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*165+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*165+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*165+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(165+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*165+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*165] 		+= temp_force*dx;
			force[2*165+1]	+= temp_force*dy;
		}

		xi = position[2*166];
		yi = position[2*166+1];

		force[2*166]	 = 0;
		force[2*166+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*166+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*166+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*166+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*166+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(166+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*166+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*166] 		+= temp_force*dx;
			force[2*166+1]	+= temp_force*dy;
		}

		xi = position[2*167];
		yi = position[2*167+1];

		force[2*167]	 = 0;
		force[2*167+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*167+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*167+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*167+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*167+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(167+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*167+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*167] 		+= temp_force*dx;
			force[2*167+1]	+= temp_force*dy;
		}

		xi = position[2*168];
		yi = position[2*168+1];

		force[2*168]	 = 0;
		force[2*168+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*168+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*168+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*168+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*168+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(168+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*168+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*168] 		+= temp_force*dx;
			force[2*168+1]	+= temp_force*dy;
		}

		xi = position[2*169];
		yi = position[2*169+1];

		force[2*169]	 = 0;
		force[2*169+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*169+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*169+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*169+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*169+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(169+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*169+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*169] 		+= temp_force*dx;
			force[2*169+1]	+= temp_force*dy;
		}

		xi = position[2*170];
		yi = position[2*170+1];

		force[2*170]	 = 0;
		force[2*170+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*170+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*170+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*170+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*170+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(170+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*170+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*170] 		+= temp_force*dx;
			force[2*170+1]	+= temp_force*dy;
		}

		xi = position[2*171];
		yi = position[2*171+1];

		force[2*171]	 = 0;
		force[2*171+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*171+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*171+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*171+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*171+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(171+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*171+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*171] 		+= temp_force*dx;
			force[2*171+1]	+= temp_force*dy;
		}

		xi = position[2*172];
		yi = position[2*172+1];

		force[2*172]	 = 0;
		force[2*172+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*172+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*172+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*172+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*172+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(172+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*172+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*172] 		+= temp_force*dx;
			force[2*172+1]	+= temp_force*dy;
		}

		xi = position[2*173];
		yi = position[2*173+1];

		force[2*173]	 = 0;
		force[2*173+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*173+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*173+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*173+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*173+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(173+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*173+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*173] 		+= temp_force*dx;
			force[2*173+1]	+= temp_force*dy;
		}

		xi = position[2*174];
		yi = position[2*174+1];

		force[2*174]	 = 0;
		force[2*174+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*174+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*174+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*174+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*174+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(174+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*174+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*174] 		+= temp_force*dx;
			force[2*174+1]	+= temp_force*dy;
		}

		xi = position[2*175];
		yi = position[2*175+1];

		force[2*175]	 = 0;
		force[2*175+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*175+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*175+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*175+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*175+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(175+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*175+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*175] 		+= temp_force*dx;
			force[2*175+1]	+= temp_force*dy;
		}

		xi = position[2*176];
		yi = position[2*176+1];

		force[2*176]	 = 0;
		force[2*176+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*176+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*176+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*176+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*176+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(176+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*176+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*176] 		+= temp_force*dx;
			force[2*176+1]	+= temp_force*dy;
		}

		xi = position[2*177];
		yi = position[2*177+1];

		force[2*177]	 = 0;
		force[2*177+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*177+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*177+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*177+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*177+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(177+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*177+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*177] 		+= temp_force*dx;
			force[2*177+1]	+= temp_force*dy;
		}

		xi = position[2*178];
		yi = position[2*178+1];

		force[2*178]	 = 0;
		force[2*178+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*178+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*178+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*178+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*178+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(178+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*178+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*178] 		+= temp_force*dx;
			force[2*178+1]	+= temp_force*dy;
		}

		xi = position[2*179];
		yi = position[2*179+1];

		force[2*179]	 = 0;
		force[2*179+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*179+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*179+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*179+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*179+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(179+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*179+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*179] 		+= temp_force*dx;
			force[2*179+1]	+= temp_force*dy;
		}

		xi = position[2*180];
		yi = position[2*180+1];

		force[2*180]	 = 0;
		force[2*180+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*180+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*180+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*180+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*180+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(180+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*180+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*180] 		+= temp_force*dx;
			force[2*180+1]	+= temp_force*dy;
		}

		xi = position[2*181];
		yi = position[2*181+1];

		force[2*181]	 = 0;
		force[2*181+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*181+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*181+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*181+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*181+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(181+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*181+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*181] 		+= temp_force*dx;
			force[2*181+1]	+= temp_force*dy;
		}

		xi = position[2*182];
		yi = position[2*182+1];

		force[2*182]	 = 0;
		force[2*182+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*182+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*182+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*182+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*182+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(182+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*182+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*182] 		+= temp_force*dx;
			force[2*182+1]	+= temp_force*dy;
		}

		xi = position[2*183];
		yi = position[2*183+1];

		force[2*183]	 = 0;
		force[2*183+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*183+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*183+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*183+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*183+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(183+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*183+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*183] 		+= temp_force*dx;
			force[2*183+1]	+= temp_force*dy;
		}

		xi = position[2*184];
		yi = position[2*184+1];

		force[2*184]	 = 0;
		force[2*184+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*184+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*184+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*184+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*184+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(184+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*184+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*184] 		+= temp_force*dx;
			force[2*184+1]	+= temp_force*dy;
		}

		xi = position[2*185];
		yi = position[2*185+1];

		force[2*185]	 = 0;
		force[2*185+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*185+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*185+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*185+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*185+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(185+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*185+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*185] 		+= temp_force*dx;
			force[2*185+1]	+= temp_force*dy;
		}

		xi = position[2*186];
		yi = position[2*186+1];

		force[2*186]	 = 0;
		force[2*186+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*186+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*186+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*186+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*186+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(186+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*186+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*186] 		+= temp_force*dx;
			force[2*186+1]	+= temp_force*dy;
		}

		xi = position[2*187];
		yi = position[2*187+1];

		force[2*187]	 = 0;
		force[2*187+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*187+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*187+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*187+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*187+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(187+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*187+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*187] 		+= temp_force*dx;
			force[2*187+1]	+= temp_force*dy;
		}

		pthread_barrier_wait(&barrier_internal);

		xi = position[2*150];
		yi = position[2*150+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*150]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*150+1] + weigh_brown_B * g2;

		position[2*150]		+= dx;
		position[2*150+1] 	+= dy;

		displacement[2*150]	+= dx;
		displacement[2*150+1] += dy;

		verlet_distance[2*150] 	+= (xi - position[2*150]);
		verlet_distance[2*150+1] 	+= (yi - position[2*150+1]);

		if ((temp = sqrt(verlet_distance[2*150]*verlet_distance[2*150] + verlet_distance[2*150+1]*verlet_distance[2*150+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*150] 	-= floor(position[2*150]/L_x)*L_x;

		xi = position[2*151];
		yi = position[2*151+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*151]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*151+1] + weigh_brown_B * g2;

		position[2*151]		+= dx;
		position[2*151+1] 	+= dy;

		displacement[2*151]	+= dx;
		displacement[2*151+1] += dy;

		verlet_distance[2*151] 	+= (xi - position[2*151]);
		verlet_distance[2*151+1] 	+= (yi - position[2*151+1]);

		if ((temp = sqrt(verlet_distance[2*151]*verlet_distance[2*151] + verlet_distance[2*151+1]*verlet_distance[2*151+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*151] 	-= floor(position[2*151]/L_x)*L_x;

		xi = position[2*152];
		yi = position[2*152+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*152]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*152+1] + weigh_brown_B * g2;

		position[2*152]		+= dx;
		position[2*152+1] 	+= dy;

		displacement[2*152]	+= dx;
		displacement[2*152+1] += dy;

		verlet_distance[2*152] 	+= (xi - position[2*152]);
		verlet_distance[2*152+1] 	+= (yi - position[2*152+1]);

		if ((temp = sqrt(verlet_distance[2*152]*verlet_distance[2*152] + verlet_distance[2*152+1]*verlet_distance[2*152+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*152] 	-= floor(position[2*152]/L_x)*L_x;

		xi = position[2*153];
		yi = position[2*153+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*153]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*153+1] + weigh_brown_B * g2;

		position[2*153]		+= dx;
		position[2*153+1] 	+= dy;

		displacement[2*153]	+= dx;
		displacement[2*153+1] += dy;

		verlet_distance[2*153] 	+= (xi - position[2*153]);
		verlet_distance[2*153+1] 	+= (yi - position[2*153+1]);

		if ((temp = sqrt(verlet_distance[2*153]*verlet_distance[2*153] + verlet_distance[2*153+1]*verlet_distance[2*153+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*153] 	-= floor(position[2*153]/L_x)*L_x;

		xi = position[2*154];
		yi = position[2*154+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*154]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*154+1] + weigh_brown_B * g2;

		position[2*154]		+= dx;
		position[2*154+1] 	+= dy;

		displacement[2*154]	+= dx;
		displacement[2*154+1] += dy;

		verlet_distance[2*154] 	+= (xi - position[2*154]);
		verlet_distance[2*154+1] 	+= (yi - position[2*154+1]);

		if ((temp = sqrt(verlet_distance[2*154]*verlet_distance[2*154] + verlet_distance[2*154+1]*verlet_distance[2*154+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*154] 	-= floor(position[2*154]/L_x)*L_x;

		xi = position[2*155];
		yi = position[2*155+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*155]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*155+1] + weigh_brown_B * g2;

		position[2*155]		+= dx;
		position[2*155+1] 	+= dy;

		displacement[2*155]	+= dx;
		displacement[2*155+1] += dy;

		verlet_distance[2*155] 	+= (xi - position[2*155]);
		verlet_distance[2*155+1] 	+= (yi - position[2*155+1]);

		if ((temp = sqrt(verlet_distance[2*155]*verlet_distance[2*155] + verlet_distance[2*155+1]*verlet_distance[2*155+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*155] 	-= floor(position[2*155]/L_x)*L_x;

		xi = position[2*156];
		yi = position[2*156+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*156]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*156+1] + weigh_brown_B * g2;

		position[2*156]		+= dx;
		position[2*156+1] 	+= dy;

		displacement[2*156]	+= dx;
		displacement[2*156+1] += dy;

		verlet_distance[2*156] 	+= (xi - position[2*156]);
		verlet_distance[2*156+1] 	+= (yi - position[2*156+1]);

		if ((temp = sqrt(verlet_distance[2*156]*verlet_distance[2*156] + verlet_distance[2*156+1]*verlet_distance[2*156+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*156] 	-= floor(position[2*156]/L_x)*L_x;

		xi = position[2*157];
		yi = position[2*157+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*157]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*157+1] + weigh_brown_B * g2;

		position[2*157]		+= dx;
		position[2*157+1] 	+= dy;

		displacement[2*157]	+= dx;
		displacement[2*157+1] += dy;

		verlet_distance[2*157] 	+= (xi - position[2*157]);
		verlet_distance[2*157+1] 	+= (yi - position[2*157+1]);

		if ((temp = sqrt(verlet_distance[2*157]*verlet_distance[2*157] + verlet_distance[2*157+1]*verlet_distance[2*157+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*157] 	-= floor(position[2*157]/L_x)*L_x;

		xi = position[2*158];
		yi = position[2*158+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*158]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*158+1] + weigh_brown_B * g2;

		position[2*158]		+= dx;
		position[2*158+1] 	+= dy;

		displacement[2*158]	+= dx;
		displacement[2*158+1] += dy;

		verlet_distance[2*158] 	+= (xi - position[2*158]);
		verlet_distance[2*158+1] 	+= (yi - position[2*158+1]);

		if ((temp = sqrt(verlet_distance[2*158]*verlet_distance[2*158] + verlet_distance[2*158+1]*verlet_distance[2*158+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*158] 	-= floor(position[2*158]/L_x)*L_x;

		xi = position[2*159];
		yi = position[2*159+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*159]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*159+1] + weigh_brown_B * g2;

		position[2*159]		+= dx;
		position[2*159+1] 	+= dy;

		displacement[2*159]	+= dx;
		displacement[2*159+1] += dy;

		verlet_distance[2*159] 	+= (xi - position[2*159]);
		verlet_distance[2*159+1] 	+= (yi - position[2*159+1]);

		if ((temp = sqrt(verlet_distance[2*159]*verlet_distance[2*159] + verlet_distance[2*159+1]*verlet_distance[2*159+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*159] 	-= floor(position[2*159]/L_x)*L_x;

		xi = position[2*160];
		yi = position[2*160+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*160]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*160+1] + weigh_brown_B * g2;

		position[2*160]		+= dx;
		position[2*160+1] 	+= dy;

		displacement[2*160]	+= dx;
		displacement[2*160+1] += dy;

		verlet_distance[2*160] 	+= (xi - position[2*160]);
		verlet_distance[2*160+1] 	+= (yi - position[2*160+1]);

		if ((temp = sqrt(verlet_distance[2*160]*verlet_distance[2*160] + verlet_distance[2*160+1]*verlet_distance[2*160+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*160] 	-= floor(position[2*160]/L_x)*L_x;

		xi = position[2*161];
		yi = position[2*161+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*161]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*161+1] + weigh_brown_B * g2;

		position[2*161]		+= dx;
		position[2*161+1] 	+= dy;

		displacement[2*161]	+= dx;
		displacement[2*161+1] += dy;

		verlet_distance[2*161] 	+= (xi - position[2*161]);
		verlet_distance[2*161+1] 	+= (yi - position[2*161+1]);

		if ((temp = sqrt(verlet_distance[2*161]*verlet_distance[2*161] + verlet_distance[2*161+1]*verlet_distance[2*161+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*161] 	-= floor(position[2*161]/L_x)*L_x;

		xi = position[2*162];
		yi = position[2*162+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*162]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*162+1] + weigh_brown_B * g2;

		position[2*162]		+= dx;
		position[2*162+1] 	+= dy;

		displacement[2*162]	+= dx;
		displacement[2*162+1] += dy;

		verlet_distance[2*162] 	+= (xi - position[2*162]);
		verlet_distance[2*162+1] 	+= (yi - position[2*162+1]);

		if ((temp = sqrt(verlet_distance[2*162]*verlet_distance[2*162] + verlet_distance[2*162+1]*verlet_distance[2*162+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*162] 	-= floor(position[2*162]/L_x)*L_x;

		xi = position[2*163];
		yi = position[2*163+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*163]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*163+1] + weigh_brown_B * g2;

		position[2*163]		+= dx;
		position[2*163+1] 	+= dy;

		displacement[2*163]	+= dx;
		displacement[2*163+1] += dy;

		verlet_distance[2*163] 	+= (xi - position[2*163]);
		verlet_distance[2*163+1] 	+= (yi - position[2*163+1]);

		if ((temp = sqrt(verlet_distance[2*163]*verlet_distance[2*163] + verlet_distance[2*163+1]*verlet_distance[2*163+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*163] 	-= floor(position[2*163]/L_x)*L_x;

		xi = position[2*164];
		yi = position[2*164+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*164]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*164+1] + weigh_brown_B * g2;

		position[2*164]		+= dx;
		position[2*164+1] 	+= dy;

		displacement[2*164]	+= dx;
		displacement[2*164+1] += dy;

		verlet_distance[2*164] 	+= (xi - position[2*164]);
		verlet_distance[2*164+1] 	+= (yi - position[2*164+1]);

		if ((temp = sqrt(verlet_distance[2*164]*verlet_distance[2*164] + verlet_distance[2*164+1]*verlet_distance[2*164+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*164] 	-= floor(position[2*164]/L_x)*L_x;

		xi = position[2*165];
		yi = position[2*165+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*165]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*165+1] + weigh_brown_B * g2;

		position[2*165]		+= dx;
		position[2*165+1] 	+= dy;

		displacement[2*165]	+= dx;
		displacement[2*165+1] += dy;

		verlet_distance[2*165] 	+= (xi - position[2*165]);
		verlet_distance[2*165+1] 	+= (yi - position[2*165+1]);

		if ((temp = sqrt(verlet_distance[2*165]*verlet_distance[2*165] + verlet_distance[2*165+1]*verlet_distance[2*165+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*165] 	-= floor(position[2*165]/L_x)*L_x;

		xi = position[2*166];
		yi = position[2*166+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*166]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*166+1] + weigh_brown_B * g2;

		position[2*166]		+= dx;
		position[2*166+1] 	+= dy;

		displacement[2*166]	+= dx;
		displacement[2*166+1] += dy;

		verlet_distance[2*166] 	+= (xi - position[2*166]);
		verlet_distance[2*166+1] 	+= (yi - position[2*166+1]);

		if ((temp = sqrt(verlet_distance[2*166]*verlet_distance[2*166] + verlet_distance[2*166+1]*verlet_distance[2*166+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*166] 	-= floor(position[2*166]/L_x)*L_x;

		xi = position[2*167];
		yi = position[2*167+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*167]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*167+1] + weigh_brown_B * g2;

		position[2*167]		+= dx;
		position[2*167+1] 	+= dy;

		displacement[2*167]	+= dx;
		displacement[2*167+1] += dy;

		verlet_distance[2*167] 	+= (xi - position[2*167]);
		verlet_distance[2*167+1] 	+= (yi - position[2*167+1]);

		if ((temp = sqrt(verlet_distance[2*167]*verlet_distance[2*167] + verlet_distance[2*167+1]*verlet_distance[2*167+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*167] 	-= floor(position[2*167]/L_x)*L_x;

		xi = position[2*168];
		yi = position[2*168+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*168]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*168+1] + weigh_brown_B * g2;

		position[2*168]		+= dx;
		position[2*168+1] 	+= dy;

		displacement[2*168]	+= dx;
		displacement[2*168+1] += dy;

		verlet_distance[2*168] 	+= (xi - position[2*168]);
		verlet_distance[2*168+1] 	+= (yi - position[2*168+1]);

		if ((temp = sqrt(verlet_distance[2*168]*verlet_distance[2*168] + verlet_distance[2*168+1]*verlet_distance[2*168+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*168] 	-= floor(position[2*168]/L_x)*L_x;

		xi = position[2*169];
		yi = position[2*169+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*169]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*169+1] + weigh_brown_B * g2;

		position[2*169]		+= dx;
		position[2*169+1] 	+= dy;

		displacement[2*169]	+= dx;
		displacement[2*169+1] += dy;

		verlet_distance[2*169] 	+= (xi - position[2*169]);
		verlet_distance[2*169+1] 	+= (yi - position[2*169+1]);

		if ((temp = sqrt(verlet_distance[2*169]*verlet_distance[2*169] + verlet_distance[2*169+1]*verlet_distance[2*169+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*169] 	-= floor(position[2*169]/L_x)*L_x;

		xi = position[2*170];
		yi = position[2*170+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*170]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*170+1] + weigh_brown_B * g2;

		position[2*170]		+= dx;
		position[2*170+1] 	+= dy;

		displacement[2*170]	+= dx;
		displacement[2*170+1] += dy;

		verlet_distance[2*170] 	+= (xi - position[2*170]);
		verlet_distance[2*170+1] 	+= (yi - position[2*170+1]);

		if ((temp = sqrt(verlet_distance[2*170]*verlet_distance[2*170] + verlet_distance[2*170+1]*verlet_distance[2*170+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*170] 	-= floor(position[2*170]/L_x)*L_x;

		xi = position[2*171];
		yi = position[2*171+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*171]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*171+1] + weigh_brown_B * g2;

		position[2*171]		+= dx;
		position[2*171+1] 	+= dy;

		displacement[2*171]	+= dx;
		displacement[2*171+1] += dy;

		verlet_distance[2*171] 	+= (xi - position[2*171]);
		verlet_distance[2*171+1] 	+= (yi - position[2*171+1]);

		if ((temp = sqrt(verlet_distance[2*171]*verlet_distance[2*171] + verlet_distance[2*171+1]*verlet_distance[2*171+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*171] 	-= floor(position[2*171]/L_x)*L_x;

		xi = position[2*172];
		yi = position[2*172+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*172]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*172+1] + weigh_brown_B * g2;

		position[2*172]		+= dx;
		position[2*172+1] 	+= dy;

		displacement[2*172]	+= dx;
		displacement[2*172+1] += dy;

		verlet_distance[2*172] 	+= (xi - position[2*172]);
		verlet_distance[2*172+1] 	+= (yi - position[2*172+1]);

		if ((temp = sqrt(verlet_distance[2*172]*verlet_distance[2*172] + verlet_distance[2*172+1]*verlet_distance[2*172+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*172] 	-= floor(position[2*172]/L_x)*L_x;

		xi = position[2*173];
		yi = position[2*173+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*173]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*173+1] + weigh_brown_B * g2;

		position[2*173]		+= dx;
		position[2*173+1] 	+= dy;

		displacement[2*173]	+= dx;
		displacement[2*173+1] += dy;

		verlet_distance[2*173] 	+= (xi - position[2*173]);
		verlet_distance[2*173+1] 	+= (yi - position[2*173+1]);

		if ((temp = sqrt(verlet_distance[2*173]*verlet_distance[2*173] + verlet_distance[2*173+1]*verlet_distance[2*173+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*173] 	-= floor(position[2*173]/L_x)*L_x;

		xi = position[2*174];
		yi = position[2*174+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*174]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*174+1] + weigh_brown_B * g2;

		position[2*174]		+= dx;
		position[2*174+1] 	+= dy;

		displacement[2*174]	+= dx;
		displacement[2*174+1] += dy;

		verlet_distance[2*174] 	+= (xi - position[2*174]);
		verlet_distance[2*174+1] 	+= (yi - position[2*174+1]);

		if ((temp = sqrt(verlet_distance[2*174]*verlet_distance[2*174] + verlet_distance[2*174+1]*verlet_distance[2*174+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*174] 	-= floor(position[2*174]/L_x)*L_x;

		xi = position[2*175];
		yi = position[2*175+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*175]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*175+1] + weigh_brown_B * g2;

		position[2*175]		+= dx;
		position[2*175+1] 	+= dy;

		displacement[2*175]	+= dx;
		displacement[2*175+1] += dy;

		verlet_distance[2*175] 	+= (xi - position[2*175]);
		verlet_distance[2*175+1] 	+= (yi - position[2*175+1]);

		if ((temp = sqrt(verlet_distance[2*175]*verlet_distance[2*175] + verlet_distance[2*175+1]*verlet_distance[2*175+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*175] 	-= floor(position[2*175]/L_x)*L_x;

		xi = position[2*176];
		yi = position[2*176+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*176]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*176+1] + weigh_brown_B * g2;

		position[2*176]		+= dx;
		position[2*176+1] 	+= dy;

		displacement[2*176]	+= dx;
		displacement[2*176+1] += dy;

		verlet_distance[2*176] 	+= (xi - position[2*176]);
		verlet_distance[2*176+1] 	+= (yi - position[2*176+1]);

		if ((temp = sqrt(verlet_distance[2*176]*verlet_distance[2*176] + verlet_distance[2*176+1]*verlet_distance[2*176+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*176] 	-= floor(position[2*176]/L_x)*L_x;

		xi = position[2*177];
		yi = position[2*177+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*177]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*177+1] + weigh_brown_B * g2;

		position[2*177]		+= dx;
		position[2*177+1] 	+= dy;

		displacement[2*177]	+= dx;
		displacement[2*177+1] += dy;

		verlet_distance[2*177] 	+= (xi - position[2*177]);
		verlet_distance[2*177+1] 	+= (yi - position[2*177+1]);

		if ((temp = sqrt(verlet_distance[2*177]*verlet_distance[2*177] + verlet_distance[2*177+1]*verlet_distance[2*177+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*177] 	-= floor(position[2*177]/L_x)*L_x;

		xi = position[2*178];
		yi = position[2*178+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*178]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*178+1] + weigh_brown_B * g2;

		position[2*178]		+= dx;
		position[2*178+1] 	+= dy;

		displacement[2*178]	+= dx;
		displacement[2*178+1] += dy;

		verlet_distance[2*178] 	+= (xi - position[2*178]);
		verlet_distance[2*178+1] 	+= (yi - position[2*178+1]);

		if ((temp = sqrt(verlet_distance[2*178]*verlet_distance[2*178] + verlet_distance[2*178+1]*verlet_distance[2*178+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*178] 	-= floor(position[2*178]/L_x)*L_x;

		xi = position[2*179];
		yi = position[2*179+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*179]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*179+1] + weigh_brown_B * g2;

		position[2*179]		+= dx;
		position[2*179+1] 	+= dy;

		displacement[2*179]	+= dx;
		displacement[2*179+1] += dy;

		verlet_distance[2*179] 	+= (xi - position[2*179]);
		verlet_distance[2*179+1] 	+= (yi - position[2*179+1]);

		if ((temp = sqrt(verlet_distance[2*179]*verlet_distance[2*179] + verlet_distance[2*179+1]*verlet_distance[2*179+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*179] 	-= floor(position[2*179]/L_x)*L_x;

		xi = position[2*180];
		yi = position[2*180+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*180]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*180+1] + weigh_brown_B * g2;

		position[2*180]		+= dx;
		position[2*180+1] 	+= dy;

		displacement[2*180]	+= dx;
		displacement[2*180+1] += dy;

		verlet_distance[2*180] 	+= (xi - position[2*180]);
		verlet_distance[2*180+1] 	+= (yi - position[2*180+1]);

		if ((temp = sqrt(verlet_distance[2*180]*verlet_distance[2*180] + verlet_distance[2*180+1]*verlet_distance[2*180+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*180] 	-= floor(position[2*180]/L_x)*L_x;

		xi = position[2*181];
		yi = position[2*181+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*181]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*181+1] + weigh_brown_B * g2;

		position[2*181]		+= dx;
		position[2*181+1] 	+= dy;

		displacement[2*181]	+= dx;
		displacement[2*181+1] += dy;

		verlet_distance[2*181] 	+= (xi - position[2*181]);
		verlet_distance[2*181+1] 	+= (yi - position[2*181+1]);

		if ((temp = sqrt(verlet_distance[2*181]*verlet_distance[2*181] + verlet_distance[2*181+1]*verlet_distance[2*181+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*181] 	-= floor(position[2*181]/L_x)*L_x;

		xi = position[2*182];
		yi = position[2*182+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*182]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*182+1] + weigh_brown_B * g2;

		position[2*182]		+= dx;
		position[2*182+1] 	+= dy;

		displacement[2*182]	+= dx;
		displacement[2*182+1] += dy;

		verlet_distance[2*182] 	+= (xi - position[2*182]);
		verlet_distance[2*182+1] 	+= (yi - position[2*182+1]);

		if ((temp = sqrt(verlet_distance[2*182]*verlet_distance[2*182] + verlet_distance[2*182+1]*verlet_distance[2*182+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*182] 	-= floor(position[2*182]/L_x)*L_x;

		xi = position[2*183];
		yi = position[2*183+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*183]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*183+1] + weigh_brown_B * g2;

		position[2*183]		+= dx;
		position[2*183+1] 	+= dy;

		displacement[2*183]	+= dx;
		displacement[2*183+1] += dy;

		verlet_distance[2*183] 	+= (xi - position[2*183]);
		verlet_distance[2*183+1] 	+= (yi - position[2*183+1]);

		if ((temp = sqrt(verlet_distance[2*183]*verlet_distance[2*183] + verlet_distance[2*183+1]*verlet_distance[2*183+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*183] 	-= floor(position[2*183]/L_x)*L_x;

		xi = position[2*184];
		yi = position[2*184+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*184]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*184+1] + weigh_brown_B * g2;

		position[2*184]		+= dx;
		position[2*184+1] 	+= dy;

		displacement[2*184]	+= dx;
		displacement[2*184+1] += dy;

		verlet_distance[2*184] 	+= (xi - position[2*184]);
		verlet_distance[2*184+1] 	+= (yi - position[2*184+1]);

		if ((temp = sqrt(verlet_distance[2*184]*verlet_distance[2*184] + verlet_distance[2*184+1]*verlet_distance[2*184+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*184] 	-= floor(position[2*184]/L_x)*L_x;

		xi = position[2*185];
		yi = position[2*185+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*185]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*185+1] + weigh_brown_B * g2;

		position[2*185]		+= dx;
		position[2*185+1] 	+= dy;

		displacement[2*185]	+= dx;
		displacement[2*185+1] += dy;

		verlet_distance[2*185] 	+= (xi - position[2*185]);
		verlet_distance[2*185+1] 	+= (yi - position[2*185+1]);

		if ((temp = sqrt(verlet_distance[2*185]*verlet_distance[2*185] + verlet_distance[2*185+1]*verlet_distance[2*185+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*185] 	-= floor(position[2*185]/L_x)*L_x;

		xi = position[2*186];
		yi = position[2*186+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*186]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*186+1] + weigh_brown_B * g2;

		position[2*186]		+= dx;
		position[2*186+1] 	+= dy;

		displacement[2*186]	+= dx;
		displacement[2*186+1] += dy;

		verlet_distance[2*186] 	+= (xi - position[2*186]);
		verlet_distance[2*186+1] 	+= (yi - position[2*186+1]);

		if ((temp = sqrt(verlet_distance[2*186]*verlet_distance[2*186] + verlet_distance[2*186+1]*verlet_distance[2*186+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*186] 	-= floor(position[2*186]/L_x)*L_x;

		xi = position[2*187];
		yi = position[2*187+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*187]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*187+1] + weigh_brown_B * g2;

		position[2*187]		+= dx;
		position[2*187+1] 	+= dy;

		displacement[2*187]	+= dx;
		displacement[2*187+1] += dy;

		verlet_distance[2*187] 	+= (xi - position[2*187]);
		verlet_distance[2*187+1] 	+= (yi - position[2*187+1]);

		if ((temp = sqrt(verlet_distance[2*187]*verlet_distance[2*187] + verlet_distance[2*187+1]*verlet_distance[2*187+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*187] 	-= floor(position[2*187]/L_x)*L_x;

		// Signal to main thread, that all threads have finished their iteration
		pthread_barrier_wait(&barrier_main_one);

		// thread is done with one iteration, it will now wait for a signal from the main thread to continue
		pthread_barrier_wait(&barrier_main_two);
	}

	return NULL;
}



static void *iteration_5 (int *no) {
	// Define necessary variables
	double xi, xj, yi, yj;
	double dx, dy;
	double r_squared;

	double temp_force;
	double temp;

	double u1, u2, g1, g2;

	int iterate;
	int j;

	double m_i_A = 1.0;
	double m_i_B = m;

	double D_kT_A = D_Brown_A/kT;
	double D_kT_B = D_Brown_B/kT;

	double shear_A = v_A * delta_t;

	double m_j;

	double dist_bottom;
	double dist_top;

	// cutoff for the walls force, this will lead to a smoother transition
	double dist_cutoff 			= 1.0;
	double prefactor			= 3.0;
	double force_wall_cutoff 	= 1.0*prefactor *(kappa/dist_cutoff + 1.0/(dist_cutoff*dist_cutoff))*exp(-kappa*dist_cutoff);

	while(cont == 1) {
		xi = position[2*188];
		yi = position[2*188+1];

		force[2*188]	 = 0;
		force[2*188+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*188+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*188+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*188+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*188+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(188+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*188+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*188] 		+= temp_force*dx;
			force[2*188+1]	+= temp_force*dy;
		}

		xi = position[2*189];
		yi = position[2*189+1];

		force[2*189]	 = 0;
		force[2*189+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*189+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*189+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*189+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*189+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(189+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*189+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*189] 		+= temp_force*dx;
			force[2*189+1]	+= temp_force*dy;
		}

		xi = position[2*190];
		yi = position[2*190+1];

		force[2*190]	 = 0;
		force[2*190+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*190+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*190+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*190+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*190+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(190+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*190+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*190] 		+= temp_force*dx;
			force[2*190+1]	+= temp_force*dy;
		}

		xi = position[2*191];
		yi = position[2*191+1];

		force[2*191]	 = 0;
		force[2*191+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*191+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*191+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*191+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*191+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(191+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*191+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*191] 		+= temp_force*dx;
			force[2*191+1]	+= temp_force*dy;
		}

		xi = position[2*192];
		yi = position[2*192+1];

		force[2*192]	 = 0;
		force[2*192+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*192+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*192+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*192+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*192+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(192+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*192+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*192] 		+= temp_force*dx;
			force[2*192+1]	+= temp_force*dy;
		}

		xi = position[2*193];
		yi = position[2*193+1];

		force[2*193]	 = 0;
		force[2*193+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*193+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*193+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*193+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*193+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(193+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*193+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*193] 		+= temp_force*dx;
			force[2*193+1]	+= temp_force*dy;
		}

		xi = position[2*194];
		yi = position[2*194+1];

		force[2*194]	 = 0;
		force[2*194+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*194+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*194+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*194+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*194+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(194+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*194+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*194] 		+= temp_force*dx;
			force[2*194+1]	+= temp_force*dy;
		}

		xi = position[2*195];
		yi = position[2*195+1];

		force[2*195]	 = 0;
		force[2*195+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*195+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*195+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*195+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*195+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(195+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*195+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*195] 		+= temp_force*dx;
			force[2*195+1]	+= temp_force*dy;
		}

		xi = position[2*196];
		yi = position[2*196+1];

		force[2*196]	 = 0;
		force[2*196+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*196+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*196+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*196+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*196+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(196+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*196+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*196] 		+= temp_force*dx;
			force[2*196+1]	+= temp_force*dy;
		}

		xi = position[2*197];
		yi = position[2*197+1];

		force[2*197]	 = 0;
		force[2*197+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*197+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*197+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*197+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*197+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(197+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*197+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*197] 		+= temp_force*dx;
			force[2*197+1]	+= temp_force*dy;
		}

		xi = position[2*198];
		yi = position[2*198+1];

		force[2*198]	 = 0;
		force[2*198+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*198+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*198+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*198+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*198+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(198+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*198+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*198] 		+= temp_force*dx;
			force[2*198+1]	+= temp_force*dy;
		}

		xi = position[2*199];
		yi = position[2*199+1];

		force[2*199]	 = 0;
		force[2*199+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*199+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*199+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*199+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*199+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(199+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*199+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*199] 		+= temp_force*dx;
			force[2*199+1]	+= temp_force*dy;
		}

		xi = position[2*200];
		yi = position[2*200+1];

		force[2*200]	 = 0;
		force[2*200+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*200+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*200+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*200+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*200+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(200+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*200+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*200] 		+= temp_force*dx;
			force[2*200+1]	+= temp_force*dy;
		}

		xi = position[2*201];
		yi = position[2*201+1];

		force[2*201]	 = 0;
		force[2*201+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*201+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*201+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*201+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*201+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(201+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*201+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*201] 		+= temp_force*dx;
			force[2*201+1]	+= temp_force*dy;
		}

		xi = position[2*202];
		yi = position[2*202+1];

		force[2*202]	 = 0;
		force[2*202+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*202+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*202+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*202+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*202+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(202+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*202+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*202] 		+= temp_force*dx;
			force[2*202+1]	+= temp_force*dy;
		}

		xi = position[2*203];
		yi = position[2*203+1];

		force[2*203]	 = 0;
		force[2*203+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*203+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*203+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*203+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*203+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(203+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*203+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*203] 		+= temp_force*dx;
			force[2*203+1]	+= temp_force*dy;
		}

		xi = position[2*204];
		yi = position[2*204+1];

		force[2*204]	 = 0;
		force[2*204+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*204+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*204+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*204+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*204+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(204+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*204+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*204] 		+= temp_force*dx;
			force[2*204+1]	+= temp_force*dy;
		}

		xi = position[2*205];
		yi = position[2*205+1];

		force[2*205]	 = 0;
		force[2*205+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*205+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*205+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*205+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*205+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(205+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*205+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*205] 		+= temp_force*dx;
			force[2*205+1]	+= temp_force*dy;
		}

		xi = position[2*206];
		yi = position[2*206+1];

		force[2*206]	 = 0;
		force[2*206+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*206+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*206+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*206+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*206+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(206+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*206+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*206] 		+= temp_force*dx;
			force[2*206+1]	+= temp_force*dy;
		}

		xi = position[2*207];
		yi = position[2*207+1];

		force[2*207]	 = 0;
		force[2*207+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*207+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*207+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*207+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*207+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(207+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*207+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*207] 		+= temp_force*dx;
			force[2*207+1]	+= temp_force*dy;
		}

		xi = position[2*208];
		yi = position[2*208+1];

		force[2*208]	 = 0;
		force[2*208+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*208+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*208+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*208+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*208+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(208+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*208+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*208] 		+= temp_force*dx;
			force[2*208+1]	+= temp_force*dy;
		}

		xi = position[2*209];
		yi = position[2*209+1];

		force[2*209]	 = 0;
		force[2*209+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*209+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*209+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*209+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*209+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(209+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*209+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*209] 		+= temp_force*dx;
			force[2*209+1]	+= temp_force*dy;
		}

		xi = position[2*210];
		yi = position[2*210+1];

		force[2*210]	 = 0;
		force[2*210+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*210+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*210+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*210+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*210+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(210+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*210+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*210] 		+= temp_force*dx;
			force[2*210+1]	+= temp_force*dy;
		}

		xi = position[2*211];
		yi = position[2*211+1];

		force[2*211]	 = 0;
		force[2*211+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*211+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*211+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*211+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*211+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(211+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*211+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*211] 		+= temp_force*dx;
			force[2*211+1]	+= temp_force*dy;
		}

		xi = position[2*212];
		yi = position[2*212+1];

		force[2*212]	 = 0;
		force[2*212+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*212+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*212+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*212+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*212+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(212+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*212+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*212] 		+= temp_force*dx;
			force[2*212+1]	+= temp_force*dy;
		}

		xi = position[2*213];
		yi = position[2*213+1];

		force[2*213]	 = 0;
		force[2*213+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*213+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*213+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*213+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*213+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(213+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*213+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*213] 		+= temp_force*dx;
			force[2*213+1]	+= temp_force*dy;
		}

		xi = position[2*214];
		yi = position[2*214+1];

		force[2*214]	 = 0;
		force[2*214+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*214+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*214+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*214+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*214+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(214+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*214+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*214] 		+= temp_force*dx;
			force[2*214+1]	+= temp_force*dy;
		}

		xi = position[2*215];
		yi = position[2*215+1];

		force[2*215]	 = 0;
		force[2*215+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*215+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*215+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*215+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*215+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(215+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*215+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*215] 		+= temp_force*dx;
			force[2*215+1]	+= temp_force*dy;
		}

		xi = position[2*216];
		yi = position[2*216+1];

		force[2*216]	 = 0;
		force[2*216+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*216+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*216+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*216+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*216+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(216+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*216+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*216] 		+= temp_force*dx;
			force[2*216+1]	+= temp_force*dy;
		}

		xi = position[2*217];
		yi = position[2*217+1];

		force[2*217]	 = 0;
		force[2*217+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*217+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*217+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*217+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*217+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(217+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*217+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*217] 		+= temp_force*dx;
			force[2*217+1]	+= temp_force*dy;
		}

		xi = position[2*218];
		yi = position[2*218+1];

		force[2*218]	 = 0;
		force[2*218+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*218+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*218+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*218+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*218+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(218+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*218+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*218] 		+= temp_force*dx;
			force[2*218+1]	+= temp_force*dy;
		}

		xi = position[2*219];
		yi = position[2*219+1];

		force[2*219]	 = 0;
		force[2*219+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*219+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*219+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*219+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*219+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(219+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*219+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*219] 		+= temp_force*dx;
			force[2*219+1]	+= temp_force*dy;
		}

		xi = position[2*220];
		yi = position[2*220+1];

		force[2*220]	 = 0;
		force[2*220+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*220+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*220+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*220+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*220+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(220+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*220+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*220] 		+= temp_force*dx;
			force[2*220+1]	+= temp_force*dy;
		}

		xi = position[2*221];
		yi = position[2*221+1];

		force[2*221]	 = 0;
		force[2*221+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*221+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*221+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*221+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*221+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(221+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*221+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*221] 		+= temp_force*dx;
			force[2*221+1]	+= temp_force*dy;
		}

		xi = position[2*222];
		yi = position[2*222+1];

		force[2*222]	 = 0;
		force[2*222+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*222+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*222+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*222+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*222+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(222+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*222+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*222] 		+= temp_force*dx;
			force[2*222+1]	+= temp_force*dy;
		}

		xi = position[2*223];
		yi = position[2*223+1];

		force[2*223]	 = 0;
		force[2*223+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*223+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*223+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*223+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*223+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(223+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*223+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*223] 		+= temp_force*dx;
			force[2*223+1]	+= temp_force*dy;
		}

		xi = position[2*224];
		yi = position[2*224+1];

		force[2*224]	 = 0;
		force[2*224+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*224+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*224+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*224+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*224+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(224+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*224+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*224] 		+= temp_force*dx;
			force[2*224+1]	+= temp_force*dy;
		}

		pthread_barrier_wait(&barrier_internal);

		xi = position[2*188];
		yi = position[2*188+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*188]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*188+1] + weigh_brown_B * g2;

		position[2*188]		+= dx;
		position[2*188+1] 	+= dy;

		displacement[2*188]	+= dx;
		displacement[2*188+1] += dy;

		verlet_distance[2*188] 	+= (xi - position[2*188]);
		verlet_distance[2*188+1] 	+= (yi - position[2*188+1]);

		if ((temp = sqrt(verlet_distance[2*188]*verlet_distance[2*188] + verlet_distance[2*188+1]*verlet_distance[2*188+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*188] 	-= floor(position[2*188]/L_x)*L_x;

		xi = position[2*189];
		yi = position[2*189+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*189]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*189+1] + weigh_brown_B * g2;

		position[2*189]		+= dx;
		position[2*189+1] 	+= dy;

		displacement[2*189]	+= dx;
		displacement[2*189+1] += dy;

		verlet_distance[2*189] 	+= (xi - position[2*189]);
		verlet_distance[2*189+1] 	+= (yi - position[2*189+1]);

		if ((temp = sqrt(verlet_distance[2*189]*verlet_distance[2*189] + verlet_distance[2*189+1]*verlet_distance[2*189+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*189] 	-= floor(position[2*189]/L_x)*L_x;

		xi = position[2*190];
		yi = position[2*190+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*190]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*190+1] + weigh_brown_B * g2;

		position[2*190]		+= dx;
		position[2*190+1] 	+= dy;

		displacement[2*190]	+= dx;
		displacement[2*190+1] += dy;

		verlet_distance[2*190] 	+= (xi - position[2*190]);
		verlet_distance[2*190+1] 	+= (yi - position[2*190+1]);

		if ((temp = sqrt(verlet_distance[2*190]*verlet_distance[2*190] + verlet_distance[2*190+1]*verlet_distance[2*190+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*190] 	-= floor(position[2*190]/L_x)*L_x;

		xi = position[2*191];
		yi = position[2*191+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*191]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*191+1] + weigh_brown_B * g2;

		position[2*191]		+= dx;
		position[2*191+1] 	+= dy;

		displacement[2*191]	+= dx;
		displacement[2*191+1] += dy;

		verlet_distance[2*191] 	+= (xi - position[2*191]);
		verlet_distance[2*191+1] 	+= (yi - position[2*191+1]);

		if ((temp = sqrt(verlet_distance[2*191]*verlet_distance[2*191] + verlet_distance[2*191+1]*verlet_distance[2*191+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*191] 	-= floor(position[2*191]/L_x)*L_x;

		xi = position[2*192];
		yi = position[2*192+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*192]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*192+1] + weigh_brown_B * g2;

		position[2*192]		+= dx;
		position[2*192+1] 	+= dy;

		displacement[2*192]	+= dx;
		displacement[2*192+1] += dy;

		verlet_distance[2*192] 	+= (xi - position[2*192]);
		verlet_distance[2*192+1] 	+= (yi - position[2*192+1]);

		if ((temp = sqrt(verlet_distance[2*192]*verlet_distance[2*192] + verlet_distance[2*192+1]*verlet_distance[2*192+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*192] 	-= floor(position[2*192]/L_x)*L_x;

		xi = position[2*193];
		yi = position[2*193+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*193]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*193+1] + weigh_brown_B * g2;

		position[2*193]		+= dx;
		position[2*193+1] 	+= dy;

		displacement[2*193]	+= dx;
		displacement[2*193+1] += dy;

		verlet_distance[2*193] 	+= (xi - position[2*193]);
		verlet_distance[2*193+1] 	+= (yi - position[2*193+1]);

		if ((temp = sqrt(verlet_distance[2*193]*verlet_distance[2*193] + verlet_distance[2*193+1]*verlet_distance[2*193+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*193] 	-= floor(position[2*193]/L_x)*L_x;

		xi = position[2*194];
		yi = position[2*194+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*194]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*194+1] + weigh_brown_B * g2;

		position[2*194]		+= dx;
		position[2*194+1] 	+= dy;

		displacement[2*194]	+= dx;
		displacement[2*194+1] += dy;

		verlet_distance[2*194] 	+= (xi - position[2*194]);
		verlet_distance[2*194+1] 	+= (yi - position[2*194+1]);

		if ((temp = sqrt(verlet_distance[2*194]*verlet_distance[2*194] + verlet_distance[2*194+1]*verlet_distance[2*194+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*194] 	-= floor(position[2*194]/L_x)*L_x;

		xi = position[2*195];
		yi = position[2*195+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*195]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*195+1] + weigh_brown_B * g2;

		position[2*195]		+= dx;
		position[2*195+1] 	+= dy;

		displacement[2*195]	+= dx;
		displacement[2*195+1] += dy;

		verlet_distance[2*195] 	+= (xi - position[2*195]);
		verlet_distance[2*195+1] 	+= (yi - position[2*195+1]);

		if ((temp = sqrt(verlet_distance[2*195]*verlet_distance[2*195] + verlet_distance[2*195+1]*verlet_distance[2*195+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*195] 	-= floor(position[2*195]/L_x)*L_x;

		xi = position[2*196];
		yi = position[2*196+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*196]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*196+1] + weigh_brown_B * g2;

		position[2*196]		+= dx;
		position[2*196+1] 	+= dy;

		displacement[2*196]	+= dx;
		displacement[2*196+1] += dy;

		verlet_distance[2*196] 	+= (xi - position[2*196]);
		verlet_distance[2*196+1] 	+= (yi - position[2*196+1]);

		if ((temp = sqrt(verlet_distance[2*196]*verlet_distance[2*196] + verlet_distance[2*196+1]*verlet_distance[2*196+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*196] 	-= floor(position[2*196]/L_x)*L_x;

		xi = position[2*197];
		yi = position[2*197+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*197]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*197+1] + weigh_brown_B * g2;

		position[2*197]		+= dx;
		position[2*197+1] 	+= dy;

		displacement[2*197]	+= dx;
		displacement[2*197+1] += dy;

		verlet_distance[2*197] 	+= (xi - position[2*197]);
		verlet_distance[2*197+1] 	+= (yi - position[2*197+1]);

		if ((temp = sqrt(verlet_distance[2*197]*verlet_distance[2*197] + verlet_distance[2*197+1]*verlet_distance[2*197+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*197] 	-= floor(position[2*197]/L_x)*L_x;

		xi = position[2*198];
		yi = position[2*198+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*198]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*198+1] + weigh_brown_B * g2;

		position[2*198]		+= dx;
		position[2*198+1] 	+= dy;

		displacement[2*198]	+= dx;
		displacement[2*198+1] += dy;

		verlet_distance[2*198] 	+= (xi - position[2*198]);
		verlet_distance[2*198+1] 	+= (yi - position[2*198+1]);

		if ((temp = sqrt(verlet_distance[2*198]*verlet_distance[2*198] + verlet_distance[2*198+1]*verlet_distance[2*198+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*198] 	-= floor(position[2*198]/L_x)*L_x;

		xi = position[2*199];
		yi = position[2*199+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*199]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*199+1] + weigh_brown_B * g2;

		position[2*199]		+= dx;
		position[2*199+1] 	+= dy;

		displacement[2*199]	+= dx;
		displacement[2*199+1] += dy;

		verlet_distance[2*199] 	+= (xi - position[2*199]);
		verlet_distance[2*199+1] 	+= (yi - position[2*199+1]);

		if ((temp = sqrt(verlet_distance[2*199]*verlet_distance[2*199] + verlet_distance[2*199+1]*verlet_distance[2*199+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*199] 	-= floor(position[2*199]/L_x)*L_x;

		xi = position[2*200];
		yi = position[2*200+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*200]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*200+1] + weigh_brown_B * g2;

		position[2*200]		+= dx;
		position[2*200+1] 	+= dy;

		displacement[2*200]	+= dx;
		displacement[2*200+1] += dy;

		verlet_distance[2*200] 	+= (xi - position[2*200]);
		verlet_distance[2*200+1] 	+= (yi - position[2*200+1]);

		if ((temp = sqrt(verlet_distance[2*200]*verlet_distance[2*200] + verlet_distance[2*200+1]*verlet_distance[2*200+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*200] 	-= floor(position[2*200]/L_x)*L_x;

		xi = position[2*201];
		yi = position[2*201+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*201]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*201+1] + weigh_brown_B * g2;

		position[2*201]		+= dx;
		position[2*201+1] 	+= dy;

		displacement[2*201]	+= dx;
		displacement[2*201+1] += dy;

		verlet_distance[2*201] 	+= (xi - position[2*201]);
		verlet_distance[2*201+1] 	+= (yi - position[2*201+1]);

		if ((temp = sqrt(verlet_distance[2*201]*verlet_distance[2*201] + verlet_distance[2*201+1]*verlet_distance[2*201+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*201] 	-= floor(position[2*201]/L_x)*L_x;

		xi = position[2*202];
		yi = position[2*202+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*202]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*202+1] + weigh_brown_B * g2;

		position[2*202]		+= dx;
		position[2*202+1] 	+= dy;

		displacement[2*202]	+= dx;
		displacement[2*202+1] += dy;

		verlet_distance[2*202] 	+= (xi - position[2*202]);
		verlet_distance[2*202+1] 	+= (yi - position[2*202+1]);

		if ((temp = sqrt(verlet_distance[2*202]*verlet_distance[2*202] + verlet_distance[2*202+1]*verlet_distance[2*202+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*202] 	-= floor(position[2*202]/L_x)*L_x;

		xi = position[2*203];
		yi = position[2*203+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*203]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*203+1] + weigh_brown_B * g2;

		position[2*203]		+= dx;
		position[2*203+1] 	+= dy;

		displacement[2*203]	+= dx;
		displacement[2*203+1] += dy;

		verlet_distance[2*203] 	+= (xi - position[2*203]);
		verlet_distance[2*203+1] 	+= (yi - position[2*203+1]);

		if ((temp = sqrt(verlet_distance[2*203]*verlet_distance[2*203] + verlet_distance[2*203+1]*verlet_distance[2*203+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*203] 	-= floor(position[2*203]/L_x)*L_x;

		xi = position[2*204];
		yi = position[2*204+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*204]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*204+1] + weigh_brown_B * g2;

		position[2*204]		+= dx;
		position[2*204+1] 	+= dy;

		displacement[2*204]	+= dx;
		displacement[2*204+1] += dy;

		verlet_distance[2*204] 	+= (xi - position[2*204]);
		verlet_distance[2*204+1] 	+= (yi - position[2*204+1]);

		if ((temp = sqrt(verlet_distance[2*204]*verlet_distance[2*204] + verlet_distance[2*204+1]*verlet_distance[2*204+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*204] 	-= floor(position[2*204]/L_x)*L_x;

		xi = position[2*205];
		yi = position[2*205+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*205]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*205+1] + weigh_brown_B * g2;

		position[2*205]		+= dx;
		position[2*205+1] 	+= dy;

		displacement[2*205]	+= dx;
		displacement[2*205+1] += dy;

		verlet_distance[2*205] 	+= (xi - position[2*205]);
		verlet_distance[2*205+1] 	+= (yi - position[2*205+1]);

		if ((temp = sqrt(verlet_distance[2*205]*verlet_distance[2*205] + verlet_distance[2*205+1]*verlet_distance[2*205+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*205] 	-= floor(position[2*205]/L_x)*L_x;

		xi = position[2*206];
		yi = position[2*206+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*206]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*206+1] + weigh_brown_B * g2;

		position[2*206]		+= dx;
		position[2*206+1] 	+= dy;

		displacement[2*206]	+= dx;
		displacement[2*206+1] += dy;

		verlet_distance[2*206] 	+= (xi - position[2*206]);
		verlet_distance[2*206+1] 	+= (yi - position[2*206+1]);

		if ((temp = sqrt(verlet_distance[2*206]*verlet_distance[2*206] + verlet_distance[2*206+1]*verlet_distance[2*206+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*206] 	-= floor(position[2*206]/L_x)*L_x;

		xi = position[2*207];
		yi = position[2*207+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*207]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*207+1] + weigh_brown_B * g2;

		position[2*207]		+= dx;
		position[2*207+1] 	+= dy;

		displacement[2*207]	+= dx;
		displacement[2*207+1] += dy;

		verlet_distance[2*207] 	+= (xi - position[2*207]);
		verlet_distance[2*207+1] 	+= (yi - position[2*207+1]);

		if ((temp = sqrt(verlet_distance[2*207]*verlet_distance[2*207] + verlet_distance[2*207+1]*verlet_distance[2*207+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*207] 	-= floor(position[2*207]/L_x)*L_x;

		xi = position[2*208];
		yi = position[2*208+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*208]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*208+1] + weigh_brown_B * g2;

		position[2*208]		+= dx;
		position[2*208+1] 	+= dy;

		displacement[2*208]	+= dx;
		displacement[2*208+1] += dy;

		verlet_distance[2*208] 	+= (xi - position[2*208]);
		verlet_distance[2*208+1] 	+= (yi - position[2*208+1]);

		if ((temp = sqrt(verlet_distance[2*208]*verlet_distance[2*208] + verlet_distance[2*208+1]*verlet_distance[2*208+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*208] 	-= floor(position[2*208]/L_x)*L_x;

		xi = position[2*209];
		yi = position[2*209+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*209]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*209+1] + weigh_brown_B * g2;

		position[2*209]		+= dx;
		position[2*209+1] 	+= dy;

		displacement[2*209]	+= dx;
		displacement[2*209+1] += dy;

		verlet_distance[2*209] 	+= (xi - position[2*209]);
		verlet_distance[2*209+1] 	+= (yi - position[2*209+1]);

		if ((temp = sqrt(verlet_distance[2*209]*verlet_distance[2*209] + verlet_distance[2*209+1]*verlet_distance[2*209+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*209] 	-= floor(position[2*209]/L_x)*L_x;

		xi = position[2*210];
		yi = position[2*210+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*210]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*210+1] + weigh_brown_B * g2;

		position[2*210]		+= dx;
		position[2*210+1] 	+= dy;

		displacement[2*210]	+= dx;
		displacement[2*210+1] += dy;

		verlet_distance[2*210] 	+= (xi - position[2*210]);
		verlet_distance[2*210+1] 	+= (yi - position[2*210+1]);

		if ((temp = sqrt(verlet_distance[2*210]*verlet_distance[2*210] + verlet_distance[2*210+1]*verlet_distance[2*210+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*210] 	-= floor(position[2*210]/L_x)*L_x;

		xi = position[2*211];
		yi = position[2*211+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*211]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*211+1] + weigh_brown_B * g2;

		position[2*211]		+= dx;
		position[2*211+1] 	+= dy;

		displacement[2*211]	+= dx;
		displacement[2*211+1] += dy;

		verlet_distance[2*211] 	+= (xi - position[2*211]);
		verlet_distance[2*211+1] 	+= (yi - position[2*211+1]);

		if ((temp = sqrt(verlet_distance[2*211]*verlet_distance[2*211] + verlet_distance[2*211+1]*verlet_distance[2*211+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*211] 	-= floor(position[2*211]/L_x)*L_x;

		xi = position[2*212];
		yi = position[2*212+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*212]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*212+1] + weigh_brown_B * g2;

		position[2*212]		+= dx;
		position[2*212+1] 	+= dy;

		displacement[2*212]	+= dx;
		displacement[2*212+1] += dy;

		verlet_distance[2*212] 	+= (xi - position[2*212]);
		verlet_distance[2*212+1] 	+= (yi - position[2*212+1]);

		if ((temp = sqrt(verlet_distance[2*212]*verlet_distance[2*212] + verlet_distance[2*212+1]*verlet_distance[2*212+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*212] 	-= floor(position[2*212]/L_x)*L_x;

		xi = position[2*213];
		yi = position[2*213+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*213]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*213+1] + weigh_brown_B * g2;

		position[2*213]		+= dx;
		position[2*213+1] 	+= dy;

		displacement[2*213]	+= dx;
		displacement[2*213+1] += dy;

		verlet_distance[2*213] 	+= (xi - position[2*213]);
		verlet_distance[2*213+1] 	+= (yi - position[2*213+1]);

		if ((temp = sqrt(verlet_distance[2*213]*verlet_distance[2*213] + verlet_distance[2*213+1]*verlet_distance[2*213+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*213] 	-= floor(position[2*213]/L_x)*L_x;

		xi = position[2*214];
		yi = position[2*214+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*214]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*214+1] + weigh_brown_B * g2;

		position[2*214]		+= dx;
		position[2*214+1] 	+= dy;

		displacement[2*214]	+= dx;
		displacement[2*214+1] += dy;

		verlet_distance[2*214] 	+= (xi - position[2*214]);
		verlet_distance[2*214+1] 	+= (yi - position[2*214+1]);

		if ((temp = sqrt(verlet_distance[2*214]*verlet_distance[2*214] + verlet_distance[2*214+1]*verlet_distance[2*214+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*214] 	-= floor(position[2*214]/L_x)*L_x;

		xi = position[2*215];
		yi = position[2*215+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*215]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*215+1] + weigh_brown_B * g2;

		position[2*215]		+= dx;
		position[2*215+1] 	+= dy;

		displacement[2*215]	+= dx;
		displacement[2*215+1] += dy;

		verlet_distance[2*215] 	+= (xi - position[2*215]);
		verlet_distance[2*215+1] 	+= (yi - position[2*215+1]);

		if ((temp = sqrt(verlet_distance[2*215]*verlet_distance[2*215] + verlet_distance[2*215+1]*verlet_distance[2*215+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*215] 	-= floor(position[2*215]/L_x)*L_x;

		xi = position[2*216];
		yi = position[2*216+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*216]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*216+1] + weigh_brown_B * g2;

		position[2*216]		+= dx;
		position[2*216+1] 	+= dy;

		displacement[2*216]	+= dx;
		displacement[2*216+1] += dy;

		verlet_distance[2*216] 	+= (xi - position[2*216]);
		verlet_distance[2*216+1] 	+= (yi - position[2*216+1]);

		if ((temp = sqrt(verlet_distance[2*216]*verlet_distance[2*216] + verlet_distance[2*216+1]*verlet_distance[2*216+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*216] 	-= floor(position[2*216]/L_x)*L_x;

		xi = position[2*217];
		yi = position[2*217+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*217]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*217+1] + weigh_brown_B * g2;

		position[2*217]		+= dx;
		position[2*217+1] 	+= dy;

		displacement[2*217]	+= dx;
		displacement[2*217+1] += dy;

		verlet_distance[2*217] 	+= (xi - position[2*217]);
		verlet_distance[2*217+1] 	+= (yi - position[2*217+1]);

		if ((temp = sqrt(verlet_distance[2*217]*verlet_distance[2*217] + verlet_distance[2*217+1]*verlet_distance[2*217+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*217] 	-= floor(position[2*217]/L_x)*L_x;

		xi = position[2*218];
		yi = position[2*218+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*218]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*218+1] + weigh_brown_B * g2;

		position[2*218]		+= dx;
		position[2*218+1] 	+= dy;

		displacement[2*218]	+= dx;
		displacement[2*218+1] += dy;

		verlet_distance[2*218] 	+= (xi - position[2*218]);
		verlet_distance[2*218+1] 	+= (yi - position[2*218+1]);

		if ((temp = sqrt(verlet_distance[2*218]*verlet_distance[2*218] + verlet_distance[2*218+1]*verlet_distance[2*218+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*218] 	-= floor(position[2*218]/L_x)*L_x;

		xi = position[2*219];
		yi = position[2*219+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*219]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*219+1] + weigh_brown_B * g2;

		position[2*219]		+= dx;
		position[2*219+1] 	+= dy;

		displacement[2*219]	+= dx;
		displacement[2*219+1] += dy;

		verlet_distance[2*219] 	+= (xi - position[2*219]);
		verlet_distance[2*219+1] 	+= (yi - position[2*219+1]);

		if ((temp = sqrt(verlet_distance[2*219]*verlet_distance[2*219] + verlet_distance[2*219+1]*verlet_distance[2*219+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*219] 	-= floor(position[2*219]/L_x)*L_x;

		xi = position[2*220];
		yi = position[2*220+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*220]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*220+1] + weigh_brown_B * g2;

		position[2*220]		+= dx;
		position[2*220+1] 	+= dy;

		displacement[2*220]	+= dx;
		displacement[2*220+1] += dy;

		verlet_distance[2*220] 	+= (xi - position[2*220]);
		verlet_distance[2*220+1] 	+= (yi - position[2*220+1]);

		if ((temp = sqrt(verlet_distance[2*220]*verlet_distance[2*220] + verlet_distance[2*220+1]*verlet_distance[2*220+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*220] 	-= floor(position[2*220]/L_x)*L_x;

		xi = position[2*221];
		yi = position[2*221+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*221]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*221+1] + weigh_brown_B * g2;

		position[2*221]		+= dx;
		position[2*221+1] 	+= dy;

		displacement[2*221]	+= dx;
		displacement[2*221+1] += dy;

		verlet_distance[2*221] 	+= (xi - position[2*221]);
		verlet_distance[2*221+1] 	+= (yi - position[2*221+1]);

		if ((temp = sqrt(verlet_distance[2*221]*verlet_distance[2*221] + verlet_distance[2*221+1]*verlet_distance[2*221+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*221] 	-= floor(position[2*221]/L_x)*L_x;

		xi = position[2*222];
		yi = position[2*222+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*222]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*222+1] + weigh_brown_B * g2;

		position[2*222]		+= dx;
		position[2*222+1] 	+= dy;

		displacement[2*222]	+= dx;
		displacement[2*222+1] += dy;

		verlet_distance[2*222] 	+= (xi - position[2*222]);
		verlet_distance[2*222+1] 	+= (yi - position[2*222+1]);

		if ((temp = sqrt(verlet_distance[2*222]*verlet_distance[2*222] + verlet_distance[2*222+1]*verlet_distance[2*222+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*222] 	-= floor(position[2*222]/L_x)*L_x;

		xi = position[2*223];
		yi = position[2*223+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*223]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*223+1] + weigh_brown_B * g2;

		position[2*223]		+= dx;
		position[2*223+1] 	+= dy;

		displacement[2*223]	+= dx;
		displacement[2*223+1] += dy;

		verlet_distance[2*223] 	+= (xi - position[2*223]);
		verlet_distance[2*223+1] 	+= (yi - position[2*223+1]);

		if ((temp = sqrt(verlet_distance[2*223]*verlet_distance[2*223] + verlet_distance[2*223+1]*verlet_distance[2*223+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*223] 	-= floor(position[2*223]/L_x)*L_x;

		xi = position[2*224];
		yi = position[2*224+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*224]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*224+1] + weigh_brown_B * g2;

		position[2*224]		+= dx;
		position[2*224+1] 	+= dy;

		displacement[2*224]	+= dx;
		displacement[2*224+1] += dy;

		verlet_distance[2*224] 	+= (xi - position[2*224]);
		verlet_distance[2*224+1] 	+= (yi - position[2*224+1]);

		if ((temp = sqrt(verlet_distance[2*224]*verlet_distance[2*224] + verlet_distance[2*224+1]*verlet_distance[2*224+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*224] 	-= floor(position[2*224]/L_x)*L_x;

		// Signal to main thread, that all threads have finished their iteration
		pthread_barrier_wait(&barrier_main_one);

		// thread is done with one iteration, it will now wait for a signal from the main thread to continue
		pthread_barrier_wait(&barrier_main_two);
	}

	return NULL;
}



static void *iteration_6 (int *no) {
	// Define necessary variables
	double xi, xj, yi, yj;
	double dx, dy;
	double r_squared;

	double temp_force;
	double temp;

	double u1, u2, g1, g2;

	int iterate;
	int j;

	double m_i_A = 1.0;
	double m_i_B = m;

	double D_kT_A = D_Brown_A/kT;
	double D_kT_B = D_Brown_B/kT;

	double shear_A = v_A * delta_t;

	double m_j;

	double dist_bottom;
	double dist_top;

	// cutoff for the walls force, this will lead to a smoother transition
	double dist_cutoff 			= 1.0;
	double prefactor			= 3.0;
	double force_wall_cutoff 	= 1.0*prefactor *(kappa/dist_cutoff + 1.0/(dist_cutoff*dist_cutoff))*exp(-kappa*dist_cutoff);

	while(cont == 1) {
		xi = position[2*225];
		yi = position[2*225+1];

		force[2*225]	 = 0;
		force[2*225+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*225+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*225+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*225+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*225+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(225+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*225+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*225] 		+= temp_force*dx;
			force[2*225+1]	+= temp_force*dy;
		}

		xi = position[2*226];
		yi = position[2*226+1];

		force[2*226]	 = 0;
		force[2*226+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*226+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*226+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*226+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*226+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(226+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*226+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*226] 		+= temp_force*dx;
			force[2*226+1]	+= temp_force*dy;
		}

		xi = position[2*227];
		yi = position[2*227+1];

		force[2*227]	 = 0;
		force[2*227+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*227+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*227+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*227+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*227+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(227+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*227+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*227] 		+= temp_force*dx;
			force[2*227+1]	+= temp_force*dy;
		}

		xi = position[2*228];
		yi = position[2*228+1];

		force[2*228]	 = 0;
		force[2*228+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*228+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*228+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*228+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*228+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(228+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*228+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*228] 		+= temp_force*dx;
			force[2*228+1]	+= temp_force*dy;
		}

		xi = position[2*229];
		yi = position[2*229+1];

		force[2*229]	 = 0;
		force[2*229+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*229+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*229+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*229+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*229+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(229+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*229+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*229] 		+= temp_force*dx;
			force[2*229+1]	+= temp_force*dy;
		}

		xi = position[2*230];
		yi = position[2*230+1];

		force[2*230]	 = 0;
		force[2*230+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*230+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*230+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*230+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*230+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(230+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*230+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*230] 		+= temp_force*dx;
			force[2*230+1]	+= temp_force*dy;
		}

		xi = position[2*231];
		yi = position[2*231+1];

		force[2*231]	 = 0;
		force[2*231+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*231+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*231+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*231+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*231+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(231+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*231+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*231] 		+= temp_force*dx;
			force[2*231+1]	+= temp_force*dy;
		}

		xi = position[2*232];
		yi = position[2*232+1];

		force[2*232]	 = 0;
		force[2*232+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*232+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*232+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*232+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*232+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(232+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*232+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*232] 		+= temp_force*dx;
			force[2*232+1]	+= temp_force*dy;
		}

		xi = position[2*233];
		yi = position[2*233+1];

		force[2*233]	 = 0;
		force[2*233+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*233+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*233+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*233+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*233+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(233+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*233+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*233] 		+= temp_force*dx;
			force[2*233+1]	+= temp_force*dy;
		}

		xi = position[2*234];
		yi = position[2*234+1];

		force[2*234]	 = 0;
		force[2*234+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*234+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*234+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*234+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*234+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(234+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*234+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*234] 		+= temp_force*dx;
			force[2*234+1]	+= temp_force*dy;
		}

		xi = position[2*235];
		yi = position[2*235+1];

		force[2*235]	 = 0;
		force[2*235+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*235+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*235+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*235+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*235+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(235+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*235+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*235] 		+= temp_force*dx;
			force[2*235+1]	+= temp_force*dy;
		}

		xi = position[2*236];
		yi = position[2*236+1];

		force[2*236]	 = 0;
		force[2*236+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*236+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*236+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*236+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*236+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(236+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*236+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*236] 		+= temp_force*dx;
			force[2*236+1]	+= temp_force*dy;
		}

		xi = position[2*237];
		yi = position[2*237+1];

		force[2*237]	 = 0;
		force[2*237+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*237+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*237+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*237+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*237+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(237+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*237+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*237] 		+= temp_force*dx;
			force[2*237+1]	+= temp_force*dy;
		}

		xi = position[2*238];
		yi = position[2*238+1];

		force[2*238]	 = 0;
		force[2*238+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*238+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*238+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*238+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*238+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(238+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*238+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*238] 		+= temp_force*dx;
			force[2*238+1]	+= temp_force*dy;
		}

		xi = position[2*239];
		yi = position[2*239+1];

		force[2*239]	 = 0;
		force[2*239+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*239+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*239+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*239+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*239+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(239+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*239+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*239] 		+= temp_force*dx;
			force[2*239+1]	+= temp_force*dy;
		}

		xi = position[2*240];
		yi = position[2*240+1];

		force[2*240]	 = 0;
		force[2*240+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*240+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*240+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*240+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*240+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(240+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*240+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*240] 		+= temp_force*dx;
			force[2*240+1]	+= temp_force*dy;
		}

		xi = position[2*241];
		yi = position[2*241+1];

		force[2*241]	 = 0;
		force[2*241+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*241+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*241+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*241+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*241+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(241+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*241+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*241] 		+= temp_force*dx;
			force[2*241+1]	+= temp_force*dy;
		}

		xi = position[2*242];
		yi = position[2*242+1];

		force[2*242]	 = 0;
		force[2*242+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*242+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*242+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*242+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*242+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(242+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*242+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*242] 		+= temp_force*dx;
			force[2*242+1]	+= temp_force*dy;
		}

		xi = position[2*243];
		yi = position[2*243+1];

		force[2*243]	 = 0;
		force[2*243+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*243+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*243+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*243+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*243+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(243+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*243+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*243] 		+= temp_force*dx;
			force[2*243+1]	+= temp_force*dy;
		}

		xi = position[2*244];
		yi = position[2*244+1];

		force[2*244]	 = 0;
		force[2*244+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*244+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*244+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*244+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*244+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(244+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*244+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*244] 		+= temp_force*dx;
			force[2*244+1]	+= temp_force*dy;
		}

		xi = position[2*245];
		yi = position[2*245+1];

		force[2*245]	 = 0;
		force[2*245+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*245+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*245+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*245+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*245+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(245+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*245+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*245] 		+= temp_force*dx;
			force[2*245+1]	+= temp_force*dy;
		}

		xi = position[2*246];
		yi = position[2*246+1];

		force[2*246]	 = 0;
		force[2*246+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*246+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*246+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*246+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*246+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(246+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*246+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*246] 		+= temp_force*dx;
			force[2*246+1]	+= temp_force*dy;
		}

		xi = position[2*247];
		yi = position[2*247+1];

		force[2*247]	 = 0;
		force[2*247+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*247+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*247+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*247+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*247+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(247+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*247+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*247] 		+= temp_force*dx;
			force[2*247+1]	+= temp_force*dy;
		}

		xi = position[2*248];
		yi = position[2*248+1];

		force[2*248]	 = 0;
		force[2*248+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*248+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*248+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*248+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*248+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(248+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*248+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*248] 		+= temp_force*dx;
			force[2*248+1]	+= temp_force*dy;
		}

		xi = position[2*249];
		yi = position[2*249+1];

		force[2*249]	 = 0;
		force[2*249+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*249+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*249+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*249+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*249+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(249+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*249+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*249] 		+= temp_force*dx;
			force[2*249+1]	+= temp_force*dy;
		}

		xi = position[2*250];
		yi = position[2*250+1];

		force[2*250]	 = 0;
		force[2*250+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*250+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*250+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*250+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*250+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(250+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*250+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*250] 		+= temp_force*dx;
			force[2*250+1]	+= temp_force*dy;
		}

		xi = position[2*251];
		yi = position[2*251+1];

		force[2*251]	 = 0;
		force[2*251+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*251+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*251+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*251+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*251+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(251+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*251+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*251] 		+= temp_force*dx;
			force[2*251+1]	+= temp_force*dy;
		}

		xi = position[2*252];
		yi = position[2*252+1];

		force[2*252]	 = 0;
		force[2*252+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*252+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*252+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*252+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*252+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(252+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*252+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*252] 		+= temp_force*dx;
			force[2*252+1]	+= temp_force*dy;
		}

		xi = position[2*253];
		yi = position[2*253+1];

		force[2*253]	 = 0;
		force[2*253+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*253+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*253+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*253+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*253+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(253+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*253+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*253] 		+= temp_force*dx;
			force[2*253+1]	+= temp_force*dy;
		}

		xi = position[2*254];
		yi = position[2*254+1];

		force[2*254]	 = 0;
		force[2*254+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*254+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*254+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*254+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*254+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(254+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*254+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*254] 		+= temp_force*dx;
			force[2*254+1]	+= temp_force*dy;
		}

		xi = position[2*255];
		yi = position[2*255+1];

		force[2*255]	 = 0;
		force[2*255+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*255+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*255+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*255+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*255+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(255+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*255+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*255] 		+= temp_force*dx;
			force[2*255+1]	+= temp_force*dy;
		}

		xi = position[2*256];
		yi = position[2*256+1];

		force[2*256]	 = 0;
		force[2*256+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*256+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*256+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*256+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*256+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(256+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*256+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*256] 		+= temp_force*dx;
			force[2*256+1]	+= temp_force*dy;
		}

		xi = position[2*257];
		yi = position[2*257+1];

		force[2*257]	 = 0;
		force[2*257+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*257+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*257+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*257+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*257+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(257+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*257+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*257] 		+= temp_force*dx;
			force[2*257+1]	+= temp_force*dy;
		}

		xi = position[2*258];
		yi = position[2*258+1];

		force[2*258]	 = 0;
		force[2*258+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*258+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*258+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*258+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*258+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(258+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*258+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*258] 		+= temp_force*dx;
			force[2*258+1]	+= temp_force*dy;
		}

		xi = position[2*259];
		yi = position[2*259+1];

		force[2*259]	 = 0;
		force[2*259+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*259+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*259+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*259+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*259+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(259+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*259+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*259] 		+= temp_force*dx;
			force[2*259+1]	+= temp_force*dy;
		}

		xi = position[2*260];
		yi = position[2*260+1];

		force[2*260]	 = 0;
		force[2*260+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*260+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*260+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*260+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*260+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(260+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*260+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*260] 		+= temp_force*dx;
			force[2*260+1]	+= temp_force*dy;
		}

		xi = position[2*261];
		yi = position[2*261+1];

		force[2*261]	 = 0;
		force[2*261+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*261+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*261+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*261+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*261+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(261+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*261+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*261] 		+= temp_force*dx;
			force[2*261+1]	+= temp_force*dy;
		}

		xi = position[2*262];
		yi = position[2*262+1];

		force[2*262]	 = 0;
		force[2*262+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*262+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*262+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*262+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*262+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(262+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*262+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*262] 		+= temp_force*dx;
			force[2*262+1]	+= temp_force*dy;
		}

		pthread_barrier_wait(&barrier_internal);

		xi = position[2*225];
		yi = position[2*225+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*225]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*225+1] + weigh_brown_B * g2;

		position[2*225]		+= dx;
		position[2*225+1] 	+= dy;

		displacement[2*225]	+= dx;
		displacement[2*225+1] += dy;

		verlet_distance[2*225] 	+= (xi - position[2*225]);
		verlet_distance[2*225+1] 	+= (yi - position[2*225+1]);

		if ((temp = sqrt(verlet_distance[2*225]*verlet_distance[2*225] + verlet_distance[2*225+1]*verlet_distance[2*225+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*225] 	-= floor(position[2*225]/L_x)*L_x;

		xi = position[2*226];
		yi = position[2*226+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*226]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*226+1] + weigh_brown_B * g2;

		position[2*226]		+= dx;
		position[2*226+1] 	+= dy;

		displacement[2*226]	+= dx;
		displacement[2*226+1] += dy;

		verlet_distance[2*226] 	+= (xi - position[2*226]);
		verlet_distance[2*226+1] 	+= (yi - position[2*226+1]);

		if ((temp = sqrt(verlet_distance[2*226]*verlet_distance[2*226] + verlet_distance[2*226+1]*verlet_distance[2*226+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*226] 	-= floor(position[2*226]/L_x)*L_x;

		xi = position[2*227];
		yi = position[2*227+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*227]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*227+1] + weigh_brown_B * g2;

		position[2*227]		+= dx;
		position[2*227+1] 	+= dy;

		displacement[2*227]	+= dx;
		displacement[2*227+1] += dy;

		verlet_distance[2*227] 	+= (xi - position[2*227]);
		verlet_distance[2*227+1] 	+= (yi - position[2*227+1]);

		if ((temp = sqrt(verlet_distance[2*227]*verlet_distance[2*227] + verlet_distance[2*227+1]*verlet_distance[2*227+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*227] 	-= floor(position[2*227]/L_x)*L_x;

		xi = position[2*228];
		yi = position[2*228+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*228]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*228+1] + weigh_brown_B * g2;

		position[2*228]		+= dx;
		position[2*228+1] 	+= dy;

		displacement[2*228]	+= dx;
		displacement[2*228+1] += dy;

		verlet_distance[2*228] 	+= (xi - position[2*228]);
		verlet_distance[2*228+1] 	+= (yi - position[2*228+1]);

		if ((temp = sqrt(verlet_distance[2*228]*verlet_distance[2*228] + verlet_distance[2*228+1]*verlet_distance[2*228+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*228] 	-= floor(position[2*228]/L_x)*L_x;

		xi = position[2*229];
		yi = position[2*229+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*229]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*229+1] + weigh_brown_B * g2;

		position[2*229]		+= dx;
		position[2*229+1] 	+= dy;

		displacement[2*229]	+= dx;
		displacement[2*229+1] += dy;

		verlet_distance[2*229] 	+= (xi - position[2*229]);
		verlet_distance[2*229+1] 	+= (yi - position[2*229+1]);

		if ((temp = sqrt(verlet_distance[2*229]*verlet_distance[2*229] + verlet_distance[2*229+1]*verlet_distance[2*229+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*229] 	-= floor(position[2*229]/L_x)*L_x;

		xi = position[2*230];
		yi = position[2*230+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*230]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*230+1] + weigh_brown_B * g2;

		position[2*230]		+= dx;
		position[2*230+1] 	+= dy;

		displacement[2*230]	+= dx;
		displacement[2*230+1] += dy;

		verlet_distance[2*230] 	+= (xi - position[2*230]);
		verlet_distance[2*230+1] 	+= (yi - position[2*230+1]);

		if ((temp = sqrt(verlet_distance[2*230]*verlet_distance[2*230] + verlet_distance[2*230+1]*verlet_distance[2*230+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*230] 	-= floor(position[2*230]/L_x)*L_x;

		xi = position[2*231];
		yi = position[2*231+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*231]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*231+1] + weigh_brown_B * g2;

		position[2*231]		+= dx;
		position[2*231+1] 	+= dy;

		displacement[2*231]	+= dx;
		displacement[2*231+1] += dy;

		verlet_distance[2*231] 	+= (xi - position[2*231]);
		verlet_distance[2*231+1] 	+= (yi - position[2*231+1]);

		if ((temp = sqrt(verlet_distance[2*231]*verlet_distance[2*231] + verlet_distance[2*231+1]*verlet_distance[2*231+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*231] 	-= floor(position[2*231]/L_x)*L_x;

		xi = position[2*232];
		yi = position[2*232+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*232]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*232+1] + weigh_brown_B * g2;

		position[2*232]		+= dx;
		position[2*232+1] 	+= dy;

		displacement[2*232]	+= dx;
		displacement[2*232+1] += dy;

		verlet_distance[2*232] 	+= (xi - position[2*232]);
		verlet_distance[2*232+1] 	+= (yi - position[2*232+1]);

		if ((temp = sqrt(verlet_distance[2*232]*verlet_distance[2*232] + verlet_distance[2*232+1]*verlet_distance[2*232+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*232] 	-= floor(position[2*232]/L_x)*L_x;

		xi = position[2*233];
		yi = position[2*233+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*233]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*233+1] + weigh_brown_B * g2;

		position[2*233]		+= dx;
		position[2*233+1] 	+= dy;

		displacement[2*233]	+= dx;
		displacement[2*233+1] += dy;

		verlet_distance[2*233] 	+= (xi - position[2*233]);
		verlet_distance[2*233+1] 	+= (yi - position[2*233+1]);

		if ((temp = sqrt(verlet_distance[2*233]*verlet_distance[2*233] + verlet_distance[2*233+1]*verlet_distance[2*233+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*233] 	-= floor(position[2*233]/L_x)*L_x;

		xi = position[2*234];
		yi = position[2*234+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*234]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*234+1] + weigh_brown_B * g2;

		position[2*234]		+= dx;
		position[2*234+1] 	+= dy;

		displacement[2*234]	+= dx;
		displacement[2*234+1] += dy;

		verlet_distance[2*234] 	+= (xi - position[2*234]);
		verlet_distance[2*234+1] 	+= (yi - position[2*234+1]);

		if ((temp = sqrt(verlet_distance[2*234]*verlet_distance[2*234] + verlet_distance[2*234+1]*verlet_distance[2*234+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*234] 	-= floor(position[2*234]/L_x)*L_x;

		xi = position[2*235];
		yi = position[2*235+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*235]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*235+1] + weigh_brown_B * g2;

		position[2*235]		+= dx;
		position[2*235+1] 	+= dy;

		displacement[2*235]	+= dx;
		displacement[2*235+1] += dy;

		verlet_distance[2*235] 	+= (xi - position[2*235]);
		verlet_distance[2*235+1] 	+= (yi - position[2*235+1]);

		if ((temp = sqrt(verlet_distance[2*235]*verlet_distance[2*235] + verlet_distance[2*235+1]*verlet_distance[2*235+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*235] 	-= floor(position[2*235]/L_x)*L_x;

		xi = position[2*236];
		yi = position[2*236+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*236]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*236+1] + weigh_brown_B * g2;

		position[2*236]		+= dx;
		position[2*236+1] 	+= dy;

		displacement[2*236]	+= dx;
		displacement[2*236+1] += dy;

		verlet_distance[2*236] 	+= (xi - position[2*236]);
		verlet_distance[2*236+1] 	+= (yi - position[2*236+1]);

		if ((temp = sqrt(verlet_distance[2*236]*verlet_distance[2*236] + verlet_distance[2*236+1]*verlet_distance[2*236+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*236] 	-= floor(position[2*236]/L_x)*L_x;

		xi = position[2*237];
		yi = position[2*237+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*237]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*237+1] + weigh_brown_B * g2;

		position[2*237]		+= dx;
		position[2*237+1] 	+= dy;

		displacement[2*237]	+= dx;
		displacement[2*237+1] += dy;

		verlet_distance[2*237] 	+= (xi - position[2*237]);
		verlet_distance[2*237+1] 	+= (yi - position[2*237+1]);

		if ((temp = sqrt(verlet_distance[2*237]*verlet_distance[2*237] + verlet_distance[2*237+1]*verlet_distance[2*237+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*237] 	-= floor(position[2*237]/L_x)*L_x;

		xi = position[2*238];
		yi = position[2*238+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*238]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*238+1] + weigh_brown_B * g2;

		position[2*238]		+= dx;
		position[2*238+1] 	+= dy;

		displacement[2*238]	+= dx;
		displacement[2*238+1] += dy;

		verlet_distance[2*238] 	+= (xi - position[2*238]);
		verlet_distance[2*238+1] 	+= (yi - position[2*238+1]);

		if ((temp = sqrt(verlet_distance[2*238]*verlet_distance[2*238] + verlet_distance[2*238+1]*verlet_distance[2*238+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*238] 	-= floor(position[2*238]/L_x)*L_x;

		xi = position[2*239];
		yi = position[2*239+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*239]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*239+1] + weigh_brown_B * g2;

		position[2*239]		+= dx;
		position[2*239+1] 	+= dy;

		displacement[2*239]	+= dx;
		displacement[2*239+1] += dy;

		verlet_distance[2*239] 	+= (xi - position[2*239]);
		verlet_distance[2*239+1] 	+= (yi - position[2*239+1]);

		if ((temp = sqrt(verlet_distance[2*239]*verlet_distance[2*239] + verlet_distance[2*239+1]*verlet_distance[2*239+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*239] 	-= floor(position[2*239]/L_x)*L_x;

		xi = position[2*240];
		yi = position[2*240+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*240]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*240+1] + weigh_brown_B * g2;

		position[2*240]		+= dx;
		position[2*240+1] 	+= dy;

		displacement[2*240]	+= dx;
		displacement[2*240+1] += dy;

		verlet_distance[2*240] 	+= (xi - position[2*240]);
		verlet_distance[2*240+1] 	+= (yi - position[2*240+1]);

		if ((temp = sqrt(verlet_distance[2*240]*verlet_distance[2*240] + verlet_distance[2*240+1]*verlet_distance[2*240+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*240] 	-= floor(position[2*240]/L_x)*L_x;

		xi = position[2*241];
		yi = position[2*241+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*241]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*241+1] + weigh_brown_B * g2;

		position[2*241]		+= dx;
		position[2*241+1] 	+= dy;

		displacement[2*241]	+= dx;
		displacement[2*241+1] += dy;

		verlet_distance[2*241] 	+= (xi - position[2*241]);
		verlet_distance[2*241+1] 	+= (yi - position[2*241+1]);

		if ((temp = sqrt(verlet_distance[2*241]*verlet_distance[2*241] + verlet_distance[2*241+1]*verlet_distance[2*241+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*241] 	-= floor(position[2*241]/L_x)*L_x;

		xi = position[2*242];
		yi = position[2*242+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*242]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*242+1] + weigh_brown_B * g2;

		position[2*242]		+= dx;
		position[2*242+1] 	+= dy;

		displacement[2*242]	+= dx;
		displacement[2*242+1] += dy;

		verlet_distance[2*242] 	+= (xi - position[2*242]);
		verlet_distance[2*242+1] 	+= (yi - position[2*242+1]);

		if ((temp = sqrt(verlet_distance[2*242]*verlet_distance[2*242] + verlet_distance[2*242+1]*verlet_distance[2*242+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*242] 	-= floor(position[2*242]/L_x)*L_x;

		xi = position[2*243];
		yi = position[2*243+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*243]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*243+1] + weigh_brown_B * g2;

		position[2*243]		+= dx;
		position[2*243+1] 	+= dy;

		displacement[2*243]	+= dx;
		displacement[2*243+1] += dy;

		verlet_distance[2*243] 	+= (xi - position[2*243]);
		verlet_distance[2*243+1] 	+= (yi - position[2*243+1]);

		if ((temp = sqrt(verlet_distance[2*243]*verlet_distance[2*243] + verlet_distance[2*243+1]*verlet_distance[2*243+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*243] 	-= floor(position[2*243]/L_x)*L_x;

		xi = position[2*244];
		yi = position[2*244+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*244]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*244+1] + weigh_brown_B * g2;

		position[2*244]		+= dx;
		position[2*244+1] 	+= dy;

		displacement[2*244]	+= dx;
		displacement[2*244+1] += dy;

		verlet_distance[2*244] 	+= (xi - position[2*244]);
		verlet_distance[2*244+1] 	+= (yi - position[2*244+1]);

		if ((temp = sqrt(verlet_distance[2*244]*verlet_distance[2*244] + verlet_distance[2*244+1]*verlet_distance[2*244+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*244] 	-= floor(position[2*244]/L_x)*L_x;

		xi = position[2*245];
		yi = position[2*245+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*245]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*245+1] + weigh_brown_B * g2;

		position[2*245]		+= dx;
		position[2*245+1] 	+= dy;

		displacement[2*245]	+= dx;
		displacement[2*245+1] += dy;

		verlet_distance[2*245] 	+= (xi - position[2*245]);
		verlet_distance[2*245+1] 	+= (yi - position[2*245+1]);

		if ((temp = sqrt(verlet_distance[2*245]*verlet_distance[2*245] + verlet_distance[2*245+1]*verlet_distance[2*245+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*245] 	-= floor(position[2*245]/L_x)*L_x;

		xi = position[2*246];
		yi = position[2*246+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*246]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*246+1] + weigh_brown_B * g2;

		position[2*246]		+= dx;
		position[2*246+1] 	+= dy;

		displacement[2*246]	+= dx;
		displacement[2*246+1] += dy;

		verlet_distance[2*246] 	+= (xi - position[2*246]);
		verlet_distance[2*246+1] 	+= (yi - position[2*246+1]);

		if ((temp = sqrt(verlet_distance[2*246]*verlet_distance[2*246] + verlet_distance[2*246+1]*verlet_distance[2*246+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*246] 	-= floor(position[2*246]/L_x)*L_x;

		xi = position[2*247];
		yi = position[2*247+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*247]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*247+1] + weigh_brown_B * g2;

		position[2*247]		+= dx;
		position[2*247+1] 	+= dy;

		displacement[2*247]	+= dx;
		displacement[2*247+1] += dy;

		verlet_distance[2*247] 	+= (xi - position[2*247]);
		verlet_distance[2*247+1] 	+= (yi - position[2*247+1]);

		if ((temp = sqrt(verlet_distance[2*247]*verlet_distance[2*247] + verlet_distance[2*247+1]*verlet_distance[2*247+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*247] 	-= floor(position[2*247]/L_x)*L_x;

		xi = position[2*248];
		yi = position[2*248+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*248]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*248+1] + weigh_brown_B * g2;

		position[2*248]		+= dx;
		position[2*248+1] 	+= dy;

		displacement[2*248]	+= dx;
		displacement[2*248+1] += dy;

		verlet_distance[2*248] 	+= (xi - position[2*248]);
		verlet_distance[2*248+1] 	+= (yi - position[2*248+1]);

		if ((temp = sqrt(verlet_distance[2*248]*verlet_distance[2*248] + verlet_distance[2*248+1]*verlet_distance[2*248+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*248] 	-= floor(position[2*248]/L_x)*L_x;

		xi = position[2*249];
		yi = position[2*249+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*249]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*249+1] + weigh_brown_B * g2;

		position[2*249]		+= dx;
		position[2*249+1] 	+= dy;

		displacement[2*249]	+= dx;
		displacement[2*249+1] += dy;

		verlet_distance[2*249] 	+= (xi - position[2*249]);
		verlet_distance[2*249+1] 	+= (yi - position[2*249+1]);

		if ((temp = sqrt(verlet_distance[2*249]*verlet_distance[2*249] + verlet_distance[2*249+1]*verlet_distance[2*249+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*249] 	-= floor(position[2*249]/L_x)*L_x;

		xi = position[2*250];
		yi = position[2*250+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*250]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*250+1] + weigh_brown_B * g2;

		position[2*250]		+= dx;
		position[2*250+1] 	+= dy;

		displacement[2*250]	+= dx;
		displacement[2*250+1] += dy;

		verlet_distance[2*250] 	+= (xi - position[2*250]);
		verlet_distance[2*250+1] 	+= (yi - position[2*250+1]);

		if ((temp = sqrt(verlet_distance[2*250]*verlet_distance[2*250] + verlet_distance[2*250+1]*verlet_distance[2*250+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*250] 	-= floor(position[2*250]/L_x)*L_x;

		xi = position[2*251];
		yi = position[2*251+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*251]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*251+1] + weigh_brown_B * g2;

		position[2*251]		+= dx;
		position[2*251+1] 	+= dy;

		displacement[2*251]	+= dx;
		displacement[2*251+1] += dy;

		verlet_distance[2*251] 	+= (xi - position[2*251]);
		verlet_distance[2*251+1] 	+= (yi - position[2*251+1]);

		if ((temp = sqrt(verlet_distance[2*251]*verlet_distance[2*251] + verlet_distance[2*251+1]*verlet_distance[2*251+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*251] 	-= floor(position[2*251]/L_x)*L_x;

		xi = position[2*252];
		yi = position[2*252+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*252]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*252+1] + weigh_brown_B * g2;

		position[2*252]		+= dx;
		position[2*252+1] 	+= dy;

		displacement[2*252]	+= dx;
		displacement[2*252+1] += dy;

		verlet_distance[2*252] 	+= (xi - position[2*252]);
		verlet_distance[2*252+1] 	+= (yi - position[2*252+1]);

		if ((temp = sqrt(verlet_distance[2*252]*verlet_distance[2*252] + verlet_distance[2*252+1]*verlet_distance[2*252+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*252] 	-= floor(position[2*252]/L_x)*L_x;

		xi = position[2*253];
		yi = position[2*253+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*253]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*253+1] + weigh_brown_B * g2;

		position[2*253]		+= dx;
		position[2*253+1] 	+= dy;

		displacement[2*253]	+= dx;
		displacement[2*253+1] += dy;

		verlet_distance[2*253] 	+= (xi - position[2*253]);
		verlet_distance[2*253+1] 	+= (yi - position[2*253+1]);

		if ((temp = sqrt(verlet_distance[2*253]*verlet_distance[2*253] + verlet_distance[2*253+1]*verlet_distance[2*253+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*253] 	-= floor(position[2*253]/L_x)*L_x;

		xi = position[2*254];
		yi = position[2*254+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*254]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*254+1] + weigh_brown_B * g2;

		position[2*254]		+= dx;
		position[2*254+1] 	+= dy;

		displacement[2*254]	+= dx;
		displacement[2*254+1] += dy;

		verlet_distance[2*254] 	+= (xi - position[2*254]);
		verlet_distance[2*254+1] 	+= (yi - position[2*254+1]);

		if ((temp = sqrt(verlet_distance[2*254]*verlet_distance[2*254] + verlet_distance[2*254+1]*verlet_distance[2*254+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*254] 	-= floor(position[2*254]/L_x)*L_x;

		xi = position[2*255];
		yi = position[2*255+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*255]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*255+1] + weigh_brown_B * g2;

		position[2*255]		+= dx;
		position[2*255+1] 	+= dy;

		displacement[2*255]	+= dx;
		displacement[2*255+1] += dy;

		verlet_distance[2*255] 	+= (xi - position[2*255]);
		verlet_distance[2*255+1] 	+= (yi - position[2*255+1]);

		if ((temp = sqrt(verlet_distance[2*255]*verlet_distance[2*255] + verlet_distance[2*255+1]*verlet_distance[2*255+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*255] 	-= floor(position[2*255]/L_x)*L_x;

		xi = position[2*256];
		yi = position[2*256+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*256]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*256+1] + weigh_brown_B * g2;

		position[2*256]		+= dx;
		position[2*256+1] 	+= dy;

		displacement[2*256]	+= dx;
		displacement[2*256+1] += dy;

		verlet_distance[2*256] 	+= (xi - position[2*256]);
		verlet_distance[2*256+1] 	+= (yi - position[2*256+1]);

		if ((temp = sqrt(verlet_distance[2*256]*verlet_distance[2*256] + verlet_distance[2*256+1]*verlet_distance[2*256+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*256] 	-= floor(position[2*256]/L_x)*L_x;

		xi = position[2*257];
		yi = position[2*257+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*257]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*257+1] + weigh_brown_B * g2;

		position[2*257]		+= dx;
		position[2*257+1] 	+= dy;

		displacement[2*257]	+= dx;
		displacement[2*257+1] += dy;

		verlet_distance[2*257] 	+= (xi - position[2*257]);
		verlet_distance[2*257+1] 	+= (yi - position[2*257+1]);

		if ((temp = sqrt(verlet_distance[2*257]*verlet_distance[2*257] + verlet_distance[2*257+1]*verlet_distance[2*257+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*257] 	-= floor(position[2*257]/L_x)*L_x;

		xi = position[2*258];
		yi = position[2*258+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*258]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*258+1] + weigh_brown_B * g2;

		position[2*258]		+= dx;
		position[2*258+1] 	+= dy;

		displacement[2*258]	+= dx;
		displacement[2*258+1] += dy;

		verlet_distance[2*258] 	+= (xi - position[2*258]);
		verlet_distance[2*258+1] 	+= (yi - position[2*258+1]);

		if ((temp = sqrt(verlet_distance[2*258]*verlet_distance[2*258] + verlet_distance[2*258+1]*verlet_distance[2*258+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*258] 	-= floor(position[2*258]/L_x)*L_x;

		xi = position[2*259];
		yi = position[2*259+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*259]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*259+1] + weigh_brown_B * g2;

		position[2*259]		+= dx;
		position[2*259+1] 	+= dy;

		displacement[2*259]	+= dx;
		displacement[2*259+1] += dy;

		verlet_distance[2*259] 	+= (xi - position[2*259]);
		verlet_distance[2*259+1] 	+= (yi - position[2*259+1]);

		if ((temp = sqrt(verlet_distance[2*259]*verlet_distance[2*259] + verlet_distance[2*259+1]*verlet_distance[2*259+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*259] 	-= floor(position[2*259]/L_x)*L_x;

		xi = position[2*260];
		yi = position[2*260+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*260]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*260+1] + weigh_brown_B * g2;

		position[2*260]		+= dx;
		position[2*260+1] 	+= dy;

		displacement[2*260]	+= dx;
		displacement[2*260+1] += dy;

		verlet_distance[2*260] 	+= (xi - position[2*260]);
		verlet_distance[2*260+1] 	+= (yi - position[2*260+1]);

		if ((temp = sqrt(verlet_distance[2*260]*verlet_distance[2*260] + verlet_distance[2*260+1]*verlet_distance[2*260+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*260] 	-= floor(position[2*260]/L_x)*L_x;

		xi = position[2*261];
		yi = position[2*261+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*261]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*261+1] + weigh_brown_B * g2;

		position[2*261]		+= dx;
		position[2*261+1] 	+= dy;

		displacement[2*261]	+= dx;
		displacement[2*261+1] += dy;

		verlet_distance[2*261] 	+= (xi - position[2*261]);
		verlet_distance[2*261+1] 	+= (yi - position[2*261+1]);

		if ((temp = sqrt(verlet_distance[2*261]*verlet_distance[2*261] + verlet_distance[2*261+1]*verlet_distance[2*261+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*261] 	-= floor(position[2*261]/L_x)*L_x;

		xi = position[2*262];
		yi = position[2*262+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*262]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*262+1] + weigh_brown_B * g2;

		position[2*262]		+= dx;
		position[2*262+1] 	+= dy;

		displacement[2*262]	+= dx;
		displacement[2*262+1] += dy;

		verlet_distance[2*262] 	+= (xi - position[2*262]);
		verlet_distance[2*262+1] 	+= (yi - position[2*262+1]);

		if ((temp = sqrt(verlet_distance[2*262]*verlet_distance[2*262] + verlet_distance[2*262+1]*verlet_distance[2*262+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*262] 	-= floor(position[2*262]/L_x)*L_x;

		// Signal to main thread, that all threads have finished their iteration
		pthread_barrier_wait(&barrier_main_one);

		// thread is done with one iteration, it will now wait for a signal from the main thread to continue
		pthread_barrier_wait(&barrier_main_two);
	}

	return NULL;
}



static void *iteration_7 (int *no) {
	// Define necessary variables
	double xi, xj, yi, yj;
	double dx, dy;
	double r_squared;

	double temp_force;
	double temp;

	double u1, u2, g1, g2;

	int iterate;
	int j;

	double m_i_A = 1.0;
	double m_i_B = m;

	double D_kT_A = D_Brown_A/kT;
	double D_kT_B = D_Brown_B/kT;

	double shear_A = v_A * delta_t;

	double m_j;

	double dist_bottom;
	double dist_top;

	// cutoff for the walls force, this will lead to a smoother transition
	double dist_cutoff 			= 1.0;
	double prefactor			= 3.0;
	double force_wall_cutoff 	= 1.0*prefactor *(kappa/dist_cutoff + 1.0/(dist_cutoff*dist_cutoff))*exp(-kappa*dist_cutoff);

	while(cont == 1) {
		xi = position[2*263];
		yi = position[2*263+1];

		force[2*263]	 = 0;
		force[2*263+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*263+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*263+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*263+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*263+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(263+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*263+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*263] 		+= temp_force*dx;
			force[2*263+1]	+= temp_force*dy;
		}

		xi = position[2*264];
		yi = position[2*264+1];

		force[2*264]	 = 0;
		force[2*264+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*264+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*264+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*264+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*264+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(264+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*264+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*264] 		+= temp_force*dx;
			force[2*264+1]	+= temp_force*dy;
		}

		xi = position[2*265];
		yi = position[2*265+1];

		force[2*265]	 = 0;
		force[2*265+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*265+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*265+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*265+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*265+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(265+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*265+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*265] 		+= temp_force*dx;
			force[2*265+1]	+= temp_force*dy;
		}

		xi = position[2*266];
		yi = position[2*266+1];

		force[2*266]	 = 0;
		force[2*266+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*266+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*266+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*266+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*266+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(266+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*266+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*266] 		+= temp_force*dx;
			force[2*266+1]	+= temp_force*dy;
		}

		xi = position[2*267];
		yi = position[2*267+1];

		force[2*267]	 = 0;
		force[2*267+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*267+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*267+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*267+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*267+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(267+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*267+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*267] 		+= temp_force*dx;
			force[2*267+1]	+= temp_force*dy;
		}

		xi = position[2*268];
		yi = position[2*268+1];

		force[2*268]	 = 0;
		force[2*268+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*268+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*268+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*268+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*268+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(268+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*268+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*268] 		+= temp_force*dx;
			force[2*268+1]	+= temp_force*dy;
		}

		xi = position[2*269];
		yi = position[2*269+1];

		force[2*269]	 = 0;
		force[2*269+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*269+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*269+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*269+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*269+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(269+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*269+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*269] 		+= temp_force*dx;
			force[2*269+1]	+= temp_force*dy;
		}

		xi = position[2*270];
		yi = position[2*270+1];

		force[2*270]	 = 0;
		force[2*270+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*270+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*270+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*270+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*270+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(270+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*270+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*270] 		+= temp_force*dx;
			force[2*270+1]	+= temp_force*dy;
		}

		xi = position[2*271];
		yi = position[2*271+1];

		force[2*271]	 = 0;
		force[2*271+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*271+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*271+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*271+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*271+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(271+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*271+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*271] 		+= temp_force*dx;
			force[2*271+1]	+= temp_force*dy;
		}

		xi = position[2*272];
		yi = position[2*272+1];

		force[2*272]	 = 0;
		force[2*272+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*272+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*272+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*272+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*272+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(272+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*272+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*272] 		+= temp_force*dx;
			force[2*272+1]	+= temp_force*dy;
		}

		xi = position[2*273];
		yi = position[2*273+1];

		force[2*273]	 = 0;
		force[2*273+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*273+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*273+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*273+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*273+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(273+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*273+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*273] 		+= temp_force*dx;
			force[2*273+1]	+= temp_force*dy;
		}

		xi = position[2*274];
		yi = position[2*274+1];

		force[2*274]	 = 0;
		force[2*274+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*274+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*274+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*274+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*274+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(274+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*274+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*274] 		+= temp_force*dx;
			force[2*274+1]	+= temp_force*dy;
		}

		xi = position[2*275];
		yi = position[2*275+1];

		force[2*275]	 = 0;
		force[2*275+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*275+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*275+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*275+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*275+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(275+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*275+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*275] 		+= temp_force*dx;
			force[2*275+1]	+= temp_force*dy;
		}

		xi = position[2*276];
		yi = position[2*276+1];

		force[2*276]	 = 0;
		force[2*276+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*276+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*276+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*276+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*276+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(276+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*276+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*276] 		+= temp_force*dx;
			force[2*276+1]	+= temp_force*dy;
		}

		xi = position[2*277];
		yi = position[2*277+1];

		force[2*277]	 = 0;
		force[2*277+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*277+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*277+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*277+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*277+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(277+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*277+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*277] 		+= temp_force*dx;
			force[2*277+1]	+= temp_force*dy;
		}

		xi = position[2*278];
		yi = position[2*278+1];

		force[2*278]	 = 0;
		force[2*278+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*278+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*278+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*278+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*278+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(278+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*278+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*278] 		+= temp_force*dx;
			force[2*278+1]	+= temp_force*dy;
		}

		xi = position[2*279];
		yi = position[2*279+1];

		force[2*279]	 = 0;
		force[2*279+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*279+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*279+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*279+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*279+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(279+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*279+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*279] 		+= temp_force*dx;
			force[2*279+1]	+= temp_force*dy;
		}

		xi = position[2*280];
		yi = position[2*280+1];

		force[2*280]	 = 0;
		force[2*280+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*280+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*280+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*280+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*280+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(280+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*280+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*280] 		+= temp_force*dx;
			force[2*280+1]	+= temp_force*dy;
		}

		xi = position[2*281];
		yi = position[2*281+1];

		force[2*281]	 = 0;
		force[2*281+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*281+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*281+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*281+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*281+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(281+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*281+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*281] 		+= temp_force*dx;
			force[2*281+1]	+= temp_force*dy;
		}

		xi = position[2*282];
		yi = position[2*282+1];

		force[2*282]	 = 0;
		force[2*282+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*282+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*282+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*282+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*282+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(282+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*282+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*282] 		+= temp_force*dx;
			force[2*282+1]	+= temp_force*dy;
		}

		xi = position[2*283];
		yi = position[2*283+1];

		force[2*283]	 = 0;
		force[2*283+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*283+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*283+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*283+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*283+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(283+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*283+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*283] 		+= temp_force*dx;
			force[2*283+1]	+= temp_force*dy;
		}

		xi = position[2*284];
		yi = position[2*284+1];

		force[2*284]	 = 0;
		force[2*284+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*284+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*284+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*284+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*284+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(284+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*284+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*284] 		+= temp_force*dx;
			force[2*284+1]	+= temp_force*dy;
		}

		xi = position[2*285];
		yi = position[2*285+1];

		force[2*285]	 = 0;
		force[2*285+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*285+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*285+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*285+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*285+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(285+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*285+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*285] 		+= temp_force*dx;
			force[2*285+1]	+= temp_force*dy;
		}

		xi = position[2*286];
		yi = position[2*286+1];

		force[2*286]	 = 0;
		force[2*286+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*286+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*286+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*286+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*286+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(286+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*286+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*286] 		+= temp_force*dx;
			force[2*286+1]	+= temp_force*dy;
		}

		xi = position[2*287];
		yi = position[2*287+1];

		force[2*287]	 = 0;
		force[2*287+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*287+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*287+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*287+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*287+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(287+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*287+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*287] 		+= temp_force*dx;
			force[2*287+1]	+= temp_force*dy;
		}

		xi = position[2*288];
		yi = position[2*288+1];

		force[2*288]	 = 0;
		force[2*288+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*288+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*288+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*288+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*288+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(288+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*288+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*288] 		+= temp_force*dx;
			force[2*288+1]	+= temp_force*dy;
		}

		xi = position[2*289];
		yi = position[2*289+1];

		force[2*289]	 = 0;
		force[2*289+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*289+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*289+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*289+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*289+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(289+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*289+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*289] 		+= temp_force*dx;
			force[2*289+1]	+= temp_force*dy;
		}

		xi = position[2*290];
		yi = position[2*290+1];

		force[2*290]	 = 0;
		force[2*290+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*290+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*290+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*290+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*290+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(290+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*290+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*290] 		+= temp_force*dx;
			force[2*290+1]	+= temp_force*dy;
		}

		xi = position[2*291];
		yi = position[2*291+1];

		force[2*291]	 = 0;
		force[2*291+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*291+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*291+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*291+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*291+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(291+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*291+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*291] 		+= temp_force*dx;
			force[2*291+1]	+= temp_force*dy;
		}

		xi = position[2*292];
		yi = position[2*292+1];

		force[2*292]	 = 0;
		force[2*292+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*292+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*292+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*292+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*292+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(292+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*292+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*292] 		+= temp_force*dx;
			force[2*292+1]	+= temp_force*dy;
		}

		xi = position[2*293];
		yi = position[2*293+1];

		force[2*293]	 = 0;
		force[2*293+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*293+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*293+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*293+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*293+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(293+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*293+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*293] 		+= temp_force*dx;
			force[2*293+1]	+= temp_force*dy;
		}

		xi = position[2*294];
		yi = position[2*294+1];

		force[2*294]	 = 0;
		force[2*294+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*294+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*294+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*294+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*294+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(294+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*294+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*294] 		+= temp_force*dx;
			force[2*294+1]	+= temp_force*dy;
		}

		xi = position[2*295];
		yi = position[2*295+1];

		force[2*295]	 = 0;
		force[2*295+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*295+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*295+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*295+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*295+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(295+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*295+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*295] 		+= temp_force*dx;
			force[2*295+1]	+= temp_force*dy;
		}

		xi = position[2*296];
		yi = position[2*296+1];

		force[2*296]	 = 0;
		force[2*296+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*296+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*296+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*296+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*296+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(296+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*296+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*296] 		+= temp_force*dx;
			force[2*296+1]	+= temp_force*dy;
		}

		xi = position[2*297];
		yi = position[2*297+1];

		force[2*297]	 = 0;
		force[2*297+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*297+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*297+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*297+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*297+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(297+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*297+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*297] 		+= temp_force*dx;
			force[2*297+1]	+= temp_force*dy;
		}

		xi = position[2*298];
		yi = position[2*298+1];

		force[2*298]	 = 0;
		force[2*298+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*298+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*298+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*298+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*298+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(298+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*298+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*298] 		+= temp_force*dx;
			force[2*298+1]	+= temp_force*dy;
		}

		xi = position[2*299];
		yi = position[2*299+1];

		force[2*299]	 = 0;
		force[2*299+1] = 0;

		if (yi <= 1e-12) {
			dist_bottom = 1e-12;
			force[2*299+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}
		else if (yi <= dist_cutoff) {
			dist_bottom = yi;
			force[2*299+1] += 1.0*prefactor *(kappa/(dist_bottom) + 1.0/(dist_bottom*dist_bottom))*exp(-kappa*yi) - force_wall_cutoff;
		}

			if (L_y - yi <= 1e-12) {
				dist_top = 1e-12;
				force[2*299+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}
			else if (L_y - yi <= dist_cutoff) {
				dist_top = L_y-yi;
				force[2*299+1] -= 1.0*prefactor *(kappa/(dist_top) + 1.0/(dist_top*dist_top))*exp(-kappa*dist_top) - force_wall_cutoff;
			}

		iterate = verlet[N*(299+1)-1];

			for (int k=0; k<iterate; k+=1) {

			j  = verlet[N*299+k];
			xj = position[2*j];
			yj = position[2*j+1];

			m_j = (j%2)*m + (j+1)%2;

			dx = xi - xj;
			dy = yi - yj;

			dx 		-= dround(dx/L_x)*L_x;

			r_squared = dx*dx + dy*dy;
			temp_force = (m_i_B*m_j/sqrt(r_squared))*(3*Gamma_A/(r_squared*r_squared)-force_cutoff);

			force[2*299] 		+= temp_force*dx;
			force[2*299+1]	+= temp_force*dy;
		}

		pthread_barrier_wait(&barrier_internal);

		xi = position[2*263];
		yi = position[2*263+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*263]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*263+1] + weigh_brown_B * g2;

		position[2*263]		+= dx;
		position[2*263+1] 	+= dy;

		displacement[2*263]	+= dx;
		displacement[2*263+1] += dy;

		verlet_distance[2*263] 	+= (xi - position[2*263]);
		verlet_distance[2*263+1] 	+= (yi - position[2*263+1]);

		if ((temp = sqrt(verlet_distance[2*263]*verlet_distance[2*263] + verlet_distance[2*263+1]*verlet_distance[2*263+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*263] 	-= floor(position[2*263]/L_x)*L_x;

		xi = position[2*264];
		yi = position[2*264+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*264]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*264+1] + weigh_brown_B * g2;

		position[2*264]		+= dx;
		position[2*264+1] 	+= dy;

		displacement[2*264]	+= dx;
		displacement[2*264+1] += dy;

		verlet_distance[2*264] 	+= (xi - position[2*264]);
		verlet_distance[2*264+1] 	+= (yi - position[2*264+1]);

		if ((temp = sqrt(verlet_distance[2*264]*verlet_distance[2*264] + verlet_distance[2*264+1]*verlet_distance[2*264+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*264] 	-= floor(position[2*264]/L_x)*L_x;

		xi = position[2*265];
		yi = position[2*265+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*265]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*265+1] + weigh_brown_B * g2;

		position[2*265]		+= dx;
		position[2*265+1] 	+= dy;

		displacement[2*265]	+= dx;
		displacement[2*265+1] += dy;

		verlet_distance[2*265] 	+= (xi - position[2*265]);
		verlet_distance[2*265+1] 	+= (yi - position[2*265+1]);

		if ((temp = sqrt(verlet_distance[2*265]*verlet_distance[2*265] + verlet_distance[2*265+1]*verlet_distance[2*265+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*265] 	-= floor(position[2*265]/L_x)*L_x;

		xi = position[2*266];
		yi = position[2*266+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*266]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*266+1] + weigh_brown_B * g2;

		position[2*266]		+= dx;
		position[2*266+1] 	+= dy;

		displacement[2*266]	+= dx;
		displacement[2*266+1] += dy;

		verlet_distance[2*266] 	+= (xi - position[2*266]);
		verlet_distance[2*266+1] 	+= (yi - position[2*266+1]);

		if ((temp = sqrt(verlet_distance[2*266]*verlet_distance[2*266] + verlet_distance[2*266+1]*verlet_distance[2*266+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*266] 	-= floor(position[2*266]/L_x)*L_x;

		xi = position[2*267];
		yi = position[2*267+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*267]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*267+1] + weigh_brown_B * g2;

		position[2*267]		+= dx;
		position[2*267+1] 	+= dy;

		displacement[2*267]	+= dx;
		displacement[2*267+1] += dy;

		verlet_distance[2*267] 	+= (xi - position[2*267]);
		verlet_distance[2*267+1] 	+= (yi - position[2*267+1]);

		if ((temp = sqrt(verlet_distance[2*267]*verlet_distance[2*267] + verlet_distance[2*267+1]*verlet_distance[2*267+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*267] 	-= floor(position[2*267]/L_x)*L_x;

		xi = position[2*268];
		yi = position[2*268+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*268]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*268+1] + weigh_brown_B * g2;

		position[2*268]		+= dx;
		position[2*268+1] 	+= dy;

		displacement[2*268]	+= dx;
		displacement[2*268+1] += dy;

		verlet_distance[2*268] 	+= (xi - position[2*268]);
		verlet_distance[2*268+1] 	+= (yi - position[2*268+1]);

		if ((temp = sqrt(verlet_distance[2*268]*verlet_distance[2*268] + verlet_distance[2*268+1]*verlet_distance[2*268+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*268] 	-= floor(position[2*268]/L_x)*L_x;

		xi = position[2*269];
		yi = position[2*269+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*269]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*269+1] + weigh_brown_B * g2;

		position[2*269]		+= dx;
		position[2*269+1] 	+= dy;

		displacement[2*269]	+= dx;
		displacement[2*269+1] += dy;

		verlet_distance[2*269] 	+= (xi - position[2*269]);
		verlet_distance[2*269+1] 	+= (yi - position[2*269+1]);

		if ((temp = sqrt(verlet_distance[2*269]*verlet_distance[2*269] + verlet_distance[2*269+1]*verlet_distance[2*269+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*269] 	-= floor(position[2*269]/L_x)*L_x;

		xi = position[2*270];
		yi = position[2*270+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*270]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*270+1] + weigh_brown_B * g2;

		position[2*270]		+= dx;
		position[2*270+1] 	+= dy;

		displacement[2*270]	+= dx;
		displacement[2*270+1] += dy;

		verlet_distance[2*270] 	+= (xi - position[2*270]);
		verlet_distance[2*270+1] 	+= (yi - position[2*270+1]);

		if ((temp = sqrt(verlet_distance[2*270]*verlet_distance[2*270] + verlet_distance[2*270+1]*verlet_distance[2*270+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*270] 	-= floor(position[2*270]/L_x)*L_x;

		xi = position[2*271];
		yi = position[2*271+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*271]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*271+1] + weigh_brown_B * g2;

		position[2*271]		+= dx;
		position[2*271+1] 	+= dy;

		displacement[2*271]	+= dx;
		displacement[2*271+1] += dy;

		verlet_distance[2*271] 	+= (xi - position[2*271]);
		verlet_distance[2*271+1] 	+= (yi - position[2*271+1]);

		if ((temp = sqrt(verlet_distance[2*271]*verlet_distance[2*271] + verlet_distance[2*271+1]*verlet_distance[2*271+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*271] 	-= floor(position[2*271]/L_x)*L_x;

		xi = position[2*272];
		yi = position[2*272+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*272]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*272+1] + weigh_brown_B * g2;

		position[2*272]		+= dx;
		position[2*272+1] 	+= dy;

		displacement[2*272]	+= dx;
		displacement[2*272+1] += dy;

		verlet_distance[2*272] 	+= (xi - position[2*272]);
		verlet_distance[2*272+1] 	+= (yi - position[2*272+1]);

		if ((temp = sqrt(verlet_distance[2*272]*verlet_distance[2*272] + verlet_distance[2*272+1]*verlet_distance[2*272+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*272] 	-= floor(position[2*272]/L_x)*L_x;

		xi = position[2*273];
		yi = position[2*273+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*273]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*273+1] + weigh_brown_B * g2;

		position[2*273]		+= dx;
		position[2*273+1] 	+= dy;

		displacement[2*273]	+= dx;
		displacement[2*273+1] += dy;

		verlet_distance[2*273] 	+= (xi - position[2*273]);
		verlet_distance[2*273+1] 	+= (yi - position[2*273+1]);

		if ((temp = sqrt(verlet_distance[2*273]*verlet_distance[2*273] + verlet_distance[2*273+1]*verlet_distance[2*273+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*273] 	-= floor(position[2*273]/L_x)*L_x;

		xi = position[2*274];
		yi = position[2*274+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*274]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*274+1] + weigh_brown_B * g2;

		position[2*274]		+= dx;
		position[2*274+1] 	+= dy;

		displacement[2*274]	+= dx;
		displacement[2*274+1] += dy;

		verlet_distance[2*274] 	+= (xi - position[2*274]);
		verlet_distance[2*274+1] 	+= (yi - position[2*274+1]);

		if ((temp = sqrt(verlet_distance[2*274]*verlet_distance[2*274] + verlet_distance[2*274+1]*verlet_distance[2*274+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*274] 	-= floor(position[2*274]/L_x)*L_x;

		xi = position[2*275];
		yi = position[2*275+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*275]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*275+1] + weigh_brown_B * g2;

		position[2*275]		+= dx;
		position[2*275+1] 	+= dy;

		displacement[2*275]	+= dx;
		displacement[2*275+1] += dy;

		verlet_distance[2*275] 	+= (xi - position[2*275]);
		verlet_distance[2*275+1] 	+= (yi - position[2*275+1]);

		if ((temp = sqrt(verlet_distance[2*275]*verlet_distance[2*275] + verlet_distance[2*275+1]*verlet_distance[2*275+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*275] 	-= floor(position[2*275]/L_x)*L_x;

		xi = position[2*276];
		yi = position[2*276+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*276]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*276+1] + weigh_brown_B * g2;

		position[2*276]		+= dx;
		position[2*276+1] 	+= dy;

		displacement[2*276]	+= dx;
		displacement[2*276+1] += dy;

		verlet_distance[2*276] 	+= (xi - position[2*276]);
		verlet_distance[2*276+1] 	+= (yi - position[2*276+1]);

		if ((temp = sqrt(verlet_distance[2*276]*verlet_distance[2*276] + verlet_distance[2*276+1]*verlet_distance[2*276+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*276] 	-= floor(position[2*276]/L_x)*L_x;

		xi = position[2*277];
		yi = position[2*277+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*277]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*277+1] + weigh_brown_B * g2;

		position[2*277]		+= dx;
		position[2*277+1] 	+= dy;

		displacement[2*277]	+= dx;
		displacement[2*277+1] += dy;

		verlet_distance[2*277] 	+= (xi - position[2*277]);
		verlet_distance[2*277+1] 	+= (yi - position[2*277+1]);

		if ((temp = sqrt(verlet_distance[2*277]*verlet_distance[2*277] + verlet_distance[2*277+1]*verlet_distance[2*277+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*277] 	-= floor(position[2*277]/L_x)*L_x;

		xi = position[2*278];
		yi = position[2*278+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*278]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*278+1] + weigh_brown_B * g2;

		position[2*278]		+= dx;
		position[2*278+1] 	+= dy;

		displacement[2*278]	+= dx;
		displacement[2*278+1] += dy;

		verlet_distance[2*278] 	+= (xi - position[2*278]);
		verlet_distance[2*278+1] 	+= (yi - position[2*278+1]);

		if ((temp = sqrt(verlet_distance[2*278]*verlet_distance[2*278] + verlet_distance[2*278+1]*verlet_distance[2*278+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*278] 	-= floor(position[2*278]/L_x)*L_x;

		xi = position[2*279];
		yi = position[2*279+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*279]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*279+1] + weigh_brown_B * g2;

		position[2*279]		+= dx;
		position[2*279+1] 	+= dy;

		displacement[2*279]	+= dx;
		displacement[2*279+1] += dy;

		verlet_distance[2*279] 	+= (xi - position[2*279]);
		verlet_distance[2*279+1] 	+= (yi - position[2*279+1]);

		if ((temp = sqrt(verlet_distance[2*279]*verlet_distance[2*279] + verlet_distance[2*279+1]*verlet_distance[2*279+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*279] 	-= floor(position[2*279]/L_x)*L_x;

		xi = position[2*280];
		yi = position[2*280+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*280]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*280+1] + weigh_brown_B * g2;

		position[2*280]		+= dx;
		position[2*280+1] 	+= dy;

		displacement[2*280]	+= dx;
		displacement[2*280+1] += dy;

		verlet_distance[2*280] 	+= (xi - position[2*280]);
		verlet_distance[2*280+1] 	+= (yi - position[2*280+1]);

		if ((temp = sqrt(verlet_distance[2*280]*verlet_distance[2*280] + verlet_distance[2*280+1]*verlet_distance[2*280+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*280] 	-= floor(position[2*280]/L_x)*L_x;

		xi = position[2*281];
		yi = position[2*281+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*281]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*281+1] + weigh_brown_B * g2;

		position[2*281]		+= dx;
		position[2*281+1] 	+= dy;

		displacement[2*281]	+= dx;
		displacement[2*281+1] += dy;

		verlet_distance[2*281] 	+= (xi - position[2*281]);
		verlet_distance[2*281+1] 	+= (yi - position[2*281+1]);

		if ((temp = sqrt(verlet_distance[2*281]*verlet_distance[2*281] + verlet_distance[2*281+1]*verlet_distance[2*281+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*281] 	-= floor(position[2*281]/L_x)*L_x;

		xi = position[2*282];
		yi = position[2*282+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*282]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*282+1] + weigh_brown_B * g2;

		position[2*282]		+= dx;
		position[2*282+1] 	+= dy;

		displacement[2*282]	+= dx;
		displacement[2*282+1] += dy;

		verlet_distance[2*282] 	+= (xi - position[2*282]);
		verlet_distance[2*282+1] 	+= (yi - position[2*282+1]);

		if ((temp = sqrt(verlet_distance[2*282]*verlet_distance[2*282] + verlet_distance[2*282+1]*verlet_distance[2*282+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*282] 	-= floor(position[2*282]/L_x)*L_x;

		xi = position[2*283];
		yi = position[2*283+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*283]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*283+1] + weigh_brown_B * g2;

		position[2*283]		+= dx;
		position[2*283+1] 	+= dy;

		displacement[2*283]	+= dx;
		displacement[2*283+1] += dy;

		verlet_distance[2*283] 	+= (xi - position[2*283]);
		verlet_distance[2*283+1] 	+= (yi - position[2*283+1]);

		if ((temp = sqrt(verlet_distance[2*283]*verlet_distance[2*283] + verlet_distance[2*283+1]*verlet_distance[2*283+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*283] 	-= floor(position[2*283]/L_x)*L_x;

		xi = position[2*284];
		yi = position[2*284+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*284]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*284+1] + weigh_brown_B * g2;

		position[2*284]		+= dx;
		position[2*284+1] 	+= dy;

		displacement[2*284]	+= dx;
		displacement[2*284+1] += dy;

		verlet_distance[2*284] 	+= (xi - position[2*284]);
		verlet_distance[2*284+1] 	+= (yi - position[2*284+1]);

		if ((temp = sqrt(verlet_distance[2*284]*verlet_distance[2*284] + verlet_distance[2*284+1]*verlet_distance[2*284+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*284] 	-= floor(position[2*284]/L_x)*L_x;

		xi = position[2*285];
		yi = position[2*285+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*285]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*285+1] + weigh_brown_B * g2;

		position[2*285]		+= dx;
		position[2*285+1] 	+= dy;

		displacement[2*285]	+= dx;
		displacement[2*285+1] += dy;

		verlet_distance[2*285] 	+= (xi - position[2*285]);
		verlet_distance[2*285+1] 	+= (yi - position[2*285+1]);

		if ((temp = sqrt(verlet_distance[2*285]*verlet_distance[2*285] + verlet_distance[2*285+1]*verlet_distance[2*285+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*285] 	-= floor(position[2*285]/L_x)*L_x;

		xi = position[2*286];
		yi = position[2*286+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*286]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*286+1] + weigh_brown_B * g2;

		position[2*286]		+= dx;
		position[2*286+1] 	+= dy;

		displacement[2*286]	+= dx;
		displacement[2*286+1] += dy;

		verlet_distance[2*286] 	+= (xi - position[2*286]);
		verlet_distance[2*286+1] 	+= (yi - position[2*286+1]);

		if ((temp = sqrt(verlet_distance[2*286]*verlet_distance[2*286] + verlet_distance[2*286+1]*verlet_distance[2*286+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*286] 	-= floor(position[2*286]/L_x)*L_x;

		xi = position[2*287];
		yi = position[2*287+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*287]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*287+1] + weigh_brown_B * g2;

		position[2*287]		+= dx;
		position[2*287+1] 	+= dy;

		displacement[2*287]	+= dx;
		displacement[2*287+1] += dy;

		verlet_distance[2*287] 	+= (xi - position[2*287]);
		verlet_distance[2*287+1] 	+= (yi - position[2*287+1]);

		if ((temp = sqrt(verlet_distance[2*287]*verlet_distance[2*287] + verlet_distance[2*287+1]*verlet_distance[2*287+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*287] 	-= floor(position[2*287]/L_x)*L_x;

		xi = position[2*288];
		yi = position[2*288+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*288]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*288+1] + weigh_brown_B * g2;

		position[2*288]		+= dx;
		position[2*288+1] 	+= dy;

		displacement[2*288]	+= dx;
		displacement[2*288+1] += dy;

		verlet_distance[2*288] 	+= (xi - position[2*288]);
		verlet_distance[2*288+1] 	+= (yi - position[2*288+1]);

		if ((temp = sqrt(verlet_distance[2*288]*verlet_distance[2*288] + verlet_distance[2*288+1]*verlet_distance[2*288+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*288] 	-= floor(position[2*288]/L_x)*L_x;

		xi = position[2*289];
		yi = position[2*289+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*289]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*289+1] + weigh_brown_B * g2;

		position[2*289]		+= dx;
		position[2*289+1] 	+= dy;

		displacement[2*289]	+= dx;
		displacement[2*289+1] += dy;

		verlet_distance[2*289] 	+= (xi - position[2*289]);
		verlet_distance[2*289+1] 	+= (yi - position[2*289+1]);

		if ((temp = sqrt(verlet_distance[2*289]*verlet_distance[2*289] + verlet_distance[2*289+1]*verlet_distance[2*289+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*289] 	-= floor(position[2*289]/L_x)*L_x;

		xi = position[2*290];
		yi = position[2*290+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*290]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*290+1] + weigh_brown_B * g2;

		position[2*290]		+= dx;
		position[2*290+1] 	+= dy;

		displacement[2*290]	+= dx;
		displacement[2*290+1] += dy;

		verlet_distance[2*290] 	+= (xi - position[2*290]);
		verlet_distance[2*290+1] 	+= (yi - position[2*290+1]);

		if ((temp = sqrt(verlet_distance[2*290]*verlet_distance[2*290] + verlet_distance[2*290+1]*verlet_distance[2*290+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*290] 	-= floor(position[2*290]/L_x)*L_x;

		xi = position[2*291];
		yi = position[2*291+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*291]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*291+1] + weigh_brown_B * g2;

		position[2*291]		+= dx;
		position[2*291+1] 	+= dy;

		displacement[2*291]	+= dx;
		displacement[2*291+1] += dy;

		verlet_distance[2*291] 	+= (xi - position[2*291]);
		verlet_distance[2*291+1] 	+= (yi - position[2*291+1]);

		if ((temp = sqrt(verlet_distance[2*291]*verlet_distance[2*291] + verlet_distance[2*291+1]*verlet_distance[2*291+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*291] 	-= floor(position[2*291]/L_x)*L_x;

		xi = position[2*292];
		yi = position[2*292+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*292]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*292+1] + weigh_brown_B * g2;

		position[2*292]		+= dx;
		position[2*292+1] 	+= dy;

		displacement[2*292]	+= dx;
		displacement[2*292+1] += dy;

		verlet_distance[2*292] 	+= (xi - position[2*292]);
		verlet_distance[2*292+1] 	+= (yi - position[2*292+1]);

		if ((temp = sqrt(verlet_distance[2*292]*verlet_distance[2*292] + verlet_distance[2*292+1]*verlet_distance[2*292+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*292] 	-= floor(position[2*292]/L_x)*L_x;

		xi = position[2*293];
		yi = position[2*293+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*293]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*293+1] + weigh_brown_B * g2;

		position[2*293]		+= dx;
		position[2*293+1] 	+= dy;

		displacement[2*293]	+= dx;
		displacement[2*293+1] += dy;

		verlet_distance[2*293] 	+= (xi - position[2*293]);
		verlet_distance[2*293+1] 	+= (yi - position[2*293+1]);

		if ((temp = sqrt(verlet_distance[2*293]*verlet_distance[2*293] + verlet_distance[2*293+1]*verlet_distance[2*293+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*293] 	-= floor(position[2*293]/L_x)*L_x;

		xi = position[2*294];
		yi = position[2*294+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*294]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*294+1] + weigh_brown_B * g2;

		position[2*294]		+= dx;
		position[2*294+1] 	+= dy;

		displacement[2*294]	+= dx;
		displacement[2*294+1] += dy;

		verlet_distance[2*294] 	+= (xi - position[2*294]);
		verlet_distance[2*294+1] 	+= (yi - position[2*294+1]);

		if ((temp = sqrt(verlet_distance[2*294]*verlet_distance[2*294] + verlet_distance[2*294+1]*verlet_distance[2*294+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*294] 	-= floor(position[2*294]/L_x)*L_x;

		xi = position[2*295];
		yi = position[2*295+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*295]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*295+1] + weigh_brown_B * g2;

		position[2*295]		+= dx;
		position[2*295+1] 	+= dy;

		displacement[2*295]	+= dx;
		displacement[2*295+1] += dy;

		verlet_distance[2*295] 	+= (xi - position[2*295]);
		verlet_distance[2*295+1] 	+= (yi - position[2*295+1]);

		if ((temp = sqrt(verlet_distance[2*295]*verlet_distance[2*295] + verlet_distance[2*295+1]*verlet_distance[2*295+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*295] 	-= floor(position[2*295]/L_x)*L_x;

		xi = position[2*296];
		yi = position[2*296+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*296]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*296+1] + weigh_brown_B * g2;

		position[2*296]		+= dx;
		position[2*296+1] 	+= dy;

		displacement[2*296]	+= dx;
		displacement[2*296+1] += dy;

		verlet_distance[2*296] 	+= (xi - position[2*296]);
		verlet_distance[2*296+1] 	+= (yi - position[2*296+1]);

		if ((temp = sqrt(verlet_distance[2*296]*verlet_distance[2*296] + verlet_distance[2*296+1]*verlet_distance[2*296+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*296] 	-= floor(position[2*296]/L_x)*L_x;

		xi = position[2*297];
		yi = position[2*297+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*297]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*297+1] + weigh_brown_B * g2;

		position[2*297]		+= dx;
		position[2*297+1] 	+= dy;

		displacement[2*297]	+= dx;
		displacement[2*297+1] += dy;

		verlet_distance[2*297] 	+= (xi - position[2*297]);
		verlet_distance[2*297+1] 	+= (yi - position[2*297+1]);

		if ((temp = sqrt(verlet_distance[2*297]*verlet_distance[2*297] + verlet_distance[2*297+1]*verlet_distance[2*297+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*297] 	-= floor(position[2*297]/L_x)*L_x;

		xi = position[2*298];
		yi = position[2*298+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*298]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*298+1] + weigh_brown_B * g2;

		position[2*298]		+= dx;
		position[2*298+1] 	+= dy;

		displacement[2*298]	+= dx;
		displacement[2*298+1] += dy;

		verlet_distance[2*298] 	+= (xi - position[2*298]);
		verlet_distance[2*298+1] 	+= (yi - position[2*298+1]);

		if ((temp = sqrt(verlet_distance[2*298]*verlet_distance[2*298] + verlet_distance[2*298+1]*verlet_distance[2*298+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*298] 	-= floor(position[2*298]/L_x)*L_x;

		xi = position[2*299];
		yi = position[2*299+1];

			do {
				u1 = rand()/(double)RAND_MAX;
			} while (u1 == 0.0);
			do {
				u2 = rand()/(double)RAND_MAX;
			} while (u2 == 0.0);

			temp = sqrt(-2*log(u1));
			g1 = temp*cos(TWO_PI*u2);
			g2 = temp*sin(TWO_PI*u2);

		dx = D_kT_B*delta_t*force[2*299]   + weigh_brown_B * g1;
		dy = D_kT_B*delta_t*force[2*299+1] + weigh_brown_B * g2;

		position[2*299]		+= dx;
		position[2*299+1] 	+= dy;

		displacement[2*299]	+= dx;
		displacement[2*299+1] += dy;

		verlet_distance[2*299] 	+= (xi - position[2*299]);
		verlet_distance[2*299+1] 	+= (yi - position[2*299+1]);

		if ((temp = sqrt(verlet_distance[2*299]*verlet_distance[2*299] + verlet_distance[2*299+1]*verlet_distance[2*299+1])) > verlet_max[2*(*no)]) {
			verlet_max[2*(*no)+1] = verlet_max[2*(*no)];
			verlet_max[2*(*no)] = temp;
		} else if (temp > verlet_max[2*(*no)+1]) {
			verlet_max[2*(*no)+1] = temp;
		}

		position[2*299] 	-= floor(position[2*299]/L_x)*L_x;

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

	ret_thread 	= pthread_create(&(threads[0]), NULL, (void*)&iteration_0, &(numbers[0]));

	if (ret_thread != 0) {
		fprintf(stderr, "Thread %d could not be created. \n", 0);
		exit(EXIT_FAILURE);
	}
	ret_thread 	= pthread_create(&(threads[1]), NULL, (void*)&iteration_1, &(numbers[1]));

	if (ret_thread != 0) {
		fprintf(stderr, "Thread %d could not be created. \n", 1);
		exit(EXIT_FAILURE);
	}
	ret_thread 	= pthread_create(&(threads[2]), NULL, (void*)&iteration_2, &(numbers[2]));

	if (ret_thread != 0) {
		fprintf(stderr, "Thread %d could not be created. \n", 2);
		exit(EXIT_FAILURE);
	}
	ret_thread 	= pthread_create(&(threads[3]), NULL, (void*)&iteration_3, &(numbers[3]));

	if (ret_thread != 0) {
		fprintf(stderr, "Thread %d could not be created. \n", 3);
		exit(EXIT_FAILURE);
	}
	ret_thread 	= pthread_create(&(threads[4]), NULL, (void*)&iteration_4, &(numbers[4]));

	if (ret_thread != 0) {
		fprintf(stderr, "Thread %d could not be created. \n", 4);
		exit(EXIT_FAILURE);
	}
	ret_thread 	= pthread_create(&(threads[5]), NULL, (void*)&iteration_5, &(numbers[5]));

	if (ret_thread != 0) {
		fprintf(stderr, "Thread %d could not be created. \n", 5);
		exit(EXIT_FAILURE);
	}
	ret_thread 	= pthread_create(&(threads[6]), NULL, (void*)&iteration_6, &(numbers[6]));

	if (ret_thread != 0) {
		fprintf(stderr, "Thread %d could not be created. \n", 6);
		exit(EXIT_FAILURE);
	}
	ret_thread 	= pthread_create(&(threads[7]), NULL, (void*)&iteration_7, &(numbers[7]));

	if (ret_thread != 0) {
		fprintf(stderr, "Thread %d could not be created. \n", 7);
		exit(EXIT_FAILURE);
	}

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
