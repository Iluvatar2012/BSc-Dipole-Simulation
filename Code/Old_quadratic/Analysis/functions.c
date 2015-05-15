#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "structs.h"
#include "functions.h"

// basic variables
static int		N;

static double* 	config;
static double* 	psi4;
static double* 	psi6;
static double* 	laning;

/*----------------------------------------------------------------------------------------------------------------------------*/
void init (struct variables* var) {
	config 	= var->positions;
	N 		= var->N;
	psi4 	= var->psi4;
	psi6 	= var->psi6;
	laning 	= var->laning;
}


/*----------------------------------------------------------------------------------------------------------------------------*/
double dabs(double in) {

	// function for determining double precision absolute values
	if (in < 0)
		return -in;
	else
		return in;
}

/*----------------------------------------------------------------------------------------------------------------------------*/
void compute_psi4 (int l) {

	// basic variables
	double 	dx, dy;
	double 	r_sq;
	double 	theta;
	double	real, img, abs;

	double 	dist[4] = {1E6, 1E6, 1E6, 1E6};
	int		next[4];

	int		length  = 1;

	// iterate over all particles
	for (int i=1; i<N; i++) {

		// reinitiate both lists
		for (int k=0; k<4; k++) {
			dist[k] = 1E6;
			next[k] = N;
		}

		// get the four shortest distances to particle i
		for (int j=0; j<N; j++) {

			// ignore i=j
			if (i == j) 
				continue;

			// compute squared distance between particles i and j
			dx   = config[2*N*l+2*i]   - config[2*N*l+2*j];
			dy   = config[2*N*l+2*i+1] - config[2*N*l+2*j+1];

			r_sq = dx*dx + dy*dy;

			// continue if distance is too great
			if (r_sq > 5)
				continue;

			// order the distance and particle index within their respective lists
			if (r_sq < dist[3]) {
				dist[3] = r_sq;
				next[3] = j;

				bubble_sort(dist, next, 4);
			}
		}

		// calculate psi 4 value
		psi4[N*l+i] = psi_n(4, i, l, next);
	}
}



/*----------------------------------------------------------------------------------------------------------------------------*/
void compute_psi6(int l) {
		// basic variables
	double 	dx, dy;
	double 	r_sq;
	double 	theta;
	double	real, img, abs;

	double 	dist[6] = {1E6, 1E6, 1E6, 1E6, 1E6, 1E6};
	int		next[6];

	int		length  = 1;

	// iterate over all particles
	for (int i=0; i<N; i++) {

		// reinitiate both lists
		for (int k=0; k<6; k++) {
			dist[k] = 1E6;
			next[k] = N;
		}

		// get the four shortest distances to particle i
		for (int j=0; j<N; j++) {

			// ignore i=j
			if (i == j) 
				continue;

			// compute squared distance between particles i and j
			dx   = config[2*N*l+2*i]   - config[2*N*l+2*j];
			dy   = config[2*N*l+2*i+1] - config[2*N*l+2*j+1];

			r_sq = dx*dx + dy*dy;

			// continue if distance is too great
			if (r_sq > 5)
				continue;

			// order the distance and particle index within their respective lists
			if (r_sq < dist[5]) {
				dist[5] = r_sq;
				next[5] = j;

				bubble_sort(dist, next, 6);
			}
		}

		// calculate psi 6 value
		psi6[N*l+i] = psi_n (6, i, l, next);
	}
}


/*----------------------------------------------------------------------------------------------------------------------------*/
void compute_laning(int step) {
	
	// basic variables and boxlength
	double l, l1;
	double L = sqrt(N/2.);

	// iterate over all particles 
	for (int i=0; i<N; i++) {
		// reset the default length to one box length
		l = L;

		for (int j=0; j<N; j++) {
			// skip the particle itself and all of the same species
			if (i==j || i%2 == j%2)
				continue;

			// check whether the particles are within one lane, this is a lane of 0.25 around the particle
			if (dabs(config[2*N*step+2*i+1]-config[2*N*step+2*j+1]) < 1/8.) {

				// check whether particle j is to the right of particle i, if not, consider periodic picture
				if (config[2*N*step + 2*j] < config[2*N*step + 2*i])
					config[2*N*step + 2*j] += L;

				// store the x-distance between the two particles and the periodic pictures
				l1 = dabs(config[2*N*step+2*i]-config[2*N*step+2*j]);

				// check if a shorter lane has been found
				if (l1 < l)
					l = l1;
			}
		}

		// store the value in the laning array
		laning[N*step+i] = l;
	}
}


/*----------------------------------------------------------------------------------------------------------------------------*/
void bubble_sort (double* dist, int* next, int n) {

	// basic variables
	double 	temp;
	int		j;

	// variable for checking whether we are done
	int switched;

	// iterate until we are done
	do {
		// reset checking variable
		switched = 0;

		// iterate over the entire list
		for (int i=0; i<(n-1); i++) {

			// check whether an element is larger than it's follower
			if (dist[i] > dist[i+1]) {
				// switch the elements in both lists
				temp 	= dist[i+1];
				j		= next[i+1];

				dist[i+1] = dist[i];
				next[i+1] = next[i];

				dist[i] = temp;
				next[i] = j;

				// we need at least one more run
				switched = 1;
			}
		}

	} while (switched == 1);

}


/*----------------------------------------------------------------------------------------------------------------------------*/
double psi_n (int n, int i, int j, int* next) {

	// basic variables
	double dx, dy;
	double theta;

	// initiate real and imaginary part
	double real = 0;
	double img  = 0;

	// iterate over the saved particles
	for (int k=0; k<n; k++) {
			dx = config[2*N*j+2*next[k]]   - config[2*N*j+2*i];
			dy = config[2*N*j+2*next[k]+1] - config[2*N*j+2*i+1];

			theta	 = atan(dy/dx);

			if (dx < 0)
				theta += PI;

			real	+= cos(n*theta);
			img 	+= sin(n*theta);
	}

	// normalize real and imaginary part and calculate absolute value
	real /= n;
	img  /= n;

	return sqrt(real*real + img*img);
}
