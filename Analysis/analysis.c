#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdf5_read.h"
#include "picture.h"
#include "sort.h"

#define	cutoff_psi_4 0.9;
#define	cutoff_psi_6 0.9;

// basic variables
static double* 	config;
static double* 	psi4;
static double* 	psi6;

static int		N;



/*----------------------------------------------------------------------------------------------------------------------------*/
int compute_psi4 () {

	// basic variables
	double 	dx, dy;
	double 	r_sq;
	double 	theta;
	double	real, img, abs;

	double 	dist[4] = {1E6, 1E6, 1E6, 1E6};
	int		next[4];

	// allocate the necessary memory, check if successful
	psi4 = malloc(N*sizeof(double));

	if (psi4 == NULL) {
		return EXIT_FAILURE;
	}

	// iterate over all particles
	for (int i=0; i<N; i++) {

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
			dx   = config[2*i]   - config[2*j];
			dy   = config[2*i+1] - config[2*j+1];

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

		// initiate real and imaginary part
		real = 0;
		img  = 0;

		// iterate over the four saved particles
		for (int k=0; k<4; k++) {
			dx = config[2*next[k]] 	 - config[2*i];
			dy = config[2*next[k]+1] - config[2*i+1];

			theta	 = atan(dy/dx);
			real	+= cos(4*theta);
			img 	+= sin(4*theta);
		}

		// normalize real and imaginary part and calculate absolute value
		real /= 4.0;
		img  /= 4.0;

		psi4[i] = sqrt(real*real + img*img);
	}

	// return to caller
	return EXIT_SUCCESS;
}



/*----------------------------------------------------------------------------------------------------------------------------*/
int compute_psi6() {
		// basic variables
	double 	dx, dy;
	double 	r_sq;
	double 	theta;
	double	real, img, abs;

	double 	dist[6] = {1E6, 1E6, 1E6, 1E6, 1E6, 1E6};
	int		next[6];

	// allocate the necessary memory, check if successful
	psi6 = malloc(N*sizeof(double));

	if (psi6 == NULL) {
		return EXIT_FAILURE;
	}

	// iterate over all particles
	for (int i=0; i<N; i++) {

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
			dx   = config[2*i]   - config[2*j];
			dy   = config[2*i+1] - config[2*j+1];

			r_sq = dx*dx + dy*dy;

			// continue if distance is too great
			if (r_sq > 9)
				continue;

			// order the distance and particle index within their respective lists
			if (r_sq < dist[5]) {
				dist[5] = r_sq;
				next[5] = j;

				bubble_sort(dist, next, 6);
			}
		}

		// initiate real and imaginary part
		real = 0;
		img  = 0;

		// iterate over the four saved particles
		for (int k=0; k<6; k++) {
			dx = config[2*next[k]] 	 - config[2*i];
			dy = config[2*next[k]+1] - config[2*i+1];

			theta	 = atan(dy/dx);
			real	+= cos(6*theta);
			img 	+= sin(6*theta);
		}

		// normalize real and imaginary part and calculate absolute value
		real /= 6.0;
		img  /= 6.0;

		psi6[i] = sqrt(real*real + img*img);
	}

	// return to caller
	return EXIT_SUCCESS;
}



/*----------------------------------------------------------------------------------------------------------------------------*/
int main (int argcount, char** argvektor) {

	// basic variables
	int 	check;
	char* 	dot;
	size_t 	length;

	char	file[1024];

	// check whether user provided a filename and timestep to compute
	if (argcount != 3) {
		fprintf(stderr, "Please provide a filename (1) and timestep (2) (number, \"first\" or \"last\") to read.\n"
						"Program will now terminate. \n");
		return EXIT_SUCCESS;
	}

	// read the expected configuration from file, check whether successful
	config = read_config(argvektor[1], argvektor[2], &N);

	if (config == NULL)
		return EXIT_FAILURE;

	// compute psi-4 and psi-6 and check if successful
	if ((check = compute_psi4()) != EXIT_SUCCESS) {
		fprintf(stderr, "Computing Psi 4 failed, process will terminate. \n");
		return EXIT_FAILURE;
	}

	if ((check = compute_psi6()) != EXIT_SUCCESS) {
		fprintf(stderr, "Computing Psi 6 failed, process will terminate. \n");
		return EXIT_FAILURE;
	}

	// compute the filename and produce a picture, check whether successful
	strncpy(file, argvektor[1], 1024);

	dot		= strrchr(file, '.');
	*dot	= '\0';

	length 	= strlen(argvektor[1]);
	strncat(file, "__analysis.bmp", 1024-length);

	if ((check = graphicOutput (file, config, psi4, psi6, N)) != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

	// return to caller
	return EXIT_SUCCESS;
}