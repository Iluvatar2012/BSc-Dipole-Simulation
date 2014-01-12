#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdf5_read.h"
#include "picture.h"
#include "assist.h"

#define	cutoff_psi_4 0.9
#define	cutoff_psi_6 0.9

// basic variables
static double* 	config;
static int* 	psi4;
static int* 	psi6;

static int		N;
static int 		counter_4 = 0;
static int 		counter_6 = 0;



/*----------------------------------------------------------------------------------------------------------------------------*/
int compute_psi4 () {

	// basic variables
	double 	dx, dy;
	double 	r_sq;
	double 	theta;
	double	real, img, abs;

	double 	dist[4] = {1E6, 1E6, 1E6, 1E6};
	int		next[4];

	int		length  = 1;

	// allocate the necessary memory, check if successful
	psi4 = malloc(5*sizeof(int));

	if (psi4 == NULL) {
		return EXIT_FAILURE;
	}

	// iterate over all particles
	for (int i=1; i<N; i+=2) {

		// reinitiate both lists
		for (int k=0; k<4; k++) {
			dist[k] = 1E6;
			next[k] = N;
		}

		// get the four shortest distances to particle i
		for (int j=0; j<N; j+=2) {

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

		// calculate psi 4 value and check whether it shall be accepted
		if (psi_n(4, i, config, next) > cutoff_psi_4) {

			// check whether there is enough space in the list
			if (counter_4 >= length) {

				// double the length of the array
				length *= 2;

				// list is too short, reallocate memory
				psi4 = realloc(psi4, length*5*sizeof(int));

			}

			// save the particle i
			psi4[5*counter_4] = i;

			// save nearest neighbors
			for (int k=0; k<4; k++) {
				psi4[5*counter_4 + (k+1)] = next[k];
			}

			// increment counter
			counter_4++;
		}
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

	int		length  = 1;

	// allocate the necessary memory, check if successful
	psi6 = malloc(7*sizeof(int));

	if (psi6 == NULL) {
		return EXIT_FAILURE;
	}

	// iterate over all particles
	for (int i=0; i<N; i+=2) {

		// reinitiate both lists
		for (int k=0; k<6; k++) {
			dist[k] = 1E6;
			next[k] = N;
		}

		// get the four shortest distances to particle i
		for (int j=0; j<N; j+=2) {

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
			if (r_sq < dist[5]) {
				dist[5] = r_sq;
				next[5] = j;

				bubble_sort(dist, next, 6);
			}
		}

		// calculate psi 6 value and check whether it shall be accepted
		if (psi_n (6, i, config, next) > cutoff_psi_6) {

			// check whether there is enough space in the list
			if (counter_6 >= length) {

				// double the length of the array
				length *= 2;

				// list is too short, reallocate memory
				psi6 = realloc(psi6, length*7*sizeof(int));

			}

			// save the particle i
			psi6[7*counter_6] = i;

			// save nearest neighbors
			for (int k=0; k<6; k++) {
				psi6[7*counter_6 + (k+1)] = next[k];
			}

			// increment counter
			counter_6++;
		}
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

	if ((check = graphicOutput (file, config, psi4, psi6, N, counter_4, counter_6)) != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

	// return to caller
	return EXIT_SUCCESS;
}