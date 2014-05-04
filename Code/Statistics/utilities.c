#include <stdlib.h>
#include <stdio.h>

#include <math.h>

#include "structs.h"

// static variables for all parameters
static int 		N;
static int 		steps;

static double* 	config;
static double* 	disp;

static double* 	laning;
static double* 	psi4;
static double* 	psi6;

static double*	mean_lane;
static double*	mean_psi4;
static double* 	mean_psi6;

static double*	stddev_lane;
static double* 	stddev_psi4;
static double* 	stddev_psi6;

static double*	max_lane;
static double* 	y_lane;


/*----------------------------------------------------------------------------------------------------------------------------*/
double dabs(double in) {

	// function for determining double precision absolute values
	if (in < 0)
		return -in;
	else
		return in;
}


/*----------------------------------------------------------------------------------------------------------------------------*/
void stat_laning () {

	// basic variables
	double 	sum, sum_sq, mean,std_dev;
	int 	max, band, bands;

	int*	mark;

	// compute the number of bands
	bands = ceil(2*sqrt(N/2.));

	// allocate new memory
	mean_lane 	= malloc((steps+1)*sizeof(double));
	stddev_lane	= malloc((steps+1)*sizeof(double));
	max_lane	= malloc((steps+1)*bands*sizeof(double));
	y_lane		= malloc((steps+1)*bands*sizeof(double));

	// iterate over all time steps
	for (int i=0; i<=steps; i++) {

		// reset lane array
		for (int k=0; k<bands; k++) {
			max_lane[i*bands + k] = 0;
		}

		// reset both mean, standard deviation and sums
		mean 	= 0;
		std_dev = 0;
		sum 	= 0;
		sum_sq	= 0;
		max 	= 0;

		// iterate over all particles
		for (int j=0; j<N; j++) {
			sum 	+= laning[N*i+j];
			sum_sq	+= laning[N*i+j]*laning[N*i+j];

			// iterate the maximum lane length
			if (laning[N*i+j] > laning[N*i+max])
				max = j;

			// compute band which the particle lies in
			band = floor(2*config[2*N*i + 2*j+1]);

			// check whether a longer lane can be found in the specific band
			if (laning[N*i + j] > max_lane[i*bands + band]) {
				max_lane[i*bands+band] 	= laning[N*i + j];
				y_lane[i*bands+band]	= config[2*N*i + 2*j+1];
			}


		}		

		// calculate mean and standard deviation
		mean 	= sum/N;
		std_dev	= sqrt(sum_sq/N - mean*mean);

		// store in array
		mean_lane[i] 	= mean;
		stddev_lane[i] 	= std_dev;
	}
}


/*----------------------------------------------------------------------------------------------------------------------------*/
void stat_psi4 () {

	// basic variables
	double 	sum, sum_sq, mean,std_dev;

	// allocate new memory
	mean_psi4 	= malloc((steps+1)*sizeof(double));
	stddev_psi4	= malloc((steps+1)*sizeof(double));

	// iterate over all time steps
	for (int i=0; i<=steps; i++) {

		// reset both mean, standard deviation and sums
		mean 	= 0;
		std_dev = 0;
		sum 	= 0;
		sum_sq	= 0;

		// iterate over all particles
		for (int j=0; j<N; j++) {
			sum 	+= psi4[N*i+j];
			sum_sq	+= psi4[N*i+j]*psi4[N*i+j];
		}

		// calculate mean and standard deviation
		mean 	= sum/N;
		std_dev	= sqrt(sum_sq/N - mean*mean);

		// store in array
		mean_psi4[i] 	= mean;
		stddev_psi4[i] 	= std_dev;
	}
}


/*----------------------------------------------------------------------------------------------------------------------------*/
void stat_psi6 () {

	// basic variables
	double 	sum, sum_sq, mean,std_dev;

	// allocate new memory
	mean_psi6 	= malloc((steps+1)*sizeof(double));
	stddev_psi6	= malloc((steps+1)*sizeof(double));

	// iterate over all time steps
	for (int i=0; i<=steps; i++) {

		// reset both mean, standard deviation and sums
		mean 	= 0;
		std_dev = 0;
		sum 	= 0;
		sum_sq	= 0;

		// iterate over all particles
		for (int j=0; j<N; j++) {
			sum 	+= psi6[N*i+j];
			sum_sq	+= psi6[N*i+j]*psi6[N*i+j];
		}

		// calculate mean and standard deviation
		mean 	= sum/N;
		std_dev	= sqrt(sum_sq/N - mean*mean);

		// store in array
		mean_psi6[i] 	= mean;
		stddev_psi6[i] 	= std_dev;
	}
}


/*----------------------------------------------------------------------------------------------------------------------------*/
struct results* stat_analysis (struct parameters* param) {

	// variables from struct
	N 		= param->N;
	steps	= param->steps;

	config	= param->positions;
	disp	= param->displacement;

	laning	= param->laning;
	psi4	= param->psi4;
	psi6	= param->psi6;

	// call analysis functions
	stat_laning();
	stat_psi4();
	stat_psi6();

	// build struct to return to caller
	struct results* res = malloc(sizeof(struct results));

	res->mean_lane 		= mean_lane;
	res->mean_psi4 		= mean_psi4;
	res->mean_psi6 		= mean_psi6;

	res->stddev_lane	= stddev_lane;
	res->stddev_psi4	= stddev_psi4;
	res->stddev_psi6	= stddev_psi6;

	res->max_lane		= max_lane;
	res->y_lane			= y_lane;

	// return to caller
	return res;	
}