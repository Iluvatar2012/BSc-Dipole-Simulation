#include <math.h>

#include <stdio.h>
#include <stdlib.h>

#include "assist.h"



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
double psi_n (int n, int i, double* config, int* next) {

	// basic variables
	double dx, dy;
	double theta;

	// initiate real and imaginary part
	double real = 0;
	double img  = 0;

	// iterate over the saved particles
	for (int k=0; k<n; k++) {
			dx = config[2*next[k]] 	 - config[2*i];
			dy = config[2*next[k]+1] - config[2*i+1];

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



/*----------------------------------------------------------------------------------------------------------------------------*/
double bond_length_n (int n, int i, double* config, int* next) {
	
	// basic variables
	double dx, dy;
	double l_med = 0;
	double l;
	double b = 0;

	// calculate the medium distance
	for (int k=0; i<n; k++) {
		dx = config[2*next[k]] 	 - config[2*i];
		dy = config[2*next[k]+1] - config[2*i+1];

		l_med += sqrt(dx*dx + dy*dy);
	}

	l_med /= n;

	// calculate the bond_length_deviation
	for (int k=0; k<n; k++) {
		dx = config[2*next[k]] 	 - config[2*i];
		dy = config[2*next[k]+1] - config[2*i+1];

		l_med += sqrt(dx*dx + dy*dy);

		b += abs(l-l_med)/l_med;
	}

	b /= n;

	fprintf(stderr, "i: %d bld: %lf, x: %lf, y: %lf\n", i, b, config[2*i], config[2*i+1]);

	// return to caller
	return b;
}