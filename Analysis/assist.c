#include <math.h>

#include <stdio.h>
#include <stdlib.h>



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
			real	+= cos(n*theta);
			img 	+= sin(n*theta);
	}

	// normalize real and imaginary part and calculate absolute value
	real /= n;
	img  /= n;

	return sqrt(real*real + img*img);
}