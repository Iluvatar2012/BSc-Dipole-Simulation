#include <stdio.h>
#include <stdlib.h>



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