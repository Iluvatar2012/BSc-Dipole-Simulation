#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <math.h>

#include "structs.h"
#include "hdf5.h"

// define the timestep, the writeout step of the simulation and the steplength to output data
#define	dt 			0.000001
#define	write_step 	1000

#define steplength 	1


// arrays holding the positions of the particles
static double* 	new_positions;
static double* 	old_positions;

// arrays for holding the speeds of the particles
static double* 	speeds;

// arrays for later sorting the speeds into isolines
static int* 	y_values;
static double* 	temp_speeds;
static double* 	temp_positions;

// variable for later sorting the speeds
static double 	y_step;

// variables holding basic values of the simulation
static int 		steps;
static int 		N;

static double 	L;

// output filename
static char		outfile[999];



/*----------------------------------------------------------------------------------------------------------------------------*/
double fast_abs(double d) {

	// create a union between a double and 64 bit integer (same memory space)
	union {
		double d;
		uint64_t u64;
	} u = {d};
	
	// shift the one bit all the way to the left and AND it with the union, this way the double value gets positiv
	u.u64 &= ~( (uint64_t) 1 << 63);

	// return the unions double value
	return u.d;
}



/*----------------------------------------------------------------------------------------------------------------------------*/
void calcSpeed (int step) {

	// basic variables
	double dx, dy;

	// read the current data from the provided hdf5 file
	hdf5_read (step, old_positions);
	hdf5_read (step+1, new_positions);

	// iterate over all particles
	for (int i=0; i<N; i++) {

		// calculate the difference in height and propagation
		dx = new_positions[2*i]   - old_positions[2*i];
		dy = new_positions[2*i+1] - old_positions[2*i+1];

		// calculate the absolute value of the distances
		dx = fast_abs(dx);
		dy = fast_abs(dy);

		// calculate whether the distance is larger than L/2 (periodic boundary conditions)
		dx = dx + (int)(2.*dx/L)*(L - 2*dx);
		dy = dy + (int)(2.*dy/L)*(L - 2*dy);

		// calculate the absolute speed
		speeds[i] = sqrt(dx*dx + dy*dy)/(dt*write_step);
	}
}



/*----------------------------------------------------------------------------------------------------------------------------*/
void sortSpeeds () {

	// basic variables
	int y_counter;
	int speed_counter;
	int min_index;

	double x;
	double min_rem_x;
	double min_abs_x;

	// initiate the first counter variable
	speed_counter = 0;

	// iterate over all isolines
	for (double y=0.0; y<L; y+=y_step) {

		// reset the y counter
		y_counter = 0;

		// iterate over all particles
		for (int i=0; i<N; i++) {

			// check that the y-value of the particle is within the desired strip
			if (old_positions[2*i+1] >= y && old_positions[2*i+1] < y + y_step) {

				// copy the index of the particle to the y_values array, increment the counter
				y_values[y_counter] = i;
				y_counter++;
			}
		}

		// reset the absolute minimum value of x
		min_abs_x = -1.0;

		// iterate over all particles in the y_values array
		for (int i=0; i<y_counter; i++) {

			// reset the minimum y variable
			min_rem_x = L+1.0;

			// iterate over all particles in the y_values array
			for (int j=0; j<y_counter; j++) {

				// store the current x-position for ease
				x = old_positions[2*y_values[j]];

				// find the minimum value of 
				if (x < min_rem_x && x > min_abs_x) {

					// store the index of the value and the found values index
					min_rem_x = x;
					min_index = y_values[j];
				} 
			}

			// set the found minimum as the absolute minimum in the remaining array
			min_abs_x = min_rem_x;

			// set the found minimum as the next particle
			temp_speeds[speed_counter] 			= speeds[min_index];
			temp_positions[2*speed_counter]		= old_positions[2*min_index];
			temp_positions[2*speed_counter+1]	= old_positions[2*min_index+1];

			// increment the speed counter
			speed_counter++;
		}
	}

	// we should now have iterated over all isolines and particles, now iterate over all particles and save the sorted system
	for (int i=0; i<N; i++) {

		// copy the speeds and particle positions back to the original arrays
		speeds[i] 				= temp_speeds[i];
		old_positions[2*i]		= temp_positions[2*i];
		old_positions[2*i+1]	= temp_positions[2*i+1];
	}
}



/*----------------------------------------------------------------------------------------------------------------------------*/
void write_to_file (int number) {

	// basic variables
	char filename[1024];
	char buf[24];

	double y = y_step;

	// copy the outgoing filename into local variable
	strncpy(filename, outfile, sizeof(outfile));

	// append an underscore and the iteration number
	strncat(filename, "_", 1);
	sprintf(buf, "%d", number);
	strncat(filename, buf, 24);

	// open the file, check whether successfull
	FILE *file = fopen(filename, "w+");
	if(file == NULL) {
		fprintf(stderr, "The file \"%s\" could not be opened, process will terminate. \n", filename);
		exit(EXIT_FAILURE);
	}

	// iterate over all particles
	for (int i=0; i<N; i++) {

		// check if we have hit a new isoline
		if (old_positions[2*i+1] > y) {
			
			// write an empty line and set the checking variable to the next isoline height
			//fprintf(file, "\n");
			y+=y_step;
		}

		// write the positions and absolute speeds of all particles into the file
		fprintf(file, "%lf\t%lf\t%lf\n", old_positions[2*i], old_positions[2*i+1], speeds[i]);
	}

	// close the filestream
	fclose(file);
}


/*----------------------------------------------------------------------------------------------------------------------------*/
void readFile(int argcount, char** argvector) {

	// basic variables
	char file[1024];

	// check whether the right amount of arguments is given
	if (argcount != 3) {
		fprintf(stderr, "Please provide a filename to read from and a file to write to, process will terminate.\n");
		exit(EXIT_FAILURE);
	}

	// copy the filenames
	strncpy(file, argvector[1], sizeof(file));
	strncpy(outfile, argvector[2], sizeof(outfile));

	// read from the given file, check whether successful
	struct parameters *param = hdf5_init(file);

	// copy the arrays and values from the struct
	steps		= param->steps;
	N 			= param->N;

	// initiate the arrays
	new_positions 	= malloc(2*N*sizeof(double));
	old_positions 	= malloc(2*N*sizeof(double));
	speeds			= malloc(N*sizeof(double));

	y_values		= malloc(N*sizeof(int));
	temp_speeds		= malloc(N*sizeof(double));
	temp_positions 	= malloc(2*N*sizeof(double));

	// calculate the size of the system and the stepsize of x
	L 		= sqrt(N/2.);
	y_step 	= L/25.;
}



/*----------------------------------------------------------------------------------------------------------------------------*/
int main(int argcount, char** argvector) {

	// basic counter
	int counter		= 0;

	// read the provided file
	readFile(argcount, argvector);

	// iterate over all steps, adjust the steplength accordingly
	for (int i=0; i<steps; i += steplength) {

		// calculate the current speeds of all particles
		calcSpeed(i);

		// sort the arrays into isolines for pm3d to scan
		//sortSpeeds();

		// write the speeds to their respective files
		write_to_file(counter);

		// increment the counter
		counter++;
	}

}