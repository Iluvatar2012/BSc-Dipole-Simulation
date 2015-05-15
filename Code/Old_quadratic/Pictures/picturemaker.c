#include <stdlib.h>
#include <stdio.h>

#include <string.h>
#include <math.h>
#include <unistd.h>

#include <SDL/SDL.h>
#include <hdf5.h>


#define PI 			3.14159265358979323846264338328

static int 		N;
static int 		steps;
static double* 	positions;


/*----------------------------------------------------------------------------------------------------------------------------*/
int hdf5_read (char* file) {

	// basic variables
	double* temp;

	// identifiers for files and a status variable
	hid_t	file_id, dataset_id, tempset_id, attr_write_id, attr_N_id;
	hid_t	dataspace_id;
	herr_t	status;

	// open the file, check whether operation was successful
	file_id			= H5Fopen(file, H5F_ACC_RDONLY, H5P_DEFAULT);

	if (file_id < 0) {
		fprintf(stderr, "File could not be opened, program will now terminate!\nPath: %s\n", file);
		return EXIT_FAILURE;
	}

	// open the dataset and attributes
	dataset_id		= H5Dopen2(file_id, "/positions", H5P_DEFAULT);
	attr_N_id		= H5Aopen(dataset_id, "N", H5P_DEFAULT);
	attr_write_id	= H5Aopen(dataset_id, "Writeouts", H5P_DEFAULT);

	// read the attributes
	status			= H5Aread(attr_N_id, H5T_NATIVE_INT, &N);
	status			= H5Aread(attr_write_id, H5T_NATIVE_INT, &steps);

	// allocate memory for the temporary data and positions
	temp			= malloc(2*N*sizeof(double));
	positions 		= malloc((steps+1)*2*N*sizeof(double));

	// offset and dimension of a hyperslab
	hsize_t	offset[2]	= {0, 0};
	hsize_t	slabdim[2]	= {1, 2*N};

	// get the data space of the current dataset and initialize a buffer for reading data
	dataspace_id	= H5Dget_space(dataset_id);
	tempset_id 		= H5Screate_simple(2, slabdim, NULL);

	// iterate over all steps and copy them all into the position array
	for (int i=0; i<=steps; i++) {

		// adjust offset to the next slab
		offset[0] = i;

		// select the hyperslab and read into position array
		status 		= H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, slabdim, NULL);
		status 		= H5Dread(dataset_id, H5T_NATIVE_DOUBLE, tempset_id, dataspace_id, H5P_DEFAULT, temp);

		// copy the data from the buffer to the position array
		for (int j=0; j<2*N; j+=2) {
			positions[2*N*i + j] 	= temp[j];
			positions[2*N*i + j+1]	= temp[j+1];
		}
	}

	// close all parts of the simulation
	status = H5Aclose(attr_N_id);
	status = H5Aclose(attr_write_id);
	status = H5Sclose(tempset_id);
	status = H5Sclose(dataspace_id);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);

	// return to caller
	return EXIT_SUCCESS;
}


/*----------------------------------------------------------------------------------------------------------------------------*/
int graphicOutput (char* file, char* read) {

	// create everything we need to show the simulation
	SDL_Surface	*screen, *ball_even, *ball_uneven;
	SDL_Rect	dst_even, dst_uneven;

	// initiate the SDL video mode, terminate if there was an error
	if ((SDL_Init(SDL_INIT_VIDEO)) == -1) {
		fprintf(stderr, "Could not initiate SDL video mode: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	// tell the exit function to also terminate all SDL functions
	atexit(SDL_Quit);

	// set up the screen, this will be our frame, terminate if there was an error
	screen = SDL_SetVideoMode(600, 600, 16, SDL_HWSURFACE | SDL_DOUBLEBUF);
	SDL_FillRect(screen, NULL, SDL_MapRGB(screen->format, 255, 255, 255));
	if (screen == NULL) {
		fprintf(stderr, "Could not set video mode: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	// get the picture of the black dot
	ball_even = SDL_LoadBMP("Dots/Red_Dot_9x9px.bmp");
	if (ball_even == NULL) {
		fprintf(stderr, "Could not load image of dot: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	// get the picture of the red dot
	ball_uneven = SDL_LoadBMP("Dots/Green_Dot_5x5px.bmp");
	if (ball_uneven == NULL) {
		fprintf(stderr, "Could not load image of dot: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	// get the height and width of our image, needed to blit it to screen
	dst_even.w = ball_even->w;
	dst_even.h = ball_even->h;

	dst_uneven.w = ball_uneven->w;
	dst_uneven.h = ball_uneven->h;

	// variable for storing the step that is shall be drawn
	int step;

	// the box's length
	double L 		= sqrt(N/2.);

	// variables to hold basic values for computing absolute pixel positions
	int	scrWidth 	= screen->w;
	int scrHeight 	= screen->h;

	double posY;

	char buf[16];

	// check which step is supposed to be drawn
	if (strcmp(read, "first") == 0) 
		step = 0;
	else if (strcmp(read, "last") == 0)
		step = steps;
	else
		step = atoi(read);

	// add the step drawn and file extension to the filename
	strncat(file, "__write_", 8);
	sprintf(buf, "%d", step);
	strncat(file, buf, 16);
	strncat(file, ".bmp", 4);

	// finally draw the current frame
	for (int i=0; i<N; i++) {
		// invert the y axis, otherwise (0,0) would be in the top left corner
		posY = -positions[step*2*N+2*i+1]+L;

		// copy image to screen according to whether we need a red or green dot
		if (i%2 == 0) {
			// compute x and y position of each even dot
			dst_even.x = round((positions[step*2*N+2*i]/L)  *scrWidth - ball_even->w/2.);
			dst_even.y = round((posY/L)*scrHeight - ball_even->h/2.);
			SDL_BlitSurface(ball_even, NULL, screen, &dst_even);
		}
		else {
			// compute x and y position of each uneven dot
			dst_uneven.x = round((positions[step*2*N+2*i]/L)  *scrWidth - ball_uneven->w/2.);
			dst_uneven.y = round((posY/L)*scrHeight - ball_uneven->h/2.);
			SDL_BlitSurface(ball_uneven, NULL, screen, &dst_uneven);
		}
	}

	// Save the screen to file
	SDL_SaveBMP(screen, file);

	// clear memory of everything cluttering it
	SDL_FreeSurface(ball_even);
	SDL_FreeSurface(ball_uneven);

	// return to caller
	return EXIT_SUCCESS;
}


/*----------------------------------------------------------------------------------------------------------------------------*/
int main (int argcount, char** argvektor) {

	char 	file[1024];
	char	temp[1024];
	char 	step[8];
	char* 	ptr;

	// read user defined hdf5 file from stdin
	if(argcount == 3) {
		strncpy(file, argvektor[1], 1023);
		strncpy(step, argvektor[2], 8);
	} else {
		fprintf(stderr, "Please provide an HDF5 file to be read and the timestep ('first', 'last' or number) to be read. \n");
		return EXIT_FAILURE;
	}

	// Read file, needed for parameters for graphic output
	int check = hdf5_read(file);
	if (check != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

	// compute the name of the output file, separate the last part of the string
	ptr = strrchr(file,  '/');

	if (ptr != NULL)
		strncpy(temp, ptr+1, 1023);
	else
		strncpy(temp, file, 1023);

	// remove the last part containing the file extension
	ptr 	= strrchr(temp, '.');
	*ptr 	= '\0';

	// add the appropriate directory to outfile path
	strncpy(file, "Pictures/\0", 10);
	strncat(file, temp, 1013);
	
	// compute graphic output
	check = graphicOutput(file, step);
	if (check != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

	// return to caller
	return EXIT_SUCCESS;
}