/*
 * simulation_gui.c
 *
 *  Created on: 	April 16, 2013
 *  Last Changed:	July 26, 2013
 *  Author: 		Aiko Bernehed
 */

#include <stdlib.h>
#include <stdio.h>

#include <string.h>
#include <math.h>
#include <unistd.h>

#include <SDL/SDL.h>
#include <hdf5.h>

#include "assist_gui.h"


// #define PI 			3.14159265358979323846264338328

static int 		N;
static int 		steps;
static double* 	positions;
static double*	psi4;
static double*	psi6;

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
// TODO: write function description
int graphicOutput () {

	// create everything we need to show the simulation
	SDL_Surface	*screen, *ball_A, *ball_B, *psi_4_A, *psi_4_B, *psi_6_A, *psi_6_B;
	SDL_Rect	dst_even, dst_uneven;
	SDL_Event	event;

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
	ball_A = SDL_LoadBMP("Dots/Grey_Dot_9x9px.bmp");
	if (ball_A == NULL) {
		fprintf(stderr, "Could not load image of dot: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	// get the picture of the red dot
	ball_B = SDL_LoadBMP("Dots/Grey_Dot_5x5px.bmp");
	if (ball_B == NULL) {
		fprintf(stderr, "Could not load image of dot: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	psi_4_A = SDL_LoadBMP("Dots/Red_Dot_9x9px.bmp");
	if (psi_4_A == NULL) {
		fprintf(stderr, "Could not load image of dot: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	psi_4_B = SDL_LoadBMP("Dots/Red_Dot_5x5px.bmp");
	if (psi_4_B == NULL) {
		fprintf(stderr, "Could not load image of dot: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	psi_6_A = SDL_LoadBMP("Dots/Blue_Dot_9x9px.bmp");
	if (psi_6_A == NULL) {
		fprintf(stderr, "Could not load image of dot: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	psi_6_B = SDL_LoadBMP("Dots/Blue_Dot_5x5px.bmp");
	if (psi_6_B == NULL) {
		fprintf(stderr, "Could not load image of dot: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	// get the height and width of our image, needed to blit it to screen
	dst_even.w = ball_A->w;
	dst_even.h = ball_A->h;

	dst_uneven.w = ball_B->w;
	dst_uneven.h = ball_B->h;

	// variables for checking various states (terminating, what frame to show and whether a key was pressed)
	int	done 		= 0;

	int	play 		= 0;
	int next 		= 0;
	int last 		= 0;

	int	keyPressed 	= 0;

	// the box's length
	double L 		= sqrt(N/2.);

	// counter to set which timestep we are to draw
	int counter 	= 0;

	// variables to hold basic values for computing absolute pixel positions
	int	scrWidth 	= screen->w;
	int scrHeight 	= screen->h;

	double posY;

	// this will later hold our input from keyboard
	Uint8	*keys;

	// the entire video loop, will run until user terminates the program
	while (!done) {
		// Event loop, this will check whether any event is registered
		while (SDL_PollEvent(&event)) {
			// check what event has happened
			switch (event.type) {
			// terminate the program
			case SDL_QUIT:
				done = 1;
				break;

			// user input, get the current keyboard state
			case SDL_KEYDOWN:
				keys = SDL_GetKeyState(NULL);
				keyPressed = 1;
				break;
			}
		}

		// if a key was pressed, check user input
		if (keyPressed) {
			// reset variable
			keyPressed = 0;

			// invert varibale to play all frames
			if (keys[SDLK_SPACE]) {
				if (play)
					play = 0;
				else
					play = 1;
				// reset key variable
				keys[SDLK_SPACE] = 0;
			}

			// check whether we are supposed to show the next frame
			if (keys[SDLK_RIGHT]) {
				next = 1;
				// reset key variable
				keys [SDLK_RIGHT] = 0;
			}

			// check whether we are supposed to show the last frame
			if (keys[SDLK_LEFT] && !play) {
				last = 1;
				// reset key variable
				keys[SDLK_LEFT] = 0;
			}
		}

		// set counter to next or last frame
		if (next || play) {
			// increase counter, check that we don't exceed the number of steps, reset manipulation variable
			counter++;
			counter %= (steps + 1);

			next = 0;
		}
		else if (last) {
			// decrease counter, check that we don't exceed the number of steps, reset manipulation variable
			counter--;
			counter += (steps + 1);
			counter %= (steps + 1);
			last = 0;
		}

		// Reset the screen
		SDL_FillRect(screen, NULL, SDL_MapRGB(screen->format, 255, 255, 255));

		// finally draw the current frame
		for (int i=0; i<N; i++) {
			// invert the y axis, otherwise (0,0) would be in the top left corner
			posY = -positions[counter*2*N+2*i+1]+L;

			// set the alpha values of all overlays
			SDL_SetAlpha(psi_4_A, SDL_SRCALPHA, (int)(255*psi4[counter*N+i]));
			SDL_SetAlpha(psi_4_B, SDL_SRCALPHA, (int)(255*psi4[counter*N+i]));
			SDL_SetAlpha(psi_6_A, SDL_SRCALPHA, (int)(255*psi6[counter*N+i]));
			SDL_SetAlpha(psi_6_B, SDL_SRCALPHA, (int)(255*psi6[counter*N+i]));

			// copy image to screen according to whether we need a red or black dot
			if (i%2 == 0) {
				// compute x and y position of each dot
				dst_even.x = round((positions[counter*2*N+2*i]/L) *scrWidth - ball_A->w/2.);
				dst_even.y = round((posY/L)*scrHeight - ball_A->h/2.);
				SDL_BlitSurface(ball_A, NULL, screen, &dst_even);
				SDL_BlitSurface(psi_4_A, NULL, screen, &dst_even);
				SDL_BlitSurface(psi_6_A, NULL, screen, &dst_even);
			}
			else {
				// compute x and y position of each dot
				dst_uneven.x = round((positions[counter*2*N+2*i]/L) *scrWidth - ball_B->w/2.);
				dst_uneven.y = round((posY/L)*scrHeight - ball_B->h/2.);
				SDL_BlitSurface(ball_B, NULL, screen, &dst_uneven);
				SDL_BlitSurface(psi_4_B, NULL, screen, &dst_uneven);
				SDL_BlitSurface(psi_6_B, NULL, screen, &dst_uneven);
			}
		}

		// flip screen and show other buffer
		SDL_Flip(screen);

		// add a delay of 40 ms, equals 25 fps
		SDL_Delay(40);
	}

	// clear memory of everything cluttering it
	SDL_FreeSurface(ball_A);
	SDL_FreeSurface(ball_B);

	// return to caller
	return EXIT_SUCCESS;
}

/*----------------------------------------------------------------------------------------------------------------------------*/
int computePsi() {

	// allocate memory for psi values, check whether successfull
	psi4 = malloc((steps+1)*N*sizeof(double));
	psi6 = malloc((steps+1)*N*sizeof(double));

	if (psi4 == NULL || psi6 == NULL) {
		fprintf(stderr, "Memory for psi values could not be allocated, process will terminate. \n");
		return EXIT_FAILURE;
	}

	// initiate psi calculation
	init_psi(N, positions, psi4, psi6);

	// iterate over all steps
	for(int i=0; i<=steps; i++) {
		compute_psi4(i);
		compute_psi6(i);
	}

	return EXIT_SUCCESS;
}



/*----------------------------------------------------------------------------------------------------------------------------*/
// TODO: write function description
int main (int argcount, char** argvektor) {

	// Read user defined binary file from stdin
	char infile[1024];

	if(argcount == 2) {
		strncpy(infile, argvektor[1], 1024);
	} else {
		fprintf(stderr, "Please provide a binary file to be read! Function will now exit. \n");
		return EXIT_FAILURE;
	}

	// Read file, needed for parameters for graphic output
	int check = hdf5_read(infile);
	if (check != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

	// compute psi values
	check = computePsi();
	if (check != EXIT_SUCCESS)
		return EXIT_FAILURE;

	// compute graphic output
	check = graphicOutput();
	if (check != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

	// return to caller
	return EXIT_SUCCESS;
}
