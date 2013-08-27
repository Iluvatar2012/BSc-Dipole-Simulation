/*
 * simulation_gui.c
 *
 *  Created on: 	April 16, 2013
 *  Last Changed:	July 26, 2013
 *  Author: 		Aiko Bernehed
 */

// various includes of standard libraries, for manipulating strings
#include <stdlib.h>
#include <stdio.h>

#include <string.h>
#include <math.h>
#include <unistd.h>

#include <SDL/SDL.h>

#define PI 			3.14159265358979323846264338328


static int 		N;
static int 		steps;
static int* 	timestep;
static double* 	positions;

/*----------------------------------------------------------------------------------------------------------------------------*/
// TODO: write function description
int fileIO (char* file) {
	// Open file, check if successful
	FILE* infile = fopen(file, "r");
	if (infile == NULL) {
		fprintf(stderr, "File could not be opened, function exits...\nPath: %s", file);
		return EXIT_FAILURE;
	}

	// read N from file, exit if there is an error
	unsigned int temp;
	temp = fread(&N, sizeof(int), 1, infile);

	if(temp < 1) {
		fprintf(stderr, "N could not be read, function exits... Temp: %d\n", (int)(temp));
		return EXIT_FAILURE;
	}

	// read the amount of timesteps from the file, exit if there is an error
	temp = fread(&steps, sizeof(int), 1, infile);

	if( temp < 1) {
		fprintf(stderr, "steps could not be read from file, function exits... Temp: %d\n", (int)(temp));
		return EXIT_FAILURE;
	}

	// allocate memory for position and timestep array, check whether operation went well
	timestep = 	malloc(steps*sizeof(int));
	positions = malloc(steps*2*N*sizeof(double));

	if (timestep == NULL || positions == NULL) {
		fprintf(stderr, "Memory for timesteps and/or positions could not be allocated, function exits...\n");
		return EXIT_FAILURE;
	}

	// read timesteps
	for (int i=0; i<steps; i++) {
		fread(&(timestep[i]), sizeof(int), 1, infile);
		if ((temp = fread(&(positions[2*N*i]), sizeof(double), 2*N, infile)) < 2*N) {
			fprintf(stderr, "Error reading position data, function exits... Temp: %d\n", (int)(temp));
			return EXIT_FAILURE;
		}
	}

	return EXIT_SUCCESS;
}

/*----------------------------------------------------------------------------------------------------------------------------*/
// TODO: write function description
int graphicOutput () {

	// create everything we need to show the simulation
	SDL_Surface	*screen, *ball_black, *ball_red;
	SDL_Rect	dst;
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
	ball_black = SDL_LoadBMP("Black_Dot_3x3px.bmp");
	if (ball_black == NULL) {
		fprintf(stderr, "Could not load image of dot: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	// get the picture of the red dot
	ball_red = SDL_LoadBMP("Red_Dot_3x3px.bmp");
	if (ball_red == NULL) {
		fprintf(stderr, "Could not load image of dot: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	// get the height and width of our image, needed to blit it to screen
	dst.w = ball_black->w;
	dst.h = ball_black->h;

	dst.w = ball_red->w;
	dst.h = ball_red->h;

	// variables for checking various states (terminating, what frame to show and whether a key was pressed)
	int	done 		= 0;

	int	play 		= 0;
	int next 		= 0;
	int last 		= 0;

	int	keyPressed 	= 0;

	// the box's length
	double L 		= sqrt(N);

	// counter to set which timestep we are to draw
	int counter 	= 0;

	// variables to hold basic values for computing absolute pixel positions
	int	scrWidth 	= screen->w;
	int scrHeight 	= screen->h;

	double picWidth		= ball_black->w;
	double picHeight	= ball_black->h;

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
			counter %= steps;
			next = 0;
		}
		else if (last) {
			// decrease counter, check that we don't exceed the number of steps, reset manipulation variable
			counter--;
			counter += steps;
			counter %= steps;
			last = 0;
		}

		// Reset the screen
		SDL_FillRect(screen, NULL, SDL_MapRGB(screen->format, 255, 255, 255));

		// finally draw the current frame
		for (int i=0; i<N; i++) {
			// invert the y axis, otherwise (0,0) would be in the top left corner
			posY = -positions[counter*2*N+2*i+1]+L;

			// compute x and y position of each dot
			dst.x = round((positions[counter*2*N+2*i]/L)  *scrWidth -picWidth/2.);
			dst.y = round((posY/L)*scrHeight-picHeight/2.);

			// copy image to screen according to whether we need a red or black dot
			if (i%2 == 0)
				SDL_BlitSurface(ball_black, NULL, screen, &dst);
			else
				SDL_BlitSurface(ball_red, NULL, screen, &dst);
		}

		// flip screen and show other buffer
		SDL_Flip(screen);

		// add a delay of 40 ms, equals 25 fps
		SDL_Delay(40);
	}

	// clear memory of everything cluttering it
	SDL_FreeSurface(ball_black);
	SDL_FreeSurface(ball_red);

	// return to caller
	return EXIT_SUCCESS;
}

/*----------------------------------------------------------------------------------------------------------------------------*/
// TODO: write function description
int main (int argcount, char** argvektor) {

	// Read user defined binary file from stdin
	char infile[1024];
	getcwd(infile, sizeof(infile));
	strncat(infile, "/", 1);
	unsigned int length = sizeof(infile) - strlen(infile);

	if(argcount == 2) {
		strncat(infile, argvektor[1], length);
	} else {
		fprintf(stderr, "Please provide a binary file to be read! Function will now exit. \n");
		return EXIT_FAILURE;
	}

	// Read file, needed for parameters for graphic output
	int check = fileIO(infile);
	if (check != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

	// compute graphic output
	check = graphicOutput();
	if (check != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

	// return to caller
	return EXIT_SUCCESS;
}
