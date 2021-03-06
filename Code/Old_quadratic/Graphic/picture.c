#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <SDL/SDL.h>

#include "structs.h"



// basic variables
static double 	L;
static int 		N;

static double* 	positions;
static double* 	psi4;
static double* 	psi6;

static int 		scrWidth;
static int 		scrHeight;

// SDL variables 
static SDL_Surface 	*screen;
static SDL_Surface 	*ball_A;
static SDL_Surface 	*ball_B;
static SDL_Surface 	*psi_4_A;
static SDL_Surface 	*psi_4_B;
static SDL_Surface 	*psi_6_A;
static SDL_Surface 	*psi_6_B;

static SDL_Rect 	dst_even;
static SDL_Rect 	dst_uneven;

static SDL_Event	event;



/*----------------------------------------------------------------------------------------------------------------------------*/
int initiate (struct parameters* param) {

	// get all necessary information from the struct
	N 			= param->N;

	positions 	= param->positions;
	psi4 		= param->psi4;
	psi6 		= param->psi6;

	// compute boxlength
	L = sqrt(N/2.);

	// initiate the SDL video mode, terminate if there was an error
	if ((SDL_Init(SDL_INIT_VIDEO)) == -1) {
		fprintf(stderr, "Could not initiate SDL video mode: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	// tell the exit function to also terminate all SDL functions
	atexit(SDL_Quit);

	// set up the screen, this will be our frame, terminate if there was an error
	screen = SDL_SetVideoMode(1200, 1200, 16, SDL_HWSURFACE | SDL_DOUBLEBUF);
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

	scrWidth 	= screen->w;
	scrHeight 	= screen->h;

	// return to caller
	return EXIT_SUCCESS;
}



/*----------------------------------------------------------------------------------------------------------------------------*/
void destroy () {

	// clear memory of everything cluttering it
	SDL_FreeSurface(screen);
	SDL_FreeSurface(ball_A);
	SDL_FreeSurface(ball_B);
	SDL_FreeSurface(psi_4_A);
	SDL_FreeSurface(psi_4_B);
	SDL_FreeSurface(psi_6_A);
	SDL_FreeSurface(psi_6_B);
}



/*----------------------------------------------------------------------------------------------------------------------------*/
void draw_picture(int step, char* file) {

	// basic variables
	double 	posY;

	// fill with a blank screen
	SDL_FillRect(screen, NULL, SDL_MapRGB(screen->format, 255, 255, 255));


	// finally draw the current frame
	for (int i=0; i<N; i++) {
		// invert the y axis, otherwise (0,0) would be in the top left corner
		posY = -positions[2*i+1]+L;

		// set the alpha values of all overlays
		SDL_SetAlpha(psi_4_A, SDL_SRCALPHA, (int)(255*psi4[i]));
		SDL_SetAlpha(psi_4_B, SDL_SRCALPHA, (int)(255*psi4[i]));
		SDL_SetAlpha(psi_6_A, SDL_SRCALPHA, (int)(255*psi6[i]));
		SDL_SetAlpha(psi_6_B, SDL_SRCALPHA, (int)(255*psi6[i]));

		// copy image to screen according to whether we need a red or black dot
		if (i%2 == 0) {
			// compute x and y position of each dot
			dst_even.x = round((positions[2*i]/L) *scrWidth - ball_A->w/2.);
			dst_even.y = round((posY/L)*scrHeight - ball_A->h/2.);
			SDL_BlitSurface(ball_A, NULL, screen, &dst_even);
			SDL_BlitSurface(psi_4_A, NULL, screen, &dst_even);
			SDL_BlitSurface(psi_6_A, NULL, screen, &dst_even);
		}
		else {
			// compute x and y position of each dot
			dst_uneven.x = round((positions[2*i]/L) *scrWidth - ball_B->w/2.);
			dst_uneven.y = round((posY/L)*scrHeight - ball_B->h/2.);
			SDL_BlitSurface(ball_B, NULL, screen, &dst_uneven);
			SDL_BlitSurface(psi_4_B, NULL, screen, &dst_uneven);
			SDL_BlitSurface(psi_6_B, NULL, screen, &dst_uneven);
		}
	}

	// Save the screen to file
	SDL_SaveBMP(screen, file);
}