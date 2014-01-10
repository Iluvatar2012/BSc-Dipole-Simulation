#include <math.h>
#include <SDL/SDL.h>
#include <stdlib.h>



/*----------------------------------------------------------------------------------------------------------------------------*/
int graphicOutput (char* file, double* config, double* psi_4, double* psi_6, int N) {

	// basic variables
	double dx, dy;
	double r_sq;


	// create everything we need to show the simulation
	SDL_Surface	*screen, *ball_even, *ball_uneven, *psi_4_ball, *psi_6_ball;
	SDL_Rect	dst_even, dst_uneven, dst_psi_4, dst_psi_6;

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

	// get the picture of the red dot
	ball_even = SDL_LoadBMP("Dots/Red_Dot_9x9px.bmp");
	if (ball_even == NULL) {
		fprintf(stderr, "Could not load image of dot: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	// get the picture of the green dot
	ball_uneven = SDL_LoadBMP("Dots/Green_Dot_5x5px.bmp");
	if (ball_uneven == NULL) {
		fprintf(stderr, "Could not load image of dot: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	// get the picture of the green dot
	psi_4_ball = SDL_LoadBMP("Dots/Orange_Dot_5x5px.bmp");
	if (psi_4_ball == NULL) {
		fprintf(stderr, "Could not load image of dot: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	// get the picture of the green dot
	psi_6_ball = SDL_LoadBMP("Dots/Blue_Dot_9x9px.bmp");
	if (psi_6_ball == NULL) {
		fprintf(stderr, "Could not load image of dot: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	// get the height and width of our image, needed to blit it to screen
	dst_even.w = ball_even->w;
	dst_even.h = ball_even->h;

	dst_uneven.w = ball_uneven->w;
	dst_uneven.h = ball_uneven->h;

	dst_psi_4.w = psi_4_ball->w;
	dst_psi_4.h = psi_4_ball->h;

	dst_psi_6.w = psi_6_ball->w;
	dst_psi_6.h = psi_6_ball->h;

	// the box's length
	double L 		= sqrt(N/2.);

	// variables to hold basic values for computing absolute pixel positions
	int	scrWidth 	= screen->w;
	int scrHeight 	= screen->h;

	double posY;

	// draw the default values for the current frame
	for (int i=0; i<N; i++) {
		
		// invert the y axis, otherwise (0,0) would be in the top left corner
		posY = -config[2*i+1]+L;

		// compute x and y position of each dot
		dst_even.x 	 = round((config[2*i]/L)  *scrWidth -ball_even->w/2.);
		dst_uneven.x = round((config[2*i]/L)  *scrWidth -ball_uneven->w/2.);

		dst_even.y   = round((posY/L)*scrHeight-ball_even->h/2.);
		dst_uneven.y = round((posY/L)*scrHeight-ball_uneven->h/2.);

		// copy image to screen according to whether we need a red or green dot
		if (i%2 == 0)
			SDL_BlitSurface(ball_even, NULL, screen, &dst_even);
		else
			SDL_BlitSurface(ball_uneven, NULL, screen, &dst_uneven);
	}

	// draw the objects having to do with Psi 4 and Psi 6 values
	for (int i=0; i<N; i++) {

		// invert the y axis, otherwise (0,0) would be in the top left corner
		posY = -config[2*i+1]+L;

		// compute x and y position of each dot
		dst_psi_4.x = round((config[2*i]/L)  *scrWidth -psi_4_ball->w/2.);
		dst_psi_6.x = round((config[2*i]/L)  *scrWidth -psi_6_ball->w/2.);

		dst_psi_4.y = round((posY/L)*scrHeight-psi_4_ball->h/2.);
		dst_psi_6.y = round((posY/L)*scrHeight-psi_6_ball->h/2.);

		// decide whether any dot needs a special coloring
		if (psi_4[i] > 0.9 || psi_6[i] > 0.9) {
			if (i%2 == 0)
				SDL_BlitSurface(psi_6_ball, NULL, screen, &dst_psi_6);
			else
				SDL_BlitSurface(psi_4_ball, NULL, screen, &dst_psi_4);
		}
	}

	// save the screen to file
	SDL_SaveBMP(screen, file);

	// clear memory of everything cluttering it
	SDL_FreeSurface(ball_even);
	SDL_FreeSurface(ball_uneven);

	// return to caller
	return EXIT_SUCCESS;
}
