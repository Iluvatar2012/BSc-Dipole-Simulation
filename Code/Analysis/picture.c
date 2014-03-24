#include <math.h>
#include <SDL/SDL.h>
#include <stdlib.h>



/*----------------------------------------------------------------------------------------------------------------------------*/
int graphicOutput (char* file, double* config, double* psi_4, double* psi_6, int N) {

	// basic variables
	double dx, dy;
	double r_sq;


	// create everything we need to show the simulation
	SDL_Surface	*screen, *ball_A, *ball_B, *psi_4_A, *psi_4_B, *psi_6_A, *psi_6_B;
	SDL_Rect	dst_A, dst_B;

	// initiate the SDL video mode, terminate if there was an error
	if ((SDL_Init(SDL_INIT_VIDEO)) == -1) {
		fprintf(stderr, "Could not initiate SDL video mode: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	// tell the exit function to also terminate all SDL functions
	atexit(SDL_Quit);

	// set up the screen, this will be our frame, terminate if there was an error
	screen = SDL_SetVideoMode(600, 600, 32, SDL_HWSURFACE | SDL_DOUBLEBUF);
	SDL_FillRect(screen, NULL, SDL_MapRGB(screen->format, 255, 255, 255));
	if (screen == NULL) {
		fprintf(stderr, "Could not set video mode: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

	// get the picture of all particles
	ball_A = SDL_LoadBMP("Dots/Grey_Dot_9x9px.bmp");
	if (ball_A == NULL) {
		fprintf(stderr, "Could not load image of dot: %s\n", SDL_GetError());
		return EXIT_FAILURE;
	}

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
	dst_A.w = ball_A->w;
	dst_A.h = ball_A->h;

	dst_B.w = ball_B->w;
	dst_B.h = ball_B->h;

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
		dst_A.x = round((config[2*i]/L)  *scrWidth -ball_A->w/2.);
		dst_B.x = round((config[2*i]/L)  *scrWidth -ball_B->w/2.);

		dst_A.y = round((posY/L)*scrHeight-ball_A->h/2.);
		dst_B.y = round((posY/L)*scrHeight-ball_B->h/2.);

		// copy image to screen according to whether we need a red or green dot
		if (i%2 == 0)
			SDL_BlitSurface(ball_A, NULL, screen, &dst_A);
		else
			SDL_BlitSurface(ball_B, NULL, screen, &dst_B);

		// set the alpha values of all overlays
		SDL_SetAlpha(psi_4_A, SDL_SRCALPHA, (int)(255*psi_4[i]));
		SDL_SetAlpha(psi_4_B, SDL_SRCALPHA, (int)(255*psi_4[i]));
		SDL_SetAlpha(psi_6_A, SDL_SRCALPHA, (int)(255*psi_6[i]));
		SDL_SetAlpha(psi_6_B, SDL_SRCALPHA, (int)(255*psi_6[i]));

		// draw the psi_4 value for the current dot
		if (i%2 == 0)
			SDL_BlitSurface(psi_4_A, NULL, screen, &dst_A);
		else
			SDL_BlitSurface(psi_4_B, NULL, screen, &dst_B);		

		// draw the psi_6 value for the current dot
		if (i%2 == 0)
			SDL_BlitSurface(psi_6_A, NULL, screen, &dst_A);
		else
			SDL_BlitSurface(psi_6_B, NULL, screen, &dst_B);
	}

	// save the screen to file
	SDL_SaveBMP(screen, file);

	// clear memory of everything cluttering it
	SDL_FreeSurface(ball_A);
	SDL_FreeSurface(ball_B);
	SDL_FreeSurface(psi_4_A);
	SDL_FreeSurface(psi_4_B);
	SDL_FreeSurface(psi_6_A);
	SDL_FreeSurface(psi_6_B);
	SDL_FreeSurface(screen);

	// return to caller
	return EXIT_SUCCESS;
}
