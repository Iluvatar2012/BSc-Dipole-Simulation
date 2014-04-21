#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>

#include <string.h>
#include <pthread.h>

#include <unistd.h>

#include "functions.h"
#include "hdf5.h"
#include "structs.h"


static int 			steps, N;
static double* 		positions;

static double*		psi4, psi6, laning;

static int 			psi4_done, psi6_done, laning_done;

static pthread_t*	threads;

/*-------------------------------------------------------------------------------------------------------*/
void thread_psi4 () {
	for (int i=0; i<=(steps); i++) {
		compute_psi4(i);
	}

	psi4_done = 1;
}

/*-------------------------------------------------------------------------------------------------------*/
void thread_psi6() {
	for (int i=0; i<=(steps); i++) {
		compute_psi6(i);
	}

	psi6_done = 1;
}

/*-------------------------------------------------------------------------------------------------------*/
void thread_laning() {
	for (int i=0; i<=(steps); i++) {
		laning(i);
	}

	laning_done = 1;
}

/*-------------------------------------------------------------------------------------------------------*/
int main (int argcount, char** argvector) {

	int ret_thread;

	// check whether the right amount of arguments is given
	if (argcount != 2) {
		fprintf(stderr, "Please provide only a filename to read from, process will terminate.\n");
		return EXIT_FAILURE;
	}

	// read from the given file, check whether successful
	struct parameters *param = hdf5_read(argvector[1]);

	if (param == NULL) {
		return EXIT_FAILURE;
	}

	// copy data from struct
	N 			= param->N;
	steps		= param->steps;
	positions 	= param->positions;

	// allocate memory, check whether successful
	psi4 	= malloc((steps+1)*N*sizeof(double));
	psi6 	= malloc((steps+1)*N*sizeof(double));
	laning 	= malloc((steps+1)*N*sizeof(double));

	if (psi4 == NULL || psi6 == NULL || laning == NULL) {
		fprintf(stderr, "Memory for psi and laning values could not be allocated, process will terminate. \n");
		return EXIT_FAILURE;
	}

	// construct a new struct, initiate the helper file
	struct variables var = {.N = N, .positions = positions, .psi4 = .psi4, .psi6 = psi6, .laning = laning};
	init(&var);

	fprintf(stderr, "Starting computation. \n");

	// set waiting values
	psi4_done 	= 0;
	psi6_done 	= 0;
	laning_done = 0;

	// initiate threads, check if successful
	threads = malloc(3*sizeof(pthread_t));

	ret_thread 	= pthread_create(&(threads[0]), NULL, (void*)&thread_psi4, NULL);

	if (ret_thread != 0) {
		fprintf(stderr, "Thread %d could not be created. \n", 0);
		exit(EXIT_FAILURE);
	}
	
	ret_thread 	= pthread_create(&(threads[1]), NULL, (void*)&thread_psi6, NULL);

	if (ret_thread != 0) {
		fprintf(stderr, "Thread %d could not be created. \n", 1);
		exit(EXIT_FAILURE);
	}
	
	ret_thread 	= pthread_create(&(threads[2]), NULL, (void*)&thread_laning, NULL);

	if (ret_thread != 0) {
		fprintf(stderr, "Thread %d could not be created. \n", 2);
		exit(EXIT_FAILURE);
	}

	// wait until all threads have finished iterating
	while(psi4_done != 1 || psi6_done != 1 || laning_done != 1) {
		usleep(500*1000);
	}

	// TODO: save all data to file
}