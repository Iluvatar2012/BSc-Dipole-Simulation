#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// define variables which will probably not be altered, except for hardcode
#define 	N 			1000
#define 	kT			1.0
#define 	tau_B 		1.0
#define 	D_A 		1.0
#define 	dt 			1e-6
#define 	iter 		1000000
#define 	threads 	8
#define 	writes 		1000

/*-------------------------------------------------------------------------------------------------------*/
/*	Template for config file:
*	Number of particles (N)              	: 1000
*	Energie (k_B*T)							: 1
*	Dipole interaction relation (Gamma_A)  	: 200
*	Particle dipole ratio (m)				: 0.1
*	Shear rate A (shear_A)               	: 0
*	Brownian Diffusion Time (tau_B)       	: 1
*	Brownian Diffusion (D_Brown_A)	       	: 1
*	Time difference (delta_t)            	: 1e-6
*	Destination file (outfile)           	: Results/great_file_name.hdf5
*	Iterations                           	: 100000
*	Number of Threads                    	: 8
*	Writeouts                            	: 1000
*/

/*-------------------------------------------------------------------------------------------------------*/
// The following function defines a builder for a series of config files with differing parameters

int builder (double m, double Gamma_A, double shear) {
	// construct the filename to which config data will be written
	char config[1024] = "Config_Files/";

	// construct the filename to which results will be written
	char outfile[1024] = "Results/";
	char buf[64];

	// append the number of particles
	strncat(outfile, "N_", 2);
	sprintf(buf, "%d", N);
	strncat(config, buf, 32);
	strncat(outfile, buf, 32);

	// append the binary relation m
	strncat(outfile, "__m_", 4);
	sprintf(buf, "%.3lf", m);
	strncat(config, buf, 64);
	strncat(outfile, buf, 64);

	// append the dipole relation Gamma
	strncat(outfile, "__GammaA_", 9);
	sprintf(buf, "%lf", Gamma_A);
	strncat(config, buf, 64);
	strncat(outfile, buf, 64);

	// append the shear rate
	strncat(outfile, "__shear_", 8);
	sprintf(buf, "%lf", shear);
	strncat(config, buf, 64);
	strncat(outfile, buf, 64);

	// append the proper ending
	strncat(outfile, ".hdf5", 5);


	// open the file, check whether successfull
	FILE *file = fopen(config, "w+");
	if(file == NULL) {
		fprintf(stderr, "The file \"%s\" could not be opened.\n", config);
		return EXIT_FAILURE;
	}

	// write all parameters to the file, check whether successful
	int check;
	if ((check = fprintf(file, "Number of particles (N)              	: %d\n",  N)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Energie (k_B*T)				      		: %lf\n", kT)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Dipole interaction relation (Gamma_A)	: %lf\n", Gamma_A)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Particle dipole ratio (m)				: %lf\n", m)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Shear rate A (shear_A)               	: %lf\n", shear)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Brownian Diffusion Time (tau_B)       	: %lf\n", tau_B)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Brownian Diffusion (D_Brown_A)	    	: %lf\n", D_A)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Time difference (delta_t)            	: %lf\n", dt)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Destination file (outfile)           	: %s\n",  outfile)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Iterations                           	: %d\n",  iter)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Number of Threads                     	: %d\n",  threads)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Writeouts                             	: %d",    writes)) < 1)
		return EXIT_FAILURE;

	// close the filestream and return to caller
	fclose(file);
	return EXIT_SUCCESS;
}

/*-------------------------------------------------------------------------------------------------------*/
int main (int argcount, char** argvector) {
}