#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// define variables which will probably not be altered, except for hardcode
#define 	D_rat		1.7
#define 	dt 			1e-6
#define 	writestep	100000

// function declaration and static variables
int builder();
int read_build();

static int iter;


/*-------------------------------------------------------------------------------------------------------*/
/*	Template for config file:
*	Dipole interaction relation (Gamma_A)  	: 200
*	Particle dipole ratio (m)				: 0.1
*	Shear rate A (shear_A)               	: 200
* 	Particle Diffusion ratio (D_rat)		: 1.7
*	Time difference (delta_t)            	: 1e-6
*	Destination file (outfile)           	: Results/default.hdf5
*	Iterations                           	: 100000
*	Writeout step                          	: 1000
*/



/*-------------------------------------------------------------------------------------------------------*/
/*	This function reads from a provided builder file and builds the according files, template as follows:
*
*	iter 			1000000
*	m				0.1 		1.0 		0.1
*	Gamma 			10 			100 		10
*	shear 			0			1000 		20
*
*/
int read_build (char* infile) {

	int check = 1;

	// variables to be read
	double 	m_min, m_max, m_step;
	double 	G_min, G_max, G_step;
	double 	s_min, s_max, s_step;

	// open provided file, check if successful
	FILE *file = fopen(infile, "r");
	if(file == NULL) {
		fprintf(stderr, "The file \"%s\" could not be opened, program will now terminate. \n", infile);
		return EXIT_FAILURE;
	}

	// read parameters from provided file

	if (check != 0 && (check = fscanf(file, "iter 			%d\n", &iter)) < 1)
		check = 0;
	if (check != 0 && (check = fscanf(file, "m				%lf 		%lf 		%lf\n", &m_min, &m_max, &m_step)) < 1)
		check = 0;
	if (check != 0 && (check = fscanf(file, "Gamma 			%lf 		%lf	 		%lf\n", &G_min, &G_max, &G_step)) < 1)
		check = 0;
	if (check != 0 && (check = fscanf(file, "shear 			%lf			%lf 		%lf\n", &s_min, &s_max, &s_step)) < 1)
		check = 0;

	// check whether all reads were successful
	if (check < 1) {
		fprintf(stderr, "Variables could not be correctly read from file \"%s,\" program will now terminate. \n", infile);
		return EXIT_FAILURE;
	}

	// close the file
	fclose(file);

	// check whether the given variables actually make any sense
	if (m_min < 0 || m_min > m_max || \
		G_min < 0 || G_min > G_max || \
		s_min < 0 || s_min > s_max) {
		fprintf(stderr, "Variables either negative or minimum larger than maximum values, please check file and try again.\n");
		return EXIT_FAILURE;
	}

	// count the number of files we are writing
	int i = 1;

	// if everthing is in order, iterate over all provided variables and finally write the files
	for (double m = m_min; m <= m_max; m += m_step) {
		for (double G = G_min; G <= G_max; G += G_step) {
			for (double s = s_min; s <= s_max; s += s_step) {
				// call the file builder, increase file counter
				check = builder(m, G, s, i);
				i++;

				// check whether function worked properly
				if (check != EXIT_SUCCESS) {
					fprintf(stderr, "File could not be written, program will now terminate. \n");
					return EXIT_FAILURE;
				}
			}
		}
	}

	// return to caller
	return EXIT_SUCCESS;
}



/*-------------------------------------------------------------------------------------------------------*/
// The following function defines a builder for a series of config files with differing parameters

int builder (double m, double Gamma_A, double shear, int no) {
	// construct the filename to which config data will be written
	char config[1024] = "Config_Files/";

	// construct the filename to which results will be written
	char outfile[1024] = "Results/";
	char buf[64];

	sprintf(buf, "%d", no);
	strncat(config, buf, 32);

	// append the binary relation m
	strncat(outfile, "m_", 2);
	sprintf(buf, "%.2lf", m);
	strncat(outfile, buf, 64);

	// append the dipole relation Gamma
	strncat(outfile, "__GammaA_", 9);
	sprintf(buf, "%.1lf", Gamma_A);
	strncat(outfile, buf, 64);

	// append the shear rate
	strncat(outfile, "__shear_", 8);
	sprintf(buf, "%.0lf", shear);
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
	if ((check = fprintf(file, "Dipole interaction relation (Gamma_A)	: %lf\n", Gamma_A)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Particle dipole ratio (m)				: %lf\n", m)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Shear rate A (shear_A)               	: %lf\n", shear)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Particle Diffusion ratio (D_rat)		: %lf\n", D_rat)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Time difference (delta_t)            	: %lf\n", dt)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Destination file (outfile)           	: %s\n",  outfile)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Iterations                           	: %d\n",  iter)) < 1)
		return EXIT_FAILURE;
	if ((check = fprintf(file, "Writeout step                          	: %d\n",  writestep)) < 1)
		return EXIT_FAILURE;

	// close the filestream and return to caller
	fclose(file);
	return EXIT_SUCCESS;
}



/*-------------------------------------------------------------------------------------------------------*/
int main (int argcount, char** argvector) {

	int check;

	// check whether a builder file was provided at startup, this file is essential to altering parameters
	if (argcount != 2) {
		fprintf(stderr, "Please provide a builder file to be read from, program will now exit");
		return EXIT_FAILURE;
	}

	// read from provided file, build all files
	if ((check = read_build(argvector[1])) != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

	// return to caller
	return EXIT_FAILURE;
}
