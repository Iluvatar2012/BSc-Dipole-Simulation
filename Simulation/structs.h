/*
 * structs.h
 *
 *  Created on: Aug 21, 2013
 *      Author: Aiko Bernehed
 */

#ifndef STRUCTS_H_
#define STRUCTS_H_

/*-------------------------------------------------------------------------------------------------------*/
// define a struct for retrieving parameters from a given file
struct parameters {
	// file to be written to
	char 	outfile[2048];

	// system properties
	int 	N;
	double	kT;
	double 	Gamma_A;
	double 	m;
	double 	shear_A;
	double	tau_B;
	double 	D_Brown_A;

	// simulation and computer system properties
	double 	timestep;
	int		max_timesteps;
	int		thread_number;
	int		no_writeouts;
};
/*-------------------------------------------------------------------------------------------------------*/

#endif /* STRUCTS_H_ */