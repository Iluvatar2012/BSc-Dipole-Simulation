#ifndef STRUCTS_H_
#define STRUCTS_H_

// struct for getting coniguration out of the hdf5 file
struct parameters {
	int 	N;
	int 	steps;

	double* positions;
	double* displacement;
	double* psi4;
	double* psi6;
	double* laning;
};

#endif