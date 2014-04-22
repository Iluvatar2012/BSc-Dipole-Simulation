#ifndef STRUCTS_H_
#define STRUCTS_H_

// struct for getting coniguration out of the hdf5 file
struct parameters {
	int 	N;
	int 	steps;
	double* positions;
};

// struct for passing parameters to analysis programs
struct variables {
	int 	N;

	double* positions;
	double* psi4;
	double* psi6;
	double* laning;
};

// struct for passing final data to file
struct analysis {
	int 	N;
	int 	steps;
	char*	file;

	double* psi4;
	double* psi6;
	double* laning;
};

#endif