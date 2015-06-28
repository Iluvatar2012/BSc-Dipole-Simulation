#ifndef STRUCTS_H_
#define STRUCTS_H_

// struct for getting coniguration out of the hdf5 file
struct parameters {
	int 	N;
	int 	steps;

	double 	L_x;
	double 	L_y;
	double 	X;

	double* positions;
};

#endif