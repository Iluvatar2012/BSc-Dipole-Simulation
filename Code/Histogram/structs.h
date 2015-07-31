#ifndef STRUCTS_H_
#define STRUCTS_H_

// struct for getting coniguration out of the hdf5 file
struct parameters {
	int 	N;
	int 	steps;
	int 	last_step;

	double 	X;
	double 	L_y;

	double* positions;

};

#endif