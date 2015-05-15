#ifndef STRUCTS_H_
#define STRUCTS_H_

struct parameters {
	double* 	positions;
	double* 	displacement;

	double* 	laning;
	double* 	psi4;
	double* 	psi6;

	int 		N;
	int 		steps;
};

struct results {
	double* 	mean_lane;
	double* 	mean_psi4;
	double* 	mean_psi6;

	double*		stddev_lane;
	double* 	stddev_psi4;
	double* 	stddev_psi6;

	double* 	max_lane;
	double* 	y_lane;
};

#endif