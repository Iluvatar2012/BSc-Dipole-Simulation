#include <stdlib.h>
#include <stdio.h>

#include <hdf5.h>

#include "structs.h"


// basic variables
static double* 	positions;

static int 		N;
static int 		steps;

static double	L_x;
static double	L_y;
static double 	X;

// identifiers for files, dataspaces and -sets
static hid_t	file_id; 
static hid_t	position_set;
static hid_t	tempset_id;
static hid_t	position_space;

// data dimensions
static hsize_t	offset[2];
static hsize_t 	slabdim[2];


/*----------------------------------------------------------------------------------------------------------------------------*/
void hdf5_read (int curr_step) {
	// checking variable
	herr_t	status;

	// adjust the current offset
	offset[0] 		= curr_step;

	// select the hyperslab and read into position array
	status 	= H5Sselect_hyperslab(position_space, H5S_SELECT_SET, offset, NULL, slabdim, NULL);
	status 	= H5Dread(position_set, H5T_NATIVE_DOUBLE, tempset_id, position_space, H5P_DEFAULT, positions);
}


/*----------------------------------------------------------------------------------------------------------------------------*/
struct parameters *hdf5_init (char* file) {


	hid_t	attr_write_id, attr_N_id, attr_L_x_id, attr_L_y_id, attr_X_id;
	herr_t	status;

	// open file, check if successful
	file_id = H5Fopen(file, H5F_ACC_RDONLY, H5P_DEFAULT);

	if (file_id < 0) {
		fprintf(stderr, "File could not be opened, program will now terminate!\nPath: %s\n", file);
		return NULL;
	}

	// open the dataset and attributes
	position_set	= H5Dopen2(file_id, "/positions", H5P_DEFAULT);
	attr_N_id		= H5Aopen(position_set, "N", H5P_DEFAULT);
	attr_write_id	= H5Aopen(position_set, "Writeouts", H5P_DEFAULT);
	attr_L_x_id		= H5Aopen(position_set, "L_x", H5P_DEFAULT);
	attr_L_y_id		= H5Aopen(position_set, "L_y", H5P_DEFAULT);
	attr_X_id 		= H5Aopen(position_set, "X_A", H5P_DEFAULT);

	// read the attributes
	status			= H5Aread(attr_N_id, H5T_NATIVE_INT, &N);
	status			= H5Aread(attr_write_id, H5T_NATIVE_INT, &steps);
	status 			= H5Aread(attr_L_x_id, H5T_NATIVE_DOUBLE, &L_x);
	status 			= H5Aread(attr_L_y_id, H5T_NATIVE_DOUBLE, &L_y);
	status 			= H5Aread(attr_X_id, H5T_NATIVE_DOUBLE, &X);

	// close the attributes
	status 			= H5Aclose(attr_N_id);
	status 			= H5Aclose(attr_write_id);
	status 			= H5Aclose(attr_L_x_id);
	status 			= H5Aclose(attr_L_y_id);
	status 			= H5Aclose(attr_X_id);

	// set dimensions for all data
	offset[0]		= 0;
	offset[1]		= 0;
	slabdim[0]		= 1;
	slabdim[1] 		= 2*N;

	// get the data space of the current dataset and initialize a buffer for reading data
	position_space	= H5Dget_space(position_set);
	tempset_id 		= H5Screate_simple(2, slabdim, NULL);

	// allocate memory for data
	positions 		= malloc(2*N*sizeof(double));

	// create a struct from all known parameters
	struct parameters *param = malloc(sizeof(struct parameters));

	param->N 			= N;
	param->steps 		= steps;

	param->L_x			= L_x;
	param->L_y			= L_y;
	param->X 			= X;

	param->positions 	= positions;

	// return to caller
	return param;
}