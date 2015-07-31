#include <stdlib.h>
#include <stdio.h>

#include <hdf5.h>

#include "structs.h"

// basic variables
static int 		N;
static int 		steps;
static int 		last_step;

static double 	X;
static double 	L_y;

// variables for all arrays
static double* positions;

// dimensions of the 
static hsize_t posdim[2];
static hsize_t pos_offset[2];

// identifiers for files, dataspaces and -sets
static hid_t file_id;
static hid_t position_set;
static hid_t tempset_id;
static hid_t position_space;
static hid_t attr_write_id;
static hid_t attr_N_id;
static hid_t attr_X_id;
static hid_t attr_L_y_id;
static hid_t attr_last_id;

// status variable
static herr_t status;

/*----------------------------------------------------------------------------------------------------------------------------*/
void hdf5_read (int curr_step) {

	// adjust offset to the next slab
	pos_offset[0] = curr_step;

	// select the hyperslab and read into position array
	status 	= H5Sselect_hyperslab(position_space, H5S_SELECT_SET, pos_offset, NULL, posdim, NULL);
	status 	= H5Dread(position_set, H5T_NATIVE_DOUBLE, tempset_id, position_space, H5P_DEFAULT, positions);
}



/*----------------------------------------------------------------------------------------------------------------------------*/
struct parameters *hdf5_init (char* file) {
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
	attr_X_id		= H5Aopen(position_set, "X_A", H5P_DEFAULT);
	attr_L_y_id		= H5Aopen(position_set, "L_y", H5P_DEFAULT);
	attr_last_id	= H5Aopen(position_set, "Last_Writeout", H5P_DEFAULT);

	// read the attributes
	status			= H5Aread(attr_N_id, H5T_NATIVE_INT, &N);
	status			= H5Aread(attr_write_id, H5T_NATIVE_INT, &steps);
	status			= H5Aread(attr_X_id, H5T_NATIVE_DOUBLE, &X);
	status			= H5Aread(attr_L_y_id, H5T_NATIVE_DOUBLE, &L_y);
	status			= H5Aread(attr_last_id, H5T_NATIVE_INT, &last_step);

	// close the attributes
	status 			= H5Aclose(attr_N_id);
	status 			= H5Aclose(attr_write_id);
	status 			= H5Aclose(attr_X_id);
	status 			= H5Aclose(attr_L_y_id);
	status 			= H5Aclose(attr_last_id);

	// set dimensions for remainder of data
	pos_offset[0]	= 0;
	pos_offset[1]	= 0;
	posdim[0]		= 1;
	posdim[1]		= 2*N;

	// allocate memory for data
	positions 		= malloc(2*N*sizeof(double));

	// get the data space of the current dataset and initialize a buffer for reading data
	position_space	= H5Dget_space(position_set);
	tempset_id 		= H5Screate_simple(2, posdim, NULL);

	// create a struct from all known parameters
	struct parameters *param = malloc(sizeof(struct parameters));

	param->N 			= N;
	param->steps 		= steps;
	param->last_step	= last_step;

	param->X 			= X;
	param->L_y 			= L_y;

	param->positions 	= positions;

	// return to caller
	return param;
}