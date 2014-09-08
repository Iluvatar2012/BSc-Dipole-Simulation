#include <stdlib.h>
#include <stdio.h>

#include <hdf5.h>

#include "structs.h"


// basic variables
static double* 	positions;
static double* 	psi4;
static double* 	psi6;

static int 		N;
static int 		steps;

// identifiers for files, dataspaces and -sets
static hid_t	file_id; 
static hid_t	position_set;
static hid_t	tempset_id;
static hid_t	psi4_set;
static hid_t	psi6_set;
static hid_t	psi_tempset_id;
static hid_t	position_space;
static hid_t	psi4_space;
static hid_t	psi6_space;

// data dimensions
static hsize_t	offset[2];
static hsize_t 	slabdim[2];
static hsize_t 	psi_offset[1];
static hsize_t	psi_dim[1];


/*----------------------------------------------------------------------------------------------------------------------------*/
void hdf5_read (int curr_step) {
	// checking variable
	herr_t	status;

	// adjust the current offset
	offset[0] 		= curr_step;
	psi_offset[0]	= N*curr_step;

	// select the hyperslab and read into position array
	status 	= H5Sselect_hyperslab(position_space, H5S_SELECT_SET, offset, NULL, slabdim, NULL);
	status 	= H5Dread(position_set, H5T_NATIVE_DOUBLE, tempset_id, position_space, H5P_DEFAULT, positions);
	
	// select the hyperslab and read into psi4 array
	status	= H5Sselect_hyperslab(psi4_space, H5S_SELECT_SET, psi_offset, NULL, psi_dim, NULL);
	status	= H5Dread(psi4_set, H5T_NATIVE_DOUBLE, psi_tempset_id, psi4_space, H5P_DEFAULT, psi4);

	// select the hyperslab and read into psi6 array
	status	= H5Sselect_hyperslab(psi6_space, H5S_SELECT_SET, psi_offset, NULL, psi_dim, NULL);
	status	= H5Dread(psi6_set, H5T_NATIVE_DOUBLE, psi_tempset_id, psi6_space, H5P_DEFAULT, psi6);
}


/*----------------------------------------------------------------------------------------------------------------------------*/
struct parameters *hdf5_init (char* file) {


	hid_t	attr_write_id, attr_N_id;
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

	// read the attributes
	status			= H5Aread(attr_N_id, H5T_NATIVE_INT, &N);
	status			= H5Aread(attr_write_id, H5T_NATIVE_INT, &steps);

	// close the attributes
	status 			= H5Aclose(attr_N_id);
	status 			= H5Aclose(attr_write_id);

	// set dimensions for all data
	offset[0]		= 0;
	offset[1]		= 0;
	slabdim[0]		= 1;
	slabdim[1] 		= 2*N;

	psi_offset[0]	= 0;
	psi_dim[0]		= N;

	// get the data space of the current dataset and initialize a buffer for reading data
	position_space	= H5Dget_space(position_set);
	tempset_id 		= H5Screate_simple(2, slabdim, NULL);

	// open datasets for psi values
	psi4_set	= H5Dopen2(file_id, "/psi4", H5P_DEFAULT);
	psi6_set	= H5Dopen2(file_id, "/psi6", H5P_DEFAULT);

	// get the data space of the current psiset and initialize a buffer for reading data
	psi4_space		= H5Dget_space(psi4_set);
	psi6_space		= H5Dget_space(psi6_set);
	psi_tempset_id	= H5Screate_simple(1, psi_dim, NULL);

	// allocate memory for data
	positions 		= malloc(2*N*sizeof(double));
	psi4 			= malloc(N*sizeof(double));
	psi6 			= malloc(N*sizeof(double));

	// create a struct from all known parameters
	struct parameters *param = malloc(sizeof(struct parameters));

	param->N 			= N;
	param->steps 		= steps;

	param->positions 	= positions;
	param->psi4 		= psi4;
	param->psi6 		= psi6;

	// return to caller
	return param;
}