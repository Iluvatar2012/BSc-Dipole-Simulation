#include <stdlib.h>
#include <stdio.h>

#include <hdf5.h>

#include "structs.h"

// basic variables
static int N;
static int steps;

// variables for all arrays
static double* positions;
static double* disp;
static double* psi4;
static double* psi6;
static double* laning;

// dimensions of the 
static hsize_t posdim[2];
static hsize_t pos_offset[2];
static hsize_t psi_dim[1];
static hsize_t psi_offset[1];

// identifiers for files, dataspaces and -sets
static hid_t file_id;
static hid_t position_set;
static hid_t disp_set;
static hid_t tempset_id;
static hid_t psi4_set;
static hid_t psi6_set;
static hid_t lane_set;
static hid_t position_space;
static hid_t disp_space;
static hid_t psi4_space;
static hid_t psi6_space;
static hid_t lane_space;
static hid_t psi_temp_id;
static hid_t attr_write_id;
static hid_t attr_N_id;

// status variable
static herr_t status;

/*----------------------------------------------------------------------------------------------------------------------------*/
void hdf5_read () {
	// select the hyperslab and read into position array
	status 	= H5Sselect_hyperslab(position_space, H5S_SELECT_SET, pos_offset, NULL, posdim, NULL);
	status 	= H5Dread(position_set, H5T_NATIVE_DOUBLE, tempset_id, position_space, H5P_DEFAULT, positions);

	// select the hyperslab and read into displacement array
	status 	= H5Sselect_hyperslab(disp_space, H5S_SELECT_SET, pos_offset, NULL, posdim, NULL);
	status 	= H5Dread(disp_set, H5T_NATIVE_DOUBLE, tempset_id, disp_space, H5P_DEFAULT, disp);
	
	// select the hyperslab and read into psi4 array
	status	= H5Sselect_hyperslab(psi4_space, H5S_SELECT_SET, psi_offset, NULL, psi_dim, NULL);
	status	= H5Dread(psi4_set, H5T_NATIVE_DOUBLE, psi_temp_id, psi4_space, H5P_DEFAULT, psi4);

	// select the hyperslab and read into psi6 array
	status	= H5Sselect_hyperslab(psi6_space, H5S_SELECT_SET, psi_offset, NULL, psi_dim, NULL);
	status	= H5Dread(psi6_set, H5T_NATIVE_DOUBLE, psi_temp_id, psi6_space, H5P_DEFAULT, psi6);

	// select the hyperslab and read into laning array
	status	= H5Sselect_hyperslab(lane_space, H5S_SELECT_SET, psi_offset, NULL, psi_dim, NULL);
	status	= H5Dread(lane_set, H5T_NATIVE_DOUBLE, psi_temp_id, lane_space, H5P_DEFAULT, laning);

	// adjust offset to the next slab
	pos_offset[0]++;
	psi_offset[0] += N;
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
	disp_set		= H5Dopen2(file_id, "/displacement", H5P_DEFAULT);
	attr_N_id		= H5Aopen(position_set, "N", H5P_DEFAULT);
	attr_write_id	= H5Aopen(position_set, "Writeouts", H5P_DEFAULT);

	// read the attributes
	status			= H5Aread(attr_N_id, H5T_NATIVE_INT, &N);
	status			= H5Aread(attr_write_id, H5T_NATIVE_INT, &steps);

	// close the attributes
	status 			= H5Aclose(attr_N_id);
	status 			= H5Aclose(attr_write_id);

	// set dimensions for remainder of data
	pos_offset[0]	= 0;
	pos_offset[1]	= 0;
	posdim[0]		= 1;
	posdim[1]		= 2*N;

	psi_offset[0]	= 0;
	psi_dim[0]		= N;

	// allocate memory for data
	positions 		= malloc(2*N*sizeof(double));
	disp 			= malloc(2*N*sizeof(double));
	psi4 			= malloc(N*sizeof(double));
	psi6 			= malloc(N*sizeof(double));
	laning			= malloc(N*sizeof(double));

	// get the data space of the current dataset and initialize a buffer for reading data
	position_space	= H5Dget_space(position_set);
	disp_space		= H5Dget_space(disp_set);
	tempset_id 		= H5Screate_simple(2, posdim, NULL);

	// open datasets for psi values
	psi4_set	= H5Dopen2(file_id, "/psi4", H5P_DEFAULT);
	psi6_set	= H5Dopen2(file_id, "/psi6", H5P_DEFAULT);
	lane_set 	= H5Dopen2(file_id, "/laning", H5P_DEFAULT);

	// get the dataspaces of the previous sets
	psi4_space 	= H5Dget_space(psi4_set);
	psi6_space 	= H5Dget_space(psi6_set);
	lane_space 	= H5Dget_space(lane_set);

	// initialize a buffer for reading data
	psi_temp_id = H5Screate_simple(1, psi_dim, NULL);

	// create a struct from all known parameters
	struct parameters *param = malloc(sizeof(struct parameters));

	param->N 			= N;
	param->steps 		= steps;

	param->positions 	= positions;
	param->displacement = disp;
	param->psi4 		= psi4;
	param->psi6 		= psi6;
	param->laning		= laning;

	// return to caller
	return param;
}