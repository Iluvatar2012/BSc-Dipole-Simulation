#include <stdlib.h>
#include <stdio.h>

#include <hdf5.h>

#include "structs.h"


/*----------------------------------------------------------------------------------------------------------------------------*/
struct parameters *hdf5_read (char* file) {

	// basic variables
	double* positions;
	double*	temp;
	double*	psi4;
	double* psi6;

	int 	N, steps;

	// identifiers for files, dataspaces and -sets
	hid_t	file_id, position_set, tempset_id, psi4_set, psi6_set;
	hid_t	position_space, attr_write_id, attr_N_id;
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

	// set dimensions for remainder of data
	hsize_t	offset[2]	= {0, 0};
	hsize_t slabdim[2]	= {1, 2*N};

	// allocate memory for data
	positions 		= malloc((steps+1)*2*N*sizeof(double));
	psi4 			= malloc((steps+1)*N*sizeof(double));
	psi6 			= malloc((steps+1)*N*sizeof(double));

	temp 			= malloc(2*N*sizeof(double));


	// get the data space of the current dataset and initialize a buffer for reading data
	position_space	= H5Dget_space(position_set);
	tempset_id 		= H5Screate_simple(2, slabdim, NULL);

	// iterate over all steps and copy them all into the position array
	for (int i=0; i<=steps; i++) {

		// adjust offset to the next slab
		offset[0] = i;

		// select the hyperslab and read into position array
		status 		= H5Sselect_hyperslab(position_space, H5S_SELECT_SET, offset, NULL, slabdim, NULL);
		status 		= H5Dread(position_set, H5T_NATIVE_DOUBLE, tempset_id, position_space, H5P_DEFAULT, temp);

		// copy the data from the buffer to the position array
		for (int j=0; j<2*N; j+=2) {
			positions[2*N*i + j] 	= temp[j];
			positions[2*N*i + j+1]	= temp[j+1];
		}
	}

	// close these parts of the file, free memory not used any more
	status 		= H5Sclose(tempset_id);
	status 		= H5Sclose(position_space);
	status 		= H5Dclose(position_set);

	free(temp);

	// open datasets for psi values
	psi4_set	= H5Dopen2(file_id, "/psi4", H5P_DEFAULT);
	psi6_set	= H5Dopen2(file_id, "/psi6", H5P_DEFAULT);

	// read psi values from file
	status		= H5Dread(psi4_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, psi4);
	status		= H5Dread(psi6_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, psi6);

	// close the rest of the file
	status 		= H5Dclose(psi4_set);
	status		= H5Dclose(psi6_set);
	status		= H5Fclose(file_id);

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