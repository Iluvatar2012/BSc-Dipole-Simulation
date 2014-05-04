#include <stdlib.h>
#include <stdio.h>

#include <hdf5.h>

#include "structs.h"


/*----------------------------------------------------------------------------------------------------------------------------*/
struct parameters* hdf5_read(char* file) {

	// basic variables
	int 		N, steps;

	// variables for holding return parameters
	double* 	positions;
	double* 	displacement;

	double* 	laning;
	double*		psi4;
	double* 	psi6;

	double*		temp;

	// identifiers for files, dataspaces and -sets
	hid_t	file_id, pos_set, disp_set, tempset_id, psi4_set, psi6_set, laning_set;
	hid_t	pos_space, disp_space, attr_write_id, attr_N_id;
	herr_t	status;

	// open file, check if successful
	file_id = H5Fopen(file, H5F_ACC_RDONLY, H5P_DEFAULT);

	if (file_id < 0) {
		fprintf(stderr, "File could not be opened, program will now terminate!\nPath: %s\n", file);
		return NULL;
	}

	// open the dataset and attributes
	pos_set			= H5Dopen2(file_id, "/positions", H5P_DEFAULT);
	disp_set 		= H5Dopen2(file_id, "/displacement", H5P_DEFAULT);
	attr_N_id		= H5Aopen(pos_set, "N", H5P_DEFAULT);
	attr_write_id	= H5Aopen(pos_set, "Writeouts", H5P_DEFAULT);

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
	displacement	= malloc((steps+1)*2*N*sizeof(double));

	psi4 			= malloc((steps+1)*N*sizeof(double));
	psi6 			= malloc((steps+1)*N*sizeof(double));
	laning 			= malloc((steps+1)*N*sizeof(double));

	temp 			= malloc(2*N*sizeof(double));

	// open datasets for psi values
	psi4_set	= H5Dopen2(file_id, "/psi4", H5P_DEFAULT);
	psi6_set	= H5Dopen2(file_id, "/psi6", H5P_DEFAULT);
	laning_set	= H5Dopen2(file_id, "/laning", H5P_DEFAULT);

	// read psi and laning values, check if successful
	status		= H5Dread(psi4_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, psi4);

	if (status != 0) {
		fprintf(stderr, "Psi 4 values could not be read, program will now terminate!\n");
		return NULL;
	}

	status		= H5Dread(psi6_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, psi6);

	if (status != 0) {
		fprintf(stderr, "Psi 6 values could not be read, program will now terminate!\n");
		return NULL;
	}

	status		= H5Dread(laning_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, laning);

	if (status != 0) {
		fprintf(stderr, "Laning values could not be read, program will now terminate!\n");
		return NULL;
	}


	// close the rest of the file
	status 		= H5Dclose(psi4_set);
	status		= H5Dclose(psi6_set);
	status		= H5Dclose(laning_set);

	// get the data space of the current position dataset and initialize a buffer for reading data
	pos_space	= H5Dget_space(pos_set);
	tempset_id 	= H5Screate_simple(2, slabdim, NULL);

	// iterate over all steps and copy them all into the position array
	for (int i=0; i<=steps; i++) {

		// adjust offset to the next slab
		offset[0] = i;

		// select the hyperslab and read into position array
		status 		= H5Sselect_hyperslab(pos_space, H5S_SELECT_SET, offset, NULL, slabdim, NULL);
		status 		= H5Dread(pos_set, H5T_NATIVE_DOUBLE, tempset_id, pos_space, H5P_DEFAULT, temp);

		// copy the data from the buffer to the position array
		for (int j=0; j<2*N; j+=2) {
			positions[2*N*i + j] 	= temp[j];
			positions[2*N*i + j+1]	= temp[j+1];
		}
	}

	// close these parts of the file, free memory not used any more
	status 		= H5Sclose(pos_space);
	status 		= H5Dclose(pos_set);

	// get the data space of the current displacement dataset
	disp_space	= H5Dget_space(disp_set);

	// iterate over all steps and copy them all into the displacement array
	for (int i=0; i<=steps; i++) {

		// adjust offset to the next slab
		offset[0] = i;

		// select the hyperslab and read into displacement array
		status 		= H5Sselect_hyperslab(disp_space, H5S_SELECT_SET, offset, NULL, slabdim, NULL);
		status 		= H5Dread(disp_set, H5T_NATIVE_DOUBLE, tempset_id, disp_space, H5P_DEFAULT, temp);

		// copy the data from the buffer to the displacement array
		for (int j=0; j<2*N; j+=2) {
			displacement[2*N*i + j] 	= temp[j];
			displacement[2*N*i + j+1]	= temp[j+1];
		}
	}

	// close these parts of the file, free memory not used any more
	status 		= H5Sclose(tempset_id);
	status 		= H5Sclose(disp_space);
	status 		= H5Dclose(disp_set);

	free(temp);

	// create a struct from all known parameters
	struct parameters *param = malloc(sizeof(struct parameters));

	param->N 			= N;
	param->steps 		= steps;

	param->positions 	= positions;
	param->displacement = displacement;

	param->psi4 		= psi4;
	param->psi6 		= psi6;
	param->laning 		= laning;

	// return to caller
	return param;
}