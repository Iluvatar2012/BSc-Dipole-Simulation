#include <stdlib.h>
#include <stdio.h>

#include <hdf5.h>

#include "structs.h"

/*----------------------------------------------------------------------------------------------------------------------------*/
struct parameters *hdf5_read (char* file) {

	// basic variables
	double* temp;
	double* positions;

	int N;
	int steps;

	// identifiers for files and a status variable
	hid_t	file_id, dataset_id, tempset_id, attr_write_id, attr_N_id;
	hid_t	dataspace_id;
	herr_t	status;

	// open the file, check whether operation was successful
	file_id			= H5Fopen(file, H5F_ACC_RDONLY, H5P_DEFAULT);

	if (file_id < 0) {
		fprintf(stderr, "File could not be opened, program will now terminate!\nPath: %s\n", file);
		return NULL;
	}

	// open the dataset and attributes
	dataset_id		= H5Dopen2(file_id, "/positions", H5P_DEFAULT);
	attr_N_id		= H5Aopen(dataset_id, "N", H5P_DEFAULT);
	attr_write_id	= H5Aopen(dataset_id, "Writeouts", H5P_DEFAULT);

	// read the attributes
	status			= H5Aread(attr_N_id, H5T_NATIVE_INT, &N);
	status			= H5Aread(attr_write_id, H5T_NATIVE_INT, &steps);

	// allocate memory for the temporary data and positions
	temp			= malloc(2*N*sizeof(double));
	positions 		= malloc((steps+1)*2*N*sizeof(double));

	// offset and dimension of a hyperslab
	hsize_t	offset[2]	= {0, 0};
	hsize_t	slabdim[2]	= {1, 2*N};

	// get the data space of the current dataset and initialize a buffer for reading data
	dataspace_id	= H5Dget_space(dataset_id);
	tempset_id 		= H5Screate_simple(2, slabdim, NULL);

	// iterate over all steps and copy them all into the position array
	for (int i=0; i<=steps; i++) {

		// adjust offset to the next slab
		offset[0] = i;

		// select the hyperslab and read into position array
		status 		= H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, slabdim, NULL);
		status 		= H5Dread(dataset_id, H5T_NATIVE_DOUBLE, tempset_id, dataspace_id, H5P_DEFAULT, temp);

		// copy the data from the buffer to the position array
		for (int j=0; j<2*N; j+=2) {
			positions[2*N*i + j] 	= temp[j];
			positions[2*N*i + j+1]	= temp[j+1];
		}
	}

	// free memory
	free(temp);

	// build a new struct to store parameters in, copy variables to struct
	struct parameters *param = malloc(sizeof(struct parameters));

	param->positions 	= positions;
	param->N 			= N;
	param->steps 		= steps;

	// close all parts of the simulation
	status = H5Aclose(attr_N_id);
	status = H5Aclose(attr_write_id);
	status = H5Sclose(tempset_id);
	status = H5Sclose(dataspace_id);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);

	// return to caller
	return param;
}


/*----------------------------------------------------------------------------------------------------------------------------*/
int add_analysis (struct analysis* data) {

	// get all data from incoming struct
	int N 			= data->N;
	int steps 		= data->steps;

	char* file 		= data->file;

	double* psi4 	= data->psi4;
	double* psi6 	= data->psi6;
	double* laning 	= data->laning;	

	// identfiers
	hid_t	file_id, psi4_set, psi6_set, laning_set;
	hid_t	psi4_space, psi6_space, laning_space;
	herr_t	status;

	// dimensions
	hsize_t	dims		= (steps+1)*N;

	// open the file, check whether operation was successful
	file_id	= H5Fopen(file, H5F_ACC_RDWR, H5P_DEFAULT);

	if (file_id < 0) {
		fprintf(stderr, "File could not be opened, program will now terminate!\nPath: %s\n", file);
		return EXIT_FAILURE;
	}

	// create new dataspaces for data analysis
	psi4_space 		= H5Screate_simple(1, &dims, NULL);
	psi6_space 		= H5Screate_simple(1, &dims, NULL);
	laning_space 	= H5Screate_simple(1, &dims, NULL);

	// create new datasets for all data
	psi4_set	= H5Dcreate2(file_id, "/psi4", H5T_NATIVE_DOUBLE, psi4_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	psi6_set	= H5Dcreate2(file_id, "/psi6", H5T_NATIVE_DOUBLE, psi6_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	laning_set	= H5Dcreate2(file_id, "/laning", H5T_NATIVE_DOUBLE, laning_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// write all data
	status = H5Dwrite(psi4_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, psi4);
	status = H5Dwrite(psi6_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, psi6);
	status = H5Dwrite(laning_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, laning);

	// close everything, free struct
	H5Dclose(psi4_set);
	H5Dclose(psi6_set);
	H5Dclose(laning_set);
	H5Sclose(psi4_space);
	H5Sclose(psi6_space);
	H5Sclose(laning_space);
	H5Fclose(file_id);

	// return to caller
	return EXIT_SUCCESS;
}