#include <hdf5.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>



/*----------------------------------------------------------------------------------------------------------------------------*/
double* read_config (char* file, char* number, int* N) {

	// Basic variables
	int writes;
	int step;

	// return value
	double* config;

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
	status			= H5Aread(attr_N_id, H5T_NATIVE_INT, N);
	status			= H5Aread(attr_write_id, H5T_NATIVE_INT, &writes);

	// try to get memory from the system
	config = malloc(2* *N *sizeof(double));

	// check which configuration to open from file
	if 		(strncmp(number, "first", 5) == 0)
		step = 0;
	else if (strncmp(number, "last", 4) == 0)
		step = writes;
	else
		step = atoi(number);

		// offset and dimension of a hyperslab
	hsize_t	offset[2]	= {step, 0};
	hsize_t	slabdim[2]	= {1, 2* *N};

	// get the data space of the current dataset and initialize a buffer for reading data
	dataspace_id	= H5Dget_space(dataset_id);
	tempset_id 		= H5Screate_simple(2, slabdim, NULL);

	// select the hyperslab and read into position array
	status 		= H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, slabdim, NULL);
	status 		= H5Dread(dataset_id, H5T_NATIVE_DOUBLE, tempset_id, dataspace_id, H5P_DEFAULT, config);
	
	// close all parts of the file
	status = H5Aclose(attr_N_id);
	status = H5Aclose(attr_write_id);
	status = H5Sclose(tempset_id);
	status = H5Sclose(dataspace_id);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);

	// return to caller
	return config;
}