/*
 *hdf5_write.c
 *
 *  Created on: Oct 17, 2013
 *      Author: Aiko Bernehed
 */

#include <stdio.h>
#include <stdlib.h>
#include <hdf5/hdf5.h>

// string for the filename used by various methods
static char* 	file;

// dimensions of the dataset
static hsize_t	dims[2];

/*-------------------------------------------------------------------------------------------------------*/
int create_file (char* infile, int N) {

	// save the name of the file, this way we always use the same file throughout one run
	file = infile;

	// identifiers of all relevant spaces
	hid_t	file_id, dataspace_id, dataset_id, attr_space_id, attr_id;

	// identifier for chunking property
	hid_t	property;

	// checking value whether operations were successful
	herr_t	status;

	// set the size of all dimensions
	dims[0] = 0;
	dims[1]	= 2*N;

	// set maximum dimension size and the size of one datachunk
	hsize_t maxdims[2] 		= {H5S_UNLIMITED, 2*N};
	hsize_t chunkdims[2]	= {1, 2*N};

	// size of the attribute
	hsize_t attr_dim		= 0;
	hsize_t max_attr_dim	= H5S_UNLIMITED;

	// create new file using the identifier, check whether file has previously been opened
	if ((file_id = H5Fcreate(file, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT)) < 0 ) {
		fprintf(stderr, "File already exists, function will exit!\n");
		return EXIT_FAILURE;
	}

	// create a new dataspace for the position data of all particles
	dataspace_id 	= H5Screate_simple(2, dims, maxdims);

	// create a property list enabling chunking
	property		= H5Pcreate(H5P_DATASET_CREATE);
	status			= H5Pset_chunk(property, 2, chunkdims);

	// create a new dataset for the positions of all particles, at this point the dataset will have dimension zero
	dataset_id		= H5Dcreate2(file_id, "/positions", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, property, H5P_DEFAULT);

	// create a new dataspace for the attribute (timestep) and the attribute itself
	attr_space_id	= H5Screate_simple(1, &attr_dim, &max_attr_dim);
	attr_id			= H5Acreate2(dataset_id, "Timestep", H5T_NATIVE_INT, attr_space_id, H5P_DEFAULT, H5P_DEFAULT);

	// terminate access to all parts of the document and close the file
	status = H5Aclose(attr_id);
	status = H5Sclose(attr_space_id);
	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);
	status = H5Fclose(file_id);

	// check whether file could be closed
	if (status <0)
		fprintf(stderr, "Error closing the file, execution will continue. \n");

	// return to caller
	return EXIT_SUCCESS;
}

/*-------------------------------------------------------------------------------------------------------*/
int write_data(int timestep, double* position) {

	// variables for IDs of the file, dataset and attribute
	hid_t	file_id, dataset_id, attr_id;

	// offset in order to select the next hyperslab within a dataset
	hsize_t offset[2] = {dims[0], 0};

	// open an existing file, check whether operation was successful
	if ((file_id = H5Fopen(file, H5F_ACC_RDWR, H5P_DEFAULT)) < 0) {
		fprintf(stderr, "File could not be opened on timestep %d, program will now terminate!\n", timestep);
		return EXIT_FAILURE;
	}

	// open an existing dataset
	dataset_id = H5Dopen2(file_id, "/positions", H5P_DEFAULT);

	// increment the amount of positions already written and increase the size of the dataset
	dims[0]++;
	H5Dset_extent(dataset_id, dims);


	// return to caller
	return EXIT_SUCCESS;
}
