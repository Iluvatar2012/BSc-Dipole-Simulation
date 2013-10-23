/*
 *hdf5_write.c
 *
 *  Created on: Oct 17, 2013
 *      Author: Aiko Bernehed
 */

#include <stdio.h>
#include <stdlib.h>
#include <hdf5/hdf5.h>

/*-------------------------------------------------------------------------------------------------------*/
int create_file (char* file, int N) {

	// identifiers of all relevant spaces
	hid_t	file_id, dataspace_id, dataset_id, attr_space_id, attr_id;

	// identifier for chunking property
	hid_t	property;

	// checking value whether operations were successful
	herr_t	status;

	// size of all dimensions
	hsize_t	dims[2] 		= {0,2*N};
	hsize_t maxdims[2] 		= {H5S_UNLIMITED, 2*N};
	hsize_t chunkdims[2]	= {1, 2*N};

	// size of the attribute
	hsize_t attr_dim		= H5S_UNLIMITED;

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

	// TODO: create attributes and close all streams

	// return to caller
	return EXIT_SUCCESS;
}

/*-------------------------------------------------------------------------------------------------------*/
