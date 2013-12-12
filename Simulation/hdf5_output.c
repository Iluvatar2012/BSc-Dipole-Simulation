/*
 *hdf5_write.c
 *
 *  Created on: Oct 17, 2013
 *      Author: Aiko Bernehed
 */

#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <unistd.h>

// string for the filename used by various methods
static char* 	file;

// dimensions of the dataset and of one configuration
static hsize_t	dims[2];
static hsize_t	chunkdims[2];

// various variables
static int 	counter;
static int 	N;










/*-------------------------------------------------------------------------------------------------------*/
int create_file (char* infile, int inN, int writeouts) {

	// save the name of the file, this way we always use the same file throughout one run, initiate counter
	file 		= infile;
	N			= inN;
	counter		= -1;

	// identifiers of all relevant spaces and chunking properties
	hid_t	file_id, dataspace_id, dataset_id;
	hid_t	attr_space_id_1, attr_space_id_2, attr_space_id_3, attr_id_1, attr_id_2, attr_id_3;
	hid_t	property;

	// checking value whether operations were successful
	herr_t	status;

	// set the size of all dimensions
	dims[0] 		= 0;
	dims[1]			= 2*N;

	chunkdims[0]	= 1;
	chunkdims[1]	= 2*N;

	// set maximum dimension size and attribute size
	hsize_t maxdims[2] 		= {H5S_UNLIMITED, 2*N};
	hsize_t attr_dim		= 1;

	// check whether file already exists, if not create new file using the identifier
	if (access(file, F_OK) != -1)
		return EXIT_FAILURE;

	file_id = H5Fcreate(file, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);

	// create a new dataspace for the position data of all particles
	dataspace_id 	= H5Screate_simple(2, dims, maxdims);

	// create a property list enabling chunking
	property		= H5Pcreate(H5P_DATASET_CREATE);
	status			= H5Pset_chunk(property, 2, chunkdims);

	// create a new dataset for the positions of all particles, at this point the dataset will have dimension zero
	dataset_id		= H5Dcreate2(file_id, "/positions", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, property, H5P_DEFAULT);

	// create new dataspaces for the attributes (writeouts, N, amount of data written, and name of the config file) and the attributes themselves
	attr_space_id_1	= H5Screate_simple(1, &attr_dim, NULL);
	attr_space_id_2 = H5Screate_simple(1, &attr_dim, NULL);
	attr_space_id_3 = H5Screate_simple(1, &attr_dim, NULL);
	attr_id_1		= H5Acreate2(dataset_id, "Writeouts", H5T_NATIVE_INT, attr_space_id_1, H5P_DEFAULT, H5P_DEFAULT);
	attr_id_2		= H5Acreate2(dataset_id, "N", H5T_NATIVE_INT, attr_space_id_2, H5P_DEFAULT, H5P_DEFAULT);
	attr_id_3		= H5Acreate2(dataset_id, "Last_Writeout", H5T_NATIVE_INT, attr_space_id_3, H5P_DEFAULT, H5P_DEFAULT);

	// write the attributes
	status = H5Awrite(attr_id_1, H5T_NATIVE_INT, &writeouts);
	status = H5Awrite(attr_id_2, H5T_NATIVE_INT, &N);
	status = H5Awrite(attr_id_3, H5T_NATIVE_INT, &counter);

	// terminate access to all parts of the document and close the file
	status = H5Aclose(attr_id_1);
	status = H5Aclose(attr_id_2);
	status = H5Aclose(attr_id_3);
	status = H5Sclose(attr_space_id_1);
	status = H5Sclose(attr_space_id_2);
	status = H5Sclose(attr_space_id_3);

	status = H5Pclose(property);
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
void write_data(int timestep, double* position) {
	// variables for IDs of the file, dataset and attribute
	hid_t	file_id, dataset_id, attr_id;
	hid_t	tempset_id, dataspace_id;

	// default error variable
	herr_t	status;

	// offset, needed in order to select the next hyperslab within a dataset
	hsize_t offset[2] = {dims[0], 0};

	// open an existing file, check whether operation was successful
	if ((file_id = H5Fopen(file, H5F_ACC_RDWR, H5P_DEFAULT)) < 0) {
		fprintf(stderr, "File could not be opened on timestep %d, program will now terminate!\n", timestep);
		exit(1);
	}

	// open an existing dataset
	dataset_id = H5Dopen2(file_id, "/positions", H5P_DEFAULT);

	// increment the amount of positions already written and increase the size of the dataset
	dims[0]++;
	H5Dset_extent(dataset_id, dims);

	// create a temporary dataspace, which will be used to write the current positions to the file
	tempset_id = H5Screate_simple(2, chunkdims, NULL);

	// get the current dataspace and select one hyperslab
	dataspace_id 	= H5Dget_space(dataset_id);
	status			= H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, chunkdims, NULL);

	// write current positions to the file, check whether operation was successful
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, tempset_id, dataspace_id, H5P_DEFAULT, position);

	if (status < 0) {
		fprintf(stderr, "Writing of configuration data failed on timestep: %d, program will now terminate.\n", timestep);
		exit(1);
	}

	// open an attribute, increase the current counter and write the counter to the attribute
	counter++;
	attr_id = H5Aopen(dataset_id, "Last_Writeout", H5P_DEFAULT);
	status 	= H5Awrite(attr_id, H5T_NATIVE_INT, &counter);


	// terminate access to all parts of the file and close the latter
	status = H5Sclose(tempset_id);
	status = H5Aclose(attr_id);
	status = H5Sclose(dataspace_id);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);
}










/*-------------------------------------------------------------------------------------------------------*/
int reopen_file(double* position, int writeouts, int* count_written) {

	// various identifiers
	hid_t	file_id, dataset_id, attr_id_one, attr_id_two;
	hid_t	dataspace_id, tempspace_id;

	// default error variable
	herr_t 	status;

	// open file, dataset and an attribute
	file_id 	= H5Fopen(file, H5F_ACC_RDWR, H5P_DEFAULT);
	dataset_id 	= H5Dopen2(file_id, "/positions", H5P_DEFAULT);
	attr_id_one	= H5Aopen(dataset_id, "Last_Writeout", H5P_DEFAULT);
	attr_id_two = H5Aopen(dataset_id, "Writeouts", H5P_DEFAULT);

	// copy the attribute into counter variable
	status	= H5Aread(attr_id_one, H5T_NATIVE_INT, &counter);

	// check whether the file has been fully written, terminate if it has completed
	if (counter == writeouts) {
		fprintf(stderr, "Simulation is complete, program will terminate. \n");

		status = H5Aclose(attr_id_one);
		status = H5Aclose(attr_id_two);
		status = H5Dclose(dataset_id);
		status = H5Fclose(file_id);

		return EXIT_FAILURE;
	}

	// rewrite attribute with new writeout value, could be that the file has been reopened
	status = H5Awrite(attr_id_two, H5T_NATIVE_INT, &writeouts);

	// variables for offset and size of data copied
	hsize_t offset[2] 	= {counter, 0};
	hsize_t slabdim[2] 	= {1, 2*N};

	// select a hyperslab of the last positions stored in the file and save it to the position array
	tempspace_id = H5Screate_simple(2, slabdim, NULL);
	dataspace_id = H5Dget_space(dataset_id);

	status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, slabdim, NULL);
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, tempspace_id, dataspace_id, H5P_DEFAULT, position);

	// copy how many iterations have already been written to the file
	*count_written 	= counter;
	dims[0]			= counter+1;

	// terminate all access to and close the file
	status = H5Aclose(attr_id_one);
	status = H5Aclose(attr_id_two);
	status = H5Dclose(dataset_id);
	status = H5Sclose(tempspace_id);
	status = H5Sclose(dataspace_id);
	status = H5Fclose(file_id);

	// return to caller
	return EXIT_SUCCESS;
}
