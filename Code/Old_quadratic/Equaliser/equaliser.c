#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <hdf5.h>
#include <unistd.h>

/*-------------------------------------------------------------------------------------------------------*/
int create_file (char* file, int N, double* position) {

	// save the name of the file, this way we always use the same file throughout one run, initiate counter
	int counter		= 0;

	// identifiers of all relevant spaces and chunking properties
	hid_t	file_id, pos_dataspace_id, pos_dataset_id, tempset_id;
	hid_t	attr_space_id_1, attr_space_id_2, attr_id_1, attr_id_2;
	hid_t	property;

	// checking value whether operations were successful
	herr_t	status;

	// set the size of all dimensions
	hsize_t dims[2]		= {1, 2*N};
	hsize_t offset[2] 	= {0, 0};
	hsize_t maxdims[2] 	= {1, 2*N};

	hsize_t attr_dim	= 1;

	// check whether file already exists, if not create new file using the identifier
	if (access(file, F_OK) != -1) {
		fprintf(stderr, "File already exists, process will terminate.\n");
		return EXIT_FAILURE;
	}

	file_id = H5Fcreate(file, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);

	// create a new dataspace for the data of all particles
	pos_dataspace_id 	= H5Screate_simple(2, dims, maxdims);

	// create a property list enabling chunking
	property		= H5Pcreate(H5P_DATASET_CREATE);
	status			= H5Pset_chunk(property, 2, dims);

	// create a new datasets for all particles, at this point the dataset will have dimension zero
	pos_dataset_id	= H5Dcreate2(file_id, "/positions", H5T_NATIVE_DOUBLE, pos_dataspace_id, H5P_DEFAULT, property, H5P_DEFAULT);

	// select a hyperslab (in this case the entire available space)
	status			= H5Sselect_hyperslab(pos_dataspace_id, H5S_SELECT_SET, offset, NULL, dims, NULL);

	// create a temporary dataspace, which will be used to write the current positions to the file
	tempset_id = H5Screate_simple(2, dims, NULL);

	// write current positions to the file, check whether operation was successful
	status = H5Dwrite(pos_dataset_id, H5T_NATIVE_DOUBLE, tempset_id, pos_dataspace_id, H5P_DEFAULT, position);

	if (status < 0) {
		fprintf(stderr, "Writing of configuration data failed, process will terminate.\n");
		return EXIT_FAILURE;
	}

	// create new dataspaces for the attributes (writeouts, N, amount of data written, and name of the config file) and the attributes themselves
	attr_space_id_1 = H5Screate_simple(1, &attr_dim, NULL);
	attr_space_id_2 = H5Screate_simple(1, &attr_dim, NULL);
	attr_id_1		= H5Acreate2(pos_dataset_id, "N", H5T_NATIVE_INT, attr_space_id_1, H5P_DEFAULT, H5P_DEFAULT);
	attr_id_2		= H5Acreate2(pos_dataset_id, "Last_Writeout", H5T_NATIVE_INT, attr_space_id_2, H5P_DEFAULT, H5P_DEFAULT);

	// write the attributes
	status = H5Awrite(attr_id_1, H5T_NATIVE_INT, &N);
	status = H5Awrite(attr_id_2, H5T_NATIVE_INT, &counter);

	// terminate access to all parts of the document and close the file
	status = H5Aclose(attr_id_1);
	status = H5Aclose(attr_id_2);
	status = H5Sclose(attr_space_id_1);
	status = H5Sclose(attr_space_id_2);

	status = H5Pclose(property);
	status = H5Dclose(pos_dataset_id);	
	status = H5Sclose(pos_dataspace_id);
	status = H5Fclose(file_id);

	// return to caller
	return EXIT_SUCCESS;
}

///*-------------------------------------------------------------------------------------------------------*/
//int write_data(char* file, double* position) {
//	// variables for IDs of the file, dataset and attribute
//	hid_t	file_id, pos_dataset_id, attr_id;
//	hid_t	tempset_id, pos_dataspace_id;//

//	// default error variable
//	herr_t	status;//

//	// offset, set to {0,0}
//	hsize_t offset[2] = {0, 0};//

//	// open an existing file, check whether operation was successful
//	if ((file_id = H5Fopen(file, H5F_ACC_RDWR, H5P_DEFAULT)) < 0) {
//		fprintf(stderr, "Newly created file could not be opened, process will now terminate!\n", timestep);
//		return EXIT_FAILURE;
//	}//

//	// open an existing dataset
//	pos_dataset_id 	= H5Dopen2(file_id, "/positions", H5P_DEFAULT);//

//	// create a temporary dataspace, which will be used to write the current positions to the file
//	tempset_id = H5Screate_simple(2, chunkdims, NULL);//

//	// get the current dataspace and select one hyperslab
//	pos_dataspace_id 	= H5Dget_space(pos_dataset_id);
//	status				= H5Sselect_hyperslab(pos_dataspace_id, H5S_SELECT_SET, offset, NULL, chunkdims, NULL);//

//	// write current positions to the file, check whether operation was successful
//	status = H5Dwrite(pos_dataset_id, H5T_NATIVE_DOUBLE, tempset_id, pos_dataspace_id, H5P_DEFAULT, position);//

//	if (status < 0) {
//		fprintf(stderr, "Writing of configuration data failed on timestep: %d, program will now terminate.\n", timestep);
//		return EXIT_FAILURE;
//	}//

//	// open an attribute, increase the current counter and write the counter to the attribute
//	counter++;
//	attr_id = H5Aopen(pos_dataset_id, "Last_Writeout", H5P_DEFAULT);
//	status 	= H5Awrite(attr_id, H5T_NATIVE_INT, &counter);//
//

//	// terminate access to all parts of the file and close the latter
//	status = H5Sclose(tempset_id);
//	status = H5Aclose(attr_id);
//	status = H5Sclose(pos_dataspace_id);
//	status = H5Dclose(pos_dataset_id);
//	status = H5Fclose(file_id);//

//	// return to caller
//	return EXIT_SUCCESS;
//}

/*-------------------------------------------------------------------------------------------------------*/
int main(int argcount, char** argvector) {

	// variable which we will return and will hold the last position data
	double* position;
	int 	counter, N;

	char 	infile[1024];
	char 	outfile[1024];

	// various identifiers
	hid_t	file_id, pos_dataset_id, attr_id_one, attr_id_two;
	hid_t	pos_dataspace_id, tempspace_id;

	// default error variable
	herr_t 	status;

	// check whether everything needed is provided
	if (argcount != 3) {
		fprintf(stderr, "Please provide a file to read and a filename for the new file, process will terminate.\n");
		return EXIT_FAILURE;
	}
	// get the file to open and the file to write to
	strncpy(infile, argvector[1], 1024);
	strncpy(outfile, argvector[2], 1024);

	// open file, dataset and an attribute, check if successful
	file_id 	= H5Fopen(infile, H5F_ACC_RDWR, H5P_DEFAULT);

	if (file_id < 0) {
		fprintf(stderr, "Initial file could not be opened, process will terminate now \"%s\"\n", infile);
		return EXIT_FAILURE;
	}

	pos_dataset_id 	= H5Dopen2(file_id, "/positions", H5P_DEFAULT);
	attr_id_one		= H5Aopen(pos_dataset_id, "Last_Writeout", H5P_DEFAULT);
	attr_id_two 	= H5Aopen(pos_dataset_id, "N", H5P_DEFAULT);

	// copy the attributes into counter and N variable
	status	= H5Aread(attr_id_one, H5T_NATIVE_INT, &counter);
	status	= H5Aread(attr_id_two, H5T_NATIVE_INT, &N);

	// try to get memory from the system
	position = malloc(2*N*sizeof(double));

	if (position == NULL) {
		fprintf(stderr, "Memory for copying position data could not be allocated. \n");
		return EXIT_FAILURE;
	}

	// variables for offset and size of data copied
	hsize_t offset[2] 	= {counter, 0};
	hsize_t slabdim[2] 	= {1, 2*N};

	// select a hyperslab of the last positions stored in the file and save it to the position array
	tempspace_id 	 = H5Screate_simple(2, slabdim, NULL);
	pos_dataspace_id = H5Dget_space(pos_dataset_id);

	status = H5Sselect_hyperslab(pos_dataspace_id, H5S_SELECT_SET, offset, NULL, slabdim, NULL);
	status = H5Dread(pos_dataset_id, H5T_NATIVE_DOUBLE, tempspace_id, pos_dataspace_id, H5P_DEFAULT, position);

	// terminate all access to and close the file
	status = H5Aclose(attr_id_one);
	status = H5Aclose(attr_id_two);
	status = H5Dclose(pos_dataset_id);
	status = H5Sclose(tempspace_id);
	status = H5Sclose(pos_dataspace_id);
	status = H5Fclose(file_id);

	if (create_file(outfile, N, position) != EXIT_SUCCESS) {
		return EXIT_FAILURE;
	}

	// return to caller
	return EXIT_SUCCESS;
}