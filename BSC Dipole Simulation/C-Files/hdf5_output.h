/*
 * hdf5_write.h
 *
 *  Created on: Oct 23, 2013
 *      Author: Aiko Bernehed
 */

#ifndef HDF5_OUTPUT_H_
#define HDF5_OUTPUT_H_

extern int 	create_file (char*, int, int);
extern void	write_data(int, double*);

#endif /* HDF5_WRITE_H_ */
