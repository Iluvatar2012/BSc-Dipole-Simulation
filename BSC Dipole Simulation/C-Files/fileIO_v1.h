/*
 * fileIO_v1.h
 *
 *  	Last Changed: 	Aug 1, 2013
 *      Author: 		Aiko Bernehed
 */

#ifndef FILEIO_H_
#define FILEIO_H_

extern int 	write_file 	(char *, double *, int, int, int);
extern int 	read_file	(char*, char*, int*, double*, double*, double*, \
						 double*, double*, double*, int*, int*,	int*);

#endif /* FILEIO_H_ */
