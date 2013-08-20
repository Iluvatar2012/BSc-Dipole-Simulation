/*
 * fileIO_v1.h
 *
 *  	Last Changed: 	Aug 1, 2013
 *      Author: 		Aiko Bernehed
 */

#ifndef FILEIO_H_
#define FILEIO_H_

extern int 					write_file 	(char *, double *, int, int, int);
extern struct parameters*	read_file	(char*);

#endif /* FILEIO_H_ */
