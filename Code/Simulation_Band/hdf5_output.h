#ifndef HDF5_OUTPUT_H_
#define HDF5_OUTPUT_H_

extern int 		create_file (char*, struct attributes*);
extern void		write_data(int, double*, double*);
extern int 		reopen_file(double*, double*, int, int*);
extern double 	*read_configuration(char*);

#endif /* HDF5_WRITE_H_ */
