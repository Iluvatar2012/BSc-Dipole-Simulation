#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#define PI 		3.14159265358979323846264338328

extern void 	init (struct variables*);
extern void		compute_psi4 (int);
extern void		compute_psi6 (int);
extern void 	bubble_sort (double*, int*, int);
extern double 	psi_n (int, int, int, int*);

#endif