/*
 * Simulation.h
 *
 *  	Last Changed: 	Aug 1, 2013
 *      Author: 		Aiko Bernehed
 */

#ifndef SIMULATION_H_
#define SIMULATION_H_

#define PI 			3.14159265358979323846264338328
#define TWO_PI		6.28318530717958647692528676656

extern int 	init(void);
extern void simulation(void);
extern void	update_verlet(void);

#endif /* SIMULATION_H_ */
