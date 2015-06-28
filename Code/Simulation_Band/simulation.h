/*
 * Simulation.h
 *
 *  	Last Changed: 	Aug 1, 2013
 *      Author: 		Aiko Bernehed
 */

#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "structs.h"

extern int 	init(struct parameters*, double*);
extern void simulation(void);
extern void	update_verlet(void);

#endif /* SIMULATION_H_ */
