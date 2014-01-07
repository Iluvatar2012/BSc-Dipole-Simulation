/*
 * Simulation.h
 *
 *  	Last Changed: 	Aug 1, 2013
 *      Author: 		Aiko Bernehed
 */

#ifndef SIMULATION_H_
#define SIMULATION_H_

extern struct sim_struct 	*make_sim_struct (void);
extern int 					read_struct (char*);

extern int 	init(struct sim_struct*);
extern void simulation(void);
extern void	update_verlet(void);

#endif /* SIMULATION_H_ */
