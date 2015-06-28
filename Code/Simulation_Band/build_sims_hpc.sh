#!/bin/bash

for X in 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
	./configure $X
	scp simulation.c aiber100@hpc:~/Simulation_Band/

	ssh aiber100@hpc "module load intel/xe2013; source .bashrc; make; mkdir Work/${X}; cp Work/simulation_band Work/${X}/; exit"
done