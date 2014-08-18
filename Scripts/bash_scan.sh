#!/bin/bash

# get the current iteration number
No=$1

# iterate what G is given, we assume 5 shearing values and thus 30 simulations before a different m is reached
G=$(($No-1))
G=$(($G / 5))
G=$(($G * 5))
G=$(($G + 1))
G=$(($G % 30))

# check which m we need
e=$(($No - 1))
e=$(($e / 30))

# set the according m value
if [ "$e" = "0" ]; then 
	m=0.50
fi
if [ "$e" = "1" ]; then 
	m=0.60
fi
if [ "$e" = "2" ]; then 
	m=0.70
fi
if [ "$e" = "3" ]; then 
	m=0.80
fi
if [ "$e" = "4" ]; then 
	m=0.90
fi
if [ "$e" = "5" ]; then 
	m=1.00
fi

mid=__Gamma_
end=.0.hdf5

# start the simulation
echo "./simulation Config_Files/$No $No Equalised/Eqq_m_$m$mid$G$end"
