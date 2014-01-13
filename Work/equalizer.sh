#!/bin/bash

start="Start/Eqq_m_0.1__G_"
end=".hdf5"

./simulation Config_Files/5 5

for ((i=5; i<=45; i=i+5))

do
	((j=i+5))
	arg=$start$i$end

	./simulation Config_Files/$j $j $arg

done