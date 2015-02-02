#!/bin/bash

for i in {150..153}
do
	scp aiber100@hpc:/scratch_gs/aiber100/75782\\\[$i\\\].hpc-batch/Results/* ~/Downloads/Data/
	
done
