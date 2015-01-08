#!/bin/bash

export SSHPASS='Bumblebee_Rules'

for i in {1..180}
do
	scp aiber100@hpc:/scratch_gs/aiber100/75782\\\[$i\\\].hpc-batch/Results/* ~/Aiko/Temp/
	sshpass -e scp ~/Aiko/Temp/* admin@iluvatar.servebeer.com:/share/Bumblebee/Bachelor/Rohdaten/HPC_100_tau_shear_200/

	rm ~/Aiko/Temp/*
done

export SSHPASS=''

# 
# for i in {151..180..5}
# do
# 	echo "Processing $i"
# 	cp /scratch_gs/aiber100/51576\[$i\].hpc-batch/Results/* ~/Results/Test_2/ &
# 	cp /scratch_gs/aiber100/51629\[$i\].hpc-batch/Results/* ~/Results/Test_3/ &
# done
# 
# echo "Processing 180"
# 
# cp /scratch_gs/aiber100/51577\[180\].hpc-batch/Results/* ~/Results/Test_2/ &
# cp /scratch_gs/aiber100/51630\[180\].hpc-batch/Results/* ~/Results/Test_3/
# 
# echo "done"
