#!/bin/bash

for i in {151..180..5}
do
	echo "Processing $i"
	cp /scratch_gs/aiber100/51576\[$i\].hpc-batch/Results/* ~/Results/Test_2/ &
	cp /scratch_gs/aiber100/51629\[$i\].hpc-batch/Results/* ~/Results/Test_3/ &
done

echo "Processing 180"

cp /scratch_gs/aiber100/51577\[180\].hpc-batch/Results/* ~/Results/Test_2/ &
cp /scratch_gs/aiber100/51630\[180\].hpc-batch/Results/* ~/Results/Test_3/

echo "done"
