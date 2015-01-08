#!/bin/bash

# Define paths and filenames
EXT=../Rohdaten/HPC_100_tau_shear_500
FILES=$EXT/*.hdf5
DIR=$EXT/final_height_diff

# make a new directory for the final result
mkdir $DIR

# remove any previous analysis file
rm final_stat*

# iterate over all present HDF5 Files
for F in $FILES
do
	# extract the value of m and gamma
	M=$(echo $F | cut -d '_' -f 6)
	G=$(echo $F | cut -d '_' -f 9)

	# OUTFILE=final_stat-$M-$G

	# read the file and write to plaintext, execute ruby script
	./write_to_file $F temp 0.001 100000
	./final_step.rb temp $M $G

	# remove temporary files
	rm temp
	# rm final_stat
done

# plot the data 
./plotter.sh
