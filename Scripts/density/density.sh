#!/bin/bash

# Define paths and filenames
EXT=../Rohdaten/HPC_100_tau_shear_200
FILES=$EXT/*.hdf5
DIR=$EXT/density

# make the directory if it doesn't exist 
if [ ! -d "$DIR" ]; then
	mkdir $DIR
fi

# iterate over all files in the folder
for F in $FILES
do
	# extract the value of m and gamma
	M=$(echo $F | cut -d '_' -f 6)
	G=$(echo $F | cut -d '_' -f 9)

	# compute the file name were the final data will be stored
	OUT=$DIR/
	OUT+=m_$M
	OUT+=_Gamma_$G

	# read the file and write to plaintext, execute ruby script
	./write_to_file $F temp 0.001 100
	./density.rb temp 1000 0.1 density_stat

	# call the plotter
	./plotter.sh $OUT density_stat

	# remove the evidence
	rm temp
	rm density_stat
done
