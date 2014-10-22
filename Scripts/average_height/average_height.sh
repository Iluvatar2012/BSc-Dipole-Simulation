#!/bin/bash

FILES=./*.hdf5
DIR=mean_height_pics

# build a new directory to hold the pictures
mkdir $DIR

# iterate over all files in the given directory
for F in $FILES
do 
	# write data to a temporary file and do a short analysis
	./write_to_file $F temp 1e-5 1e4
	./average_part_a.rb temp

	# compute the basename of the input file
	TEMP=${F##*/}
	TEMP=${TEMP%.*}

	# plot the data
	./plotter.sh $DIR/$TEMP particle_a_stat

	# remove temporary files
	rm temp
	rm particle_a_stat
done