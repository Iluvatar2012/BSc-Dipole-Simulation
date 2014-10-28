#!/bin/bash

EXT=.
FILES=$EXT/*.hdf5
DIR=$EXT/mean_height_pics

# build a new directory to hold the pictures
mkdir $DIR

# iterate over all files in the given directory
for F in $FILES
do 
    # write data to a temporary file and do a short analysis
    ./write_to_file $F temp 0.001 100
    ./average_part.rb temp

    # compute the basename of the input file
    TEMP=${F##*/}
    TEMP=${TEMP%.*}

    # plot the data
    ./plotter.sh $DIR/$TEMP particle_stat

    # remove temporary files
    rm temp
    rm particle_stat
done