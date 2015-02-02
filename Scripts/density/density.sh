#!/bin/bash

# Define paths and filenames
EXT=.
FILES=$EXT/*.hdf5
DIR=$EXT

for F in $FILES
do
	# read the file and write to plaintext, execute ruby script
	./write_to_file $F temp 0.001 1000
	./density.rb temp 100 1

	# remove the evidence
	rm temp
done