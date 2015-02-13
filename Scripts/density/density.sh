#!/bin/bash

# Define paths and filenames
EXT=.
FILES=$EXT/*.hdf5
DIR=$EXT

# make the data directory if it doesn't exist 
if [ ! -d "$DIR" ]; then
	mkdir $DIR
fi

# make a temporary directory in the home folder if it doesn't exist 
if [ ! -d "/home/ayra/.temp" ]; then
	mkdir $/home/ayra/.temp
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
	./write_to_file $F /home/ayra/.temp/temp 0.001 100
	./density.rb /home/ayra/.temp/temp 1000 0.1 /home/ayra/.temp/density_stat

	# call the plotter
	./plotter.sh $OUT /home/ayra/.temp/density_stat

	# remove the evidence
	rm /home/ayra/.temp/temp
	rm /home/ayra/.temp/density_stat
done
