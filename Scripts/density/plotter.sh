#!/bin/bash

# pipe everything into gnuplot, this way command line arguments can still be given, since this is a bash script
gnuplot <<- EOF

	# set terminal and output, use input filename as base for new filename
	set terminal pngcairo size 1920,1080 enhanced font 'Verdana,30'
	set output '$1.png'

	# set view and pm3d for the color palette
	set view map
	set pm3d

	# set x and yrange and the range of the color palette
	set xr [0:100]
	set yr [0:21.466]

	set cbrange [0:1]

	# set the used borders and delete the mirrored tics, this makes for a much nicer view
	set border 10

	set xtics nomirror
	set ytics nomirror

	# plot the data 
	splot "$2" u 1:2:5 w pm3d title ""

EOF