#!/bin/bash

# pass following commands to gnuplot
gnuplot <<- EOF

# set terminal and output, use input filename as base for new filename
set terminal pngcairo size 1920,1080 enhanced font 'Verdana,30'
set output '$1.png'

# set grid
set grid

# set linestyle (a little fat in blue) and fillstyle
set style line 1 lc 3 lw 4
set style line 2 lc rgb "forest-green"

set style fill transparent solid 0.4 noborder

# set labels

# set x and yrange
set xr [0:100]
set yr [0:25]

# set labels and lines if necessary

# plot and return to caller, second argument is the actual data
plot '$2' u 1:2:(\$2-\$3) w filledcurves lines 2 title "", \
	 '$2' u 1:2:(\$2+\$3) w filledcurves lines 2 title "", \
	 '$2' u 1:2 w lines lines 1 title ""

EOF