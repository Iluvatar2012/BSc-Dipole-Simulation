#!/bin/bash

# pass following commands to gnuplot
gnuplot <<- EOF

# set terminal and output, use input filename as base for new filename
set terminal pngcairo size 1920,1080 enhanced font 'Verdana,30'
set output '$1.png'

# set grid
set grid

# set linestyle (a little fat in blue) and fillstyle
set style line 1 lc rgb '#08519c' lw 4
set style line 2 lc rgb '#3182bd'
set style line 3 lc rgb '#a50f15' lw 4
set style line 4 lc rgb '#fb6a4a'

set style fill transparent solid 0.4 noborder

# set x and yrange
set xr [0:100]
set yr [0:25]

# plot and return to caller, second argument is the actual data
plot '$2' u 1:2:(\$2-\$3) w filledcurves lines 2 title "", \
	 '$2' u 1:2:(\$2+\$3) w filledcurves lines 2 title "", \
	 '$2' u 1:4:(\$4-\$5) w filledcurves lines 4 title "", \
	 '$2' u 1:4:(\$4+\$5) w filledcurves lines 4 title "", \
	 '$2' u 1:2 w lines lines 1 title "", \
	 '$2' u 1:4 w lines lines 3 title ""

EOF