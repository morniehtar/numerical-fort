#! /usr/bin/gnuplot

set encoding iso_8859_1
set tmargin 1
set tics out
set nokey
set grid

plot "./output.dat" using 1:2 smooth csplines, \
"./output.dat" using 1:3 smooth csplines, \
"./output.dat" using 1:4 smooth csplines

pause mouse close

