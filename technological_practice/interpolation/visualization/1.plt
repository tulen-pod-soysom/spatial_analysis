#!/bin/gnuplot
set terminal png size 1024, 1024
set palette gray
set hidden3d
set pm3d at bs


set output "im.png"
# plot "im.txt" matrix with image
splot "im.txt" matrix with lines
