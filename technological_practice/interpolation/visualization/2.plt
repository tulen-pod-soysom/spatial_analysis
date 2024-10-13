#!/bin/gnuplot
set terminal png size 1024, 1024
set palette gray
set hidden3d
set pm3d at bs


set output "res.png"
# plot "res.txt" matrix with image
splot "res.txt" matrix with lines
