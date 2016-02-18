#!/usr/bin/gnuplot -persist
load '../scripts/colormap.pal'
set size ratio -1
set term postscript color
set pm3d map
set xlabel "i"
set ylabel "j"
set nokey

xMin = -300
yMin = -100
set xrange [-300 : 300];
set yrange [-100 : 100];

set output "numericalproblem3.eps"
set title "Numerical solution for problem 3"
splot "../bin/numericalproblem3" using ($1+xMin):($2+yMin):3 matrix

set output "numericalproblem3contour.eps"
set title "Numerical equipotentials for problem 3"
set isosamples 21, 21
set contour surface
set cntrparam levels auto 10
set style data lines
splot "../bin/numericalproblem3" using ($1+xMin):($2+yMin):3 matrix with lines nosurf lw 3

set output "numericalfieldproblem3.eps"
set title "Numerical field for problem 3"
set parametric;
set isosamples 100,100
plot "../bin/numericalfieldproblem3" every 20:20 using 1:2:($3*40):($4*40) with vectors
