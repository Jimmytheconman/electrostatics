#!/usr/bin/gnuplot -persist
load '../scripts/colormap.pal'
set size ratio -1
set term postscript color
set pm3d map
set xlabel "i"
set ylabel "j"
set nokey

xMin = -250
yMin = -250
set xrange [-250 : 250];
set yrange [-250 : 250];

set output "numericalproblem2.eps"
set title "Numerical solution for problem 2"
splot "../bin/numericalproblem2" using ($1+xMin):($2+yMin):3 matrix

set output "analyticalproblem2.eps"
set title "Analytical solution for problem 2"
splot "../bin/analyticalproblem2" using ($1+xMin):($2+yMin):3 matrix

set output "differenceproblem2.eps"
set title "Difference between solutions for problem 2"
splot "../bin/differenceproblem2" using ($1+xMin):($2+yMin):3 matrix

set output "numericalproblem2contour.eps"
set title "Numerical equipotentials for problem 2"
set isosamples 21, 21
set contour surface
set cntrparam levels auto 10
set style data lines
splot "../bin/numericalproblem2" using ($1+xMin):($2+yMin):3 matrix with lines nosurf lw 3

set output "numericalfieldproblem2.eps"
set title "Numerical field for problem 2"
set parametric;
set isosamples 100,100
plot "../bin/numericalfieldproblem2" every 20:20 using 1:2:($3*40):($4*40) with vectors
