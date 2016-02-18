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
set yrange [-250 : 250];
set xrange [-250 : 250];

set output "numericalproblem1.eps"
set title "Numerical solution for problem 1"
splot "../bin/numericalproblem1" using ($1+xMin):($2+yMin):3 matrix

set output "analyticalproblem1.eps"
set title "Analytical Soultion for problem 1"
splot "../bin/analyticalproblem1" using ($1+xMin):($2+yMin):3 matrix

set output "differenceproblem1.eps"
set title "Difference between analytical and numerical solutions for problem 1"
splot "../bin/differenceproblem1" using ($1+xMin):($2+yMin):3 matrix

set output "numericalproblem1contour.eps"
set title "Numerical equipotentials for problem 1"
set isosamples 21, 21
set contour surface
set cntrparam levels auto 10
set style data lines
splot "../bin/numericalproblem1" using ($1+xMin):($2+yMin):3 matrix with lines nosurf lw 3

set output "numericalfieldproblem1.eps"
set title "Numerical field for problem 1"
set parametric;
set isosamples 100,100
plot "../bin/numericalfieldproblem1" every 20:20 using 1:2:($3*40):($4*40) with vectors


