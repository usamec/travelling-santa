Lower bound calculation
================

This program calculates a 4-nearest-neighbours lowerbound.
Because of lowerbound can be greatly improved if calculated on a different cost matrix,
the program uses gradient ascent to obtain vector of values P\_i (new distance matrix D\_i,j = C\_i,j + P\_i + P\_j
and also calculate an estimate for both circles/paths lower bounds.

1. make
2. ./adjust.cpp.bin

You may play with the number of iterations/step size as current numbers are rather generic
(I was manually updating stepsize all the time but I wanted to make it easy for others to run)
