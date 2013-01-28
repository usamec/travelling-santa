Lower bound calculation using ILP
================

This program calculates a lowerbound using 4-factor for areas of points.
On the edge of the area, the edges will be "directional", in the inside, edges must
be undirected.

1. make
2. generate ../adjust/adjust.dat
3. cp ../adjust/adjust.dat
4. ./ilp\_lowerbound.cpp.bin

You may play with SOLVE\_SIZE constant. 
