Lower bound calculation of edges using ILP
================

This program calculates a lowerbound for 4-factor area by area.
On the edge of the area, the edges will be "directional", in the inside, edges must
be undirected. The program also outputs the edges. This is useful for comparing your solution
against optimal almost-4-factor and spotting "alternating cycles" between your current solution
and the lowerbound.

1. make
2. generate ../adjust/adjust.dat
3. cp ../adjust/adjust.dat
4. ./lowerbound\_edges.cpp.bin
5. use edge.dat for whatever you want

You may play with SOLVE\_SIZE constant (we were being able to go up to 4000)
and MAX\_RUNNING\_TASKS (for big tasks use >=2 only if you have plenty of ram).

Note: edge.dat is created incrementally.
