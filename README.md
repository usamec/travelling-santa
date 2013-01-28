travelling-santa
================

## How to get a solution from this mess:

1. Produce starting paths using stuff in starting\_path directory
2. Use simple (but fast) improvements in simple\_opt directory to get lenght of paths to something
around 7 milion
3. Use mega\_opt to reduce length even further. 
4. Use stuff in final\_push to get length of path even more down. But this stuff optimizes average
lenght of paths. So use mega\_opt to bring paths lenghts back together.

There is also our implementation of ILP + SAT + custom spliting pipeline. It's in directory ilp. You
don't need to use it.

##How to get a good lowerbound:

1. Produce ajdust.dat using lower\_bound/adjust (you also obtain a reasonably good lowerbound)
2. Get lowerbound for cycles and paths using lower\_bound/ilp\_lowerbound
3. Get lowerbound edges using lower\_bound/ilp\_edges
