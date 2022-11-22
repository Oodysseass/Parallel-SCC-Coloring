# Parallel-SCC-Coloring

Algorithm for detection of strongly connected components of a graph.
Graphs are given to the program as a command line argument with an .mtx file in matrix market format, which contains the graph information in COO.

There are three implementations of parallelism: one with openCilk, one with openMp and one with PThreads.
