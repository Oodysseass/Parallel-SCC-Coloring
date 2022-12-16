# Parallel-SCC-Coloring

Algorithm for detection of strongly connected components of a graph.
Graphs are given to the program as a command line argument with an .mtx file in matrix market format, which contains the graph information in COO.

There are three implementations of parallelism: one with openCilk, one with openMp and one with PThreads.

Clang compiler of the OpenCilk is required for the openCilk version which can be found here: https://www.opencilk.org/doc/users-guide/install/
Other implementations can be used with an up-to-date gcc compiler with -fopenmp and -pthread flags accordingly.

Make sure to have an .mtx file from SuiteSparse Matrix Collection, https://suitesparse-collection-website.herokuapp.com/.

To run: clone this repo, run make and use the correspoding binary with an .mtx file as argument.

e.g.

      ./ColoringSCCopenCilk celegansneural.mtx

runs the openCilk implementation with celegansneural's graph as an argument.

Accordingly,

      ./ColoringSCCSequential wiki-topcats.mtx
     
runs the sequential implementation with wiki-topcats' graph as an argument.
