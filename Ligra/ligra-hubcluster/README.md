# Source code for HubClustering 

HubClustering is a variant of frequency based clustering that 
does not require sorting (only keeps the hubs together)

## Instructions

1. The `apps` directory contains application code and the `Makefile`. 

2. Go to the `apps` directory and simply type `make`. **NOTE:** The reordering routines were
developed and tested using OpenMP. Running with Cilk might cause issues.

3. To unconditionally perform hubclustering, launch application with the `-preprocess <0/1>` argument. Eg. 
`./PageRank -preprocess 0 <path-to-ligra-formatted-graph>` (_for out-degree hub-sorting_)

4. To perform timing analysis on hubsorting for different types of graphs (directed, undirected, weighted, bipartite), 
run the `preprocessor-*` applications
