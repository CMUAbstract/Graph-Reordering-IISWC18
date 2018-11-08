# Source code for frequency-based clustering (AKA "HubSorting") 

Implementation as described in Cagra[1].

## Instructions

1. The `apps` directory contains application code and the `Makefile`. 

2. Go to the `apps` directory and simply type `make`. **NOTE:** The reordering routines were
developed and tested using OpenMP. Running with Cilk might cause issues.

3. To unconditionally perform hubsorting, launch application with the `-preprocess <0/1>` argument. Eg. 
`./PageRank -preprocess 0 <path-to-ligra-formatted-graph>` (_for out-degree hub-sorting_)

4. To perform timing analysis on hubsorting for different types of graphs (directed, undirected, weighted, bipartite), 
run the `preprocessor-*` applications


## References 

[1] Zhang, Yunming, et al. "Making caches work for graph analytics." 
Big Data (Big Data), 2017 IEEE International Conference on. IEEE, 2017
