# Source code for Selective HubSorting. 

The implementation first computes the Packing Factor and conditionally
hubsorts the graph (_if_ the packing factor of the graph is high)

## Instructions

1. The `apps` directory contains application code and the `Makefile`. (_Directory contains
PageRank as an example application_)

2. Go to the `apps` directory and simply type `make`. **NOTE:** The reordering routines were
developed and tested using OpenMP. Running with Cilk might cause issues.

3. To enable selective reordering, launch `PageRank` with the -preprocess argument. Eg. 
`./PageRank -preprocess 0 <path-to-ligra-formatted-graph>`
