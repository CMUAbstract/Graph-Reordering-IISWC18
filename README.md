# Lightweight Graph Reordering (IISWC18)

Repo for the IISWC18 paper - "When is Graph Reordering an Optimization? 
Studying the Effect of Lightweight Graph Reordering Across Applications and Input Graphs" [pre-print](https://users.ece.cmu.edu/~vigneshb/papers/IISWC2018-final-preprint.pdf)

## Contents

The repo contains source code: 

* Packing Factor computation
* Lightweight Reordering techniques (incorporated in the Ligra applications):
    * Hub Sorting (based on [1])
    * Hub Clustering

Details about finding the packing factor of graphs (must be in [ligra-format](https://github.com/jshun/ligra#input-format-for-ligra-applications-and-the-ligra-encoder)) 
and applying lightweight reordering can be found in the `Ligra` and `GAP` directories

## Tools

For running algorithms beyond those in the GAP and Ligra benchmark suites, 
there is a tool to generate a reordered graph (as an edgelist or in CSR format)
which can be used as an input for the new algorithm.

In the `GAP/*` directories, `graph-reorderer` will output reordered graphs. 
The `-l` input flag allows choosing between out-degree vs in-degree based reordering.
The `-e` flag outputs the graph in the edgelist format (note that the edgelist
is sorted by source ids). 
The `-b` flag outputs the graph in the serialized graph format used in GAP

## Requirements

Ligra and GAP require g++ version >= 5.3.0

Both benchmark suites are parallelized with OpenMP

Tested on Debian Stretch with g++ version 6.3.0

## Issues

For issues/information, please feel free to send an email at `vigneshb@alumni.cmu.edu` 


## References

[1]: Zhang, Yunming, et al. "Making caches work for graph analytics." Big Data (Big Data), 
2017 IEEE International Conference on. IEEE, 2017

