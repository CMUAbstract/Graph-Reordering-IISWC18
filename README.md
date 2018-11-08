# Lightweight Graph Reordering (IISWC18)

Repo for the IISWC18 paper - "When is Graph Reordering an Optimization? 
Studying the Effect of Lightweight Graph Reordering Across Applications and Input Graphs" [pre-print](https://users.ece.cmu.edu/~vigneshb/papers/IISWC2018-final-preprint.pdf)

## Contents

The repo contains source code: 

* Packing Factor computation
* Lightweight Reordering techniques (incorporated in the Ligra applications):
    * Hub Sorting (based on [^fn1])
    * Hub Clustering

Details about finding the packing factor of graphs (must be in ligra-format [^fn2]) 
and applying lightweight reordering can be found in the `Ligra` directory (**NOTE:**
_Extensions to the GAP benchmarks coming soon_)

## Requirements

Ligra and GAP require g++ version >= 5.3.0

Both benchmark suites are parallelized with OpenMP

Tested on Debian Stretch with g++ version 6.3.0

## Issues

For issues/information, please feel free to send an email at `vigneshb@andrew.cmu.edu` 


## References

[^fn1]: Zhang, Yunming, et al. "Making caches work for graph analytics." Big Data (Big Data), 
2017 IEEE International Conference on. IEEE, 2017

[^fn2]: [format](https://github.com/jshun/ligra#input-format-for-ligra-applications-and-the-ligra-encoder)
