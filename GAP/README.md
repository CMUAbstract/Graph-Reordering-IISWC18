# Lightweight Graph Reordering techniques for GAP

The directories contain implementation of the following LWRs:
1. Degree Sorting
2. Hub Sorting
3. Hub Clustering

The reordering function is in `*/builder.h`

Each directory contains `preprocessor-<appname>.cc` files that measure
the minimum preprocessing cost incurred for each application. For example,
GAP's implementation of pull-based PageRank does not need to build the 
complete outCSR (the application uses only out-degrees and not outgoing-neighbors)
