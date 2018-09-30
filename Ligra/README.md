# Lightweight Reordering for the Ligra Benchmark suite

## List of contents

* `ligra-hubsort`     - Unconditional Frequency-based Clustering (AKA. "HubSorting") applied to ligra applications
* `ligra-hubcluster`  - Unconditional HubClustering applied to ligra applications
* `selective-hubsort` - Selective HubSorting based on Packing Factor of the input graph

## Navigation 

Code for reordering and packing factor computation is in `ligra/IO.h` in each directory. 

Reordering is implemented in a function called `preprocessGraph`

Packing Factor computation is implemented in a function called `computePackingFactor`


