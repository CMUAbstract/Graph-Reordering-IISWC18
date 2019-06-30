// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <iostream>
#include <vector>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"


/*
  Kernel to measure the runtime of preprocessing. We have a separate file for this
  because the applications only do graph preprocessing once for multiple runs of 
  the kernel. This file allows capturing the effects of cold-caches in preprocessing.
*/


using namespace std;

const int numRuns {16};

int main(int argc, char* argv[]) {
  CLIterApp cli(argc, argv, "betweenness-centrality", 1);
  if (!cli.ParseArgs())
    return -1;
  {
    Builder b(cli);
    Graph g1 = b.MakeGraph();
    for (int i = 0; i < numRuns ; ++i) 
    {
      {
      pvector<NodeID> newIds(g1.num_nodes(), -1);
      std::cout << "[LOG] normal stat begins" << std::endl;
      Graph sorted_g = Builder::degreeSort(g1, (cli.relabel() == 0), newIds, false, false); 
      std::cout << "[LOG] normal stat ends" << std::endl;
      }
      std::cout << "***\n";
    }
  }
  return 0;
}
