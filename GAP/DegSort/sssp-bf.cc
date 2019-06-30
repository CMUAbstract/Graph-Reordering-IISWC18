// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <cinttypes>
#include <limits>
#include <iostream>
#include <queue>
#include <vector>
#include <atomic>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "platform_atomics.h"
#include "pvector.h"
#include "timer.h"
#include "bitmap.h"


/*
GAP Benchmark Suite
Kernel: Single-source Shortest Paths (SSSP)
Author: Vignesh Balaji

Returns array of distances for all vertices from given source vertex

Uses the bellman-ford algorithm to compute distances from a source. 
The implementation applies Yen's optimizations and is based on the 
Ligra implementation of Bellman Ford. Note, that this is a pull-based
implementation 
*/


using namespace std;

const WeightT kDistInf = numeric_limits<WeightT>::max()/2;
const size_t kMaxBin = numeric_limits<size_t>::max()/2;

pvector<WeightT> BellmanFord(const WGraph &g, NodeID source) {
    pvector<WeightT> dist(g.num_nodes(), kDistInf);
    dist[source] = 0;
    pvector<AtomicWeightT> actual_dist(g.num_nodes());
    //initializing atomic type
    #pragma omp parallel for 
    for (NodeID n = 0; n < g.num_nodes(); ++n) {
        //actual_dist[n] = static_cast<AtomicWeightT>(dist[n]);
        actual_dist[n].store(dist[n]);
    }
    NodeID iter(0);
    while (iter < g.num_nodes()) {
        Timer t;
        t.Start();
        NodeID numReductions(0);
        #pragma omp parallel for schedule(dynamic, 64) reduction(+ : numReductions)
        for (NodeID n = 0; n < g.num_nodes(); ++n) {
            for (auto wm : g.in_neigh(n)) {
                AtomicWeightT new_dist(actual_dist[wm.v].load() + wm.w);
                if (actual_dist[n] > new_dist) {
                    actual_dist[n].store(new_dist);
                    ++numReductions;
                }
            }
        }
        if (numReductions == 0) {
            break;
        }
        ++iter;
        t.Stop();
        PrintTime("Iter-time", t.Seconds());
    }
    //reverting back to simple type
    #pragma omp parallel for 
    for (NodeID n = 0; n < g.num_nodes(); ++n) {
        dist[n] = actual_dist[n].load();
    }
    std::cout << "Bellman Ford took " << iter << " iterations to converge " << std::endl;
    return dist;
}


void PrintSSSPStats(const WGraph &g, const pvector<WeightT> &dist) {
  auto NotInf = [](WeightT d) { return d != kDistInf; };
  int64_t num_reached = count_if(dist.begin(), dist.end(), NotInf);
  cout << "SSSP Tree reaches " << num_reached << " nodes" << endl;
}


// Compares against simple serial implementation
bool SSSPVerifier(const WGraph &g, NodeID source,
                  const pvector<WeightT> &dist_to_test) {
  // Serial Dijkstra implementation to get oracle distances
  pvector<WeightT> oracle_dist(g.num_nodes(), kDistInf);
  oracle_dist[source] = 0;
  typedef pair<WeightT, NodeID> WN;
  priority_queue<WN, vector<WN>, greater<WN>> mq;
  mq.push(make_pair(0, source));
  while (!mq.empty()) {
    WeightT td = mq.top().first;
    NodeID u = mq.top().second;
    mq.pop();
    if (td == oracle_dist[u]) {
      for (WNode wn : g.out_neigh(u)) {
        if (td + wn.w < oracle_dist[wn.v]) {
          oracle_dist[wn.v] = td + wn.w;
          mq.push(make_pair(td + wn.w, wn.v));
        }
      }
    }
  }
  // Report any mismatches
  bool all_ok = true;
  for (NodeID n : g.vertices()) {
    if (dist_to_test[n] != oracle_dist[n]) {
      cout << n << ": " << dist_to_test[n] << " != " << oracle_dist[n] << endl;
      all_ok = false;
    }
  }
  return all_ok;
}


int main(int argc, char* argv[]) {
  CLDelta<WeightT> cli(argc, argv, "single-source shortest-path");
  if (!cli.ParseArgs())
    return -1;
  WeightedBuilder b(cli);
  WGraph g = b.MakeGraph();
  pvector<NodeID> new_ids (g.num_nodes(), -1);
  WGraph sorted_g = WeightedBuilder::degreeSort_weighted(g, (cli.relabel() == 0), new_ids, false, false);
  SourcePicker<WGraph> sp(g, cli.start_vertex());
  auto SSSPBound = [&sp, &cli, &new_ids] (const WGraph &sorted_g) {
    return BellmanFord(sorted_g, new_ids[sp.PickNext()]);
  };
  SourcePicker<WGraph> vsp(g, cli.start_vertex());
  auto VerifierBound = [&vsp, &new_ids] (const WGraph &sorted_g, const pvector<WeightT> &dist) {
    return SSSPVerifier(sorted_g, new_ids[vsp.PickNext()], dist);
  };
  BenchmarkKernel(cli, sorted_g, SSSPBound, PrintSSSPStats, VerifierBound);
  return 0;
}
