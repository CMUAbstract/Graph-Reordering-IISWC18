// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef BUILDER_H_
#define BUILDER_H_

#include <algorithm>
#include <parallel/algorithm>
#include <cinttypes>
#include <fstream>
#include <functional>
#include <type_traits>
#include <utility>
#include <omp.h>
#include <cassert>
#include <vector>

#include "command_line.h"
#include "generator.h"
#include "graph.h"
#include "platform_atomics.h"
#include "pvector.h"
#include "reader.h"
#include "timer.h"
#include "util.h"
#include "sliding_queue.h"


/*
GAP Benchmark Suite
Class:  BuilderBase
Author: Scott Beamer

Given arguements from the command line (cli), returns a built graph
 - MakeGraph() will parse cli and obtain edgelist and call
   MakeGraphFromEL(edgelist) to perform actual graph construction
 - edgelist can be from file (reader) or synthetically generated (generator)
 - Common case: BuilderBase typedef'd (w/ params) to be Builder (benchmark.h)
*/


template <typename NodeID_, typename DestID_ = NodeID_,
          typename WeightT_ = NodeID_, bool invert = true>
class BuilderBase {
  typedef EdgePair<NodeID_, DestID_> Edge;
  typedef pvector<Edge> EdgeList;

  const CLBase &cli_;
  bool symmetrize_;
  bool needs_weights_;
  int64_t num_nodes_ = -1;

 public:
  explicit BuilderBase(const CLBase &cli) : cli_(cli) {
    symmetrize_ = cli_.symmetrize();
    needs_weights_ = !std::is_same<NodeID_, DestID_>::value;
  }

  DestID_ GetSource(EdgePair<NodeID_, NodeID_> e) {
    return e.u;
  }

  DestID_ GetSource(EdgePair<NodeID_, NodeWeight<NodeID_, WeightT_>> e) {
    return NodeWeight<NodeID_, WeightT_>(e.u, e.v.w);
  }

  NodeID_ FindMaxNodeID(const EdgeList &el) {
    NodeID_ max_seen = 0;
    #pragma omp parallel for reduction(max : max_seen)
    for (auto it = el.begin(); it < el.end(); it++) {
      Edge e = *it;
      max_seen = std::max(max_seen, e.u);
      max_seen = std::max(max_seen, (NodeID_) e.v);
    }
    return max_seen;
  }

  pvector<NodeID_> CountDegrees(const EdgeList &el, bool transpose) {
    pvector<NodeID_> degrees(num_nodes_, 0);
    #pragma omp parallel for
    for (auto it = el.begin(); it < el.end(); it++) {
      Edge e = *it;
      if (symmetrize_ || (!symmetrize_ && !transpose))
        fetch_and_add(degrees[e.u], 1);
      if (symmetrize_ || (!symmetrize_ && transpose))
        fetch_and_add(degrees[(NodeID_) e.v], 1);
    }
    return degrees;
  }

  static
  pvector<SGOffset> PrefixSum(const pvector<NodeID_> &degrees) {
    pvector<SGOffset> sums(degrees.size() + 1);
    SGOffset total = 0;
    for (size_t n=0; n < degrees.size(); n++) {
      sums[n] = total;
      total += degrees[n];
    }
    sums[degrees.size()] = total;
    return sums;
  }

  static
  pvector<SGOffset> ParallelPrefixSum(const pvector<NodeID_> &degrees) {
    const size_t block_size = 1<<20;
    const size_t num_blocks = (degrees.size() + block_size - 1) / block_size;
    pvector<SGOffset> local_sums(num_blocks);
    #pragma omp parallel for
    for (size_t block=0; block < num_blocks; block++) {
      SGOffset lsum = 0;
      size_t block_end = std::min((block + 1) * block_size, degrees.size());
      for (size_t i=block * block_size; i < block_end; i++)
        lsum += degrees[i];
      local_sums[block] = lsum;
    }
    pvector<SGOffset> bulk_prefix(num_blocks+1);
    SGOffset total = 0;
    for (size_t block=0; block < num_blocks; block++) {
      bulk_prefix[block] = total;
      total += local_sums[block];
    }
    bulk_prefix[num_blocks] = total;
    pvector<SGOffset> prefix(degrees.size() + 1);
    #pragma omp parallel for
    for (size_t block=0; block < num_blocks; block++) {
      SGOffset local_total = bulk_prefix[block];
      size_t block_end = std::min((block + 1) * block_size, degrees.size());
      for (size_t i=block * block_size; i < block_end; i++) {
        prefix[i] = local_total;
        local_total += degrees[i];
      }
    }
    prefix[degrees.size()] = bulk_prefix[num_blocks];
    return prefix;
  }

  // Removes self-loops and redundant edges
  // Side effect: neighbor IDs will be sorted
  void SquishCSR(const CSRGraph<NodeID_, DestID_, invert> &g, bool transpose,
                 DestID_*** sq_index, DestID_** sq_neighs) {
    pvector<NodeID_> diffs(g.num_nodes());
    DestID_ *n_start, *n_end;
    #pragma omp parallel for private(n_start, n_end)
    for (NodeID_ n=0; n < g.num_nodes(); n++) {
      if (transpose) {
        n_start = g.in_neigh(n).begin();
        n_end = g.in_neigh(n).end();
      } else {
        n_start = g.out_neigh(n).begin();
        n_end = g.out_neigh(n).end();
      }
      std::sort(n_start, n_end);
      DestID_ *new_end = std::unique(n_start, n_end);
      new_end = std::remove(n_start, new_end, n);
      diffs[n] = new_end - n_start;
    }
    pvector<SGOffset> sq_offsets = ParallelPrefixSum(diffs);
    *sq_neighs = new DestID_[sq_offsets[g.num_nodes()]];
    *sq_index = CSRGraph<NodeID_, DestID_>::GenIndex(sq_offsets, *sq_neighs);
    #pragma omp parallel for private(n_start)
    for (NodeID_ n=0; n < g.num_nodes(); n++) {
      if (transpose)
        n_start = g.in_neigh(n).begin();
      else
        n_start = g.out_neigh(n).begin();
      std::copy(n_start, n_start+diffs[n], (*sq_index)[n]);
    }
  }

  CSRGraph<NodeID_, DestID_, invert> SquishGraph(
      const CSRGraph<NodeID_, DestID_, invert> &g) {
    DestID_ **out_index, *out_neighs, **in_index, *in_neighs;
    SquishCSR(g, false, &out_index, &out_neighs);
    if (g.directed()) {
      if (invert)
        SquishCSR(g, true, &in_index, &in_neighs);
      return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), out_index,
                                                out_neighs, in_index,
                                                in_neighs);
    } else {
      return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), out_index,
                                                out_neighs);
    }
  }

  /*
  Graph Bulding Steps (for CSR):
    - Read edgelist once to determine vertex degrees (CountDegrees)
    - Determine vertex offsets by a prefix sum (ParallelPrefixSum)
    - Allocate storage and set points according to offsets (GenIndex)
    - Copy edges into storage
  */
  void MakeCSR(const EdgeList &el, bool transpose, DestID_*** index,
               DestID_** neighs) {
    pvector<NodeID_> degrees = CountDegrees(el, transpose);
    pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
    *neighs = new DestID_[offsets[num_nodes_]];
    *index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, *neighs);
    #pragma omp parallel for
    for (auto it = el.begin(); it < el.end(); it++) {
      Edge e = *it;
      if (symmetrize_ || (!symmetrize_ && !transpose))
        (*neighs)[fetch_and_add(offsets[e.u], 1)] = e.v;
      if (symmetrize_ || (!symmetrize_ && transpose))
        (*neighs)[fetch_and_add(offsets[static_cast<NodeID_>(e.v)], 1)] =
            GetSource(e);
    }
  }

  CSRGraph<NodeID_, DestID_, invert> MakeGraphFromEL(EdgeList &el) {
    DestID_ **index = nullptr, **inv_index = nullptr;
    DestID_ *neighs = nullptr, *inv_neighs = nullptr;
    Timer t;
    t.Start();
    if (num_nodes_ == -1)
      num_nodes_ = FindMaxNodeID(el)+1;
    //#if 0 //TEMP
    if (needs_weights_)
      Generator<NodeID_, DestID_, WeightT_>::InsertWeights(el);
    //#endif
    MakeCSR(el, false, &index, &neighs);
    if (!symmetrize_ && invert)
      MakeCSR(el, true, &inv_index, &inv_neighs);
    t.Stop();
    PrintTime("Build Time", t.Seconds());
    if (symmetrize_)
      return CSRGraph<NodeID_, DestID_, invert>(num_nodes_, index, neighs);
    else
      return CSRGraph<NodeID_, DestID_, invert>(num_nodes_, index, neighs,
                                                inv_index, inv_neighs);
  }
  
  #if 0 //will complete the code later
  CSRGraph<NodeID_, DestID_, invert> relabelForSpatialLocality(
      const CSRGraph<NodeID_, DestID_, invert> &g) {
    if (g.directed()) {
        // Will add support soon
    }
    else {
        Timer t;
        t.start();
        
        /* STEP I: make a map between new and old vertex labels */

        long long counter = 0; //keep track of local counts
        std::map<NodeID_, int64_t> reMap[128]; //Conservatively assuming we will never use more than 128 threads

        /* relabel vertices in parallel (using local counter) */
        #pragma omp parallel for firstprivate(count)
        for (NodeID_ v = 0; v < g.num_nodes(); v++) {
            if (reMap[omp_get_thread_num()].find(v) == reMap.end()) {
                // vertex hasn't been labelled
                reMap.insert(std::pair<NodeID_, int64_t>(v, counter));
                counter++;
            }
            for (NodeID_ u : g.in_neigh(v)) {
                if (reMap[omp_get_thread_num()].find(u) == reMap.end()) {
                    // vertex hasn't been labelled
                    reMap.insert(std::pair<NodeID_, int64_t>(u, counter));
                    counter++;
                }
            }
        }

        /* Update counts based on maximum count for each thread */
        int64_t offset = 0;
        for (int i = 0; i < 128; i++) {
            if (reMap[i].size() != 0) {
                // adding offset to all counts of current map
                std::map<NodeID_, int64_t>::iterator it, it_end;
                #pragma omp parallel for 
                for (it = reMap[i].begin(), it_end = reMap[i].end(); it != it_end; it++) {
                    it->second += offset; 
                }
                
                // finding maximum value of current set 
                int64_t maxVal = 0;
                #pragma omp parallel for reduction(max: maxVal)
                for (it = reMap[i].begin(), it_end = reMap[i].end(); it != it_end; it++) {
                    if (it->second > maxVal) {
                        maxVal = it->second;
                    }
                }
                offset = maxVal;
            }
        }
        
        /* Merge local containers */
        std::map <NodeID_, int64_t> merged_reMap; 
        for (int i = 0; i < 128; i++) {
            if (reMap[i].size() != 0) {
                merged_reMap.insert(reMap[i].begin(), reMap[i].end());
            }
        }

        /* STEP II: rewrite CSR based on this reMap */
        DestID_* neighs = new DestID_[2 * g.num_edges()];
        DestID_** index = CSRGraph<NodeID_, DestID_>::relabelIndex(offsets, neighs, reMap);


    }
  }
  #endif

  CSRGraph<NodeID_, DestID_, invert> MakeGraph() {
    CSRGraph<NodeID_, DestID_, invert> g;
    {  // extra scope to trigger earlier deletion of el (save memory)
      EdgeList el;
      if (cli_.filename() != "") {
        Reader<NodeID_, DestID_, WeightT_, invert> r(cli_.filename());
        if ((r.GetSuffix() == ".sg") || (r.GetSuffix() == ".wsg")) {
          return r.ReadSerializedGraph();
        } else {
          el = r.ReadFile(needs_weights_);
        }
      } else if (cli_.scale() != -1) {
        Generator<NodeID_, DestID_> gen(cli_.scale(), cli_.degree());
        el = gen.GenerateEL(cli_.uniform());
      }
      g = MakeGraphFromEL(el);
    }
    #if 0
    if (cli_.relabel() == 1) {
        g_new = relabelForSpatialLocality(g); 
    }
    #endif
    return SquishGraph(g);
  }

  // Relabels (and rebuilds) graph by order of decreasing degree
  static
  CSRGraph<NodeID_, DestID_, invert> RelabelByDegree(
      const CSRGraph<NodeID_, DestID_, invert> &g) {
    if (g.directed()) {
      std::cout << "Cannot relabel directed graph" << std::endl;
      std::exit(-11);
    }
    Timer t;
    t.Start();
    typedef std::pair<int64_t, NodeID_> degree_node_p;
    pvector<degree_node_p> degree_id_pairs(g.num_nodes());
    #pragma omp parallel for
    for (NodeID_ n=0; n < g.num_nodes(); n++)
      degree_id_pairs[n] = std::make_pair(g.out_degree(n), n);
    std::sort(degree_id_pairs.begin(), degree_id_pairs.end(),
              std::greater<degree_node_p>());
    pvector<NodeID_> degrees(g.num_nodes());
    pvector<NodeID_> new_ids(g.num_nodes());
    #pragma omp parallel for
    for (NodeID_ n=0; n < g.num_nodes(); n++) {
      degrees[n] = degree_id_pairs[n].first;
      new_ids[degree_id_pairs[n].second] = n;
    }
    pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
    DestID_* neighs = new DestID_[offsets[g.num_nodes()]];
    DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
    #pragma omp parallel for schedule (dynamic, 1024)
    for (NodeID_ u=0; u < g.num_nodes(); u++) {
      for (NodeID_ v : g.out_neigh(u))
        neighs[offsets[new_ids[u]]++] = new_ids[v];
      std::sort(index[new_ids[u]], index[new_ids[u]+1]);
    }
    t.Stop();
    PrintTime("Relabel", t.Seconds());
    return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs);
  }
 
  // Frequency-based Clustering as described in the CAGRA paper [BigData 17]
  static
  CSRGraph<NodeID_, DestID_, invert> degreeCluster(
      const CSRGraph<NodeID_, DestID_, invert> &g, bool outDegree, pvector<NodeID_> &new_ids, bool createOnlyDegList, bool createBothCSRs, int cutoffFactor) {
      Timer t;
      t.Start();
       
      typedef std::pair<int64_t, NodeID_> degree_node_p;
      pvector<degree_node_p> degree_id_pairs(g.num_nodes());
      NodeID_ avgDegree    = g.num_edges_directed() / g.num_nodes();
      NodeID_ degreeCutoff = avgDegree * cutoffFactor;
      NodeID_ hubCount(0);
      if (g.directed() == true) {
        /* Step I: Create a list of degrees */
        #pragma omp parallel for reduction(+ : hubCount)
        for (NodeID_ n=0; n < g.num_nodes(); n++) { 
          degree_id_pairs[n] = std::make_pair(g.out_degree(n), n);
          if (outDegree == true) {
            if (g.out_degree(n) > degreeCutoff) {
              ++hubCount;
            }
          }
          else {
            if (g.in_degree(n) > degreeCutoff) {
              ++hubCount;
            }
          }
        }

        /* Step II: Sort based on degree order */
        __gnu_parallel::sort(degree_id_pairs.begin(), degree_id_pairs.end(),
                             std::greater<degree_node_p>()); //TODO:Use parallel sort

        /* Step III: assigned remap for the hub vertices */
        #pragma omp parallel for 
        for (NodeID_ n=0; n < hubCount; n++) {
          new_ids[degree_id_pairs[n].second] = n;
        }

        /* Step IV: assigning remap for (easy) non hub vertices */
        NodeID_ numHubs = hubCount;
        SlidingQueue<NodeID_> queue(numHubs);
        int numThreads = omp_get_max_threads();
        #pragma omp parallel 
        {
            QueueBuffer<NodeID_> lqueue(queue, numHubs / numThreads); 
            #pragma omp for
            for (NodeID_ n = numHubs; n < g.num_nodes(); ++n) {
                if (new_ids[n] == -1) {
                    new_ids[n] = n;
                }
                else {
                    NodeID_ remappedTo  = new_ids[n];
                    if (new_ids[remappedTo] == -1) {
                        new_ids[remappedTo] = n; //Swap ids
                    }
                    else {
                        lqueue.push_back(n);
                    }
                }
            } 
            lqueue.flush();
        }
        queue.slide_window();
        

        /* Step V: assigning remaps for remaining non hubs */
        NodeID_ unassignedCtr = 0;
        auto q_iter = queue.begin();
        #pragma omp parallel for
        for (NodeID_ n = 0; n < numHubs; ++n) {
            if (new_ids[n] == -1) {
                NodeID_ u = *(q_iter + __sync_fetch_and_add(&unassignedCtr, 1));
                new_ids[n] = u; 
            }
        }
        
        /* Step VI: generate degree to build a new graph */
        pvector<NodeID_> degrees(g.num_nodes());
        pvector<NodeID_> inv_degrees(g.num_nodes());
        if (outDegree == true) {
          #pragma omp parallel for 
          for (NodeID_ n=0; n < g.num_nodes(); n++) {
            degrees[new_ids[n]]     = g.out_degree(n); 
            inv_degrees[new_ids[n]] = g.in_degree(n);
          }
        }
        else {
          #pragma omp parallel for 
          for (NodeID_ n=0; n < g.num_nodes(); n++) {
            degrees[new_ids[n]]     = g.in_degree(n); 
            inv_degrees[new_ids[n]] = g.out_degree(n);
          }
        }
        
        /* Graph building phase */
		pvector<SGOffset> offsets = ParallelPrefixSum(inv_degrees);
        DestID_* neighs = new DestID_[offsets[g.num_nodes()]];
        DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
        #pragma omp parallel for schedule (dynamic, 1024)
        for (NodeID_ u=0; u < g.num_nodes(); u++) {
          if (outDegree == true) {
            for (NodeID_ v : g.in_neigh(u))
              neighs[offsets[new_ids[u]]++] = new_ids[v];
          }
          else {
            for (NodeID_ v : g.out_neigh(u))
              neighs[offsets[new_ids[u]]++] = new_ids[v];
          }
        }
        DestID_* inv_neighs(nullptr);
        DestID_** inv_index(nullptr);
        if (createOnlyDegList == true || createBothCSRs == true) {
          // making the inverse list (in-degrees in this case)
          pvector<SGOffset> inv_offsets = ParallelPrefixSum(degrees);
          inv_neighs = new DestID_[inv_offsets[g.num_nodes()]];
          inv_index = CSRGraph<NodeID_, DestID_>::GenIndex(inv_offsets, inv_neighs);
          if (createBothCSRs == true) {
            #pragma omp parallel for schedule(dynamic, 1024)
            for (NodeID_ u=0; u < g.num_nodes(); u++) {
              if (outDegree == true) {
                for (NodeID_ v : g.out_neigh(u))
                  inv_neighs[inv_offsets[new_ids[u]]++] = new_ids[v];
              }
              else {
                for (NodeID_ v : g.in_neigh(u))
                  inv_neighs[inv_offsets[new_ids[u]]++] = new_ids[v];
              }
            }
          }
        }
        t.Stop();
        PrintTime("HubSort time", t.Seconds());
        if (outDegree == true) {
          return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), inv_index, inv_neighs, index, neighs);
        }
        else {
          return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs, inv_index, inv_neighs);
        }
      }
      else {
        /* Undirected graphs - no need to make separate lists for in and out degree */
        /* Step I: Create a list of degrees */
        #pragma omp parallel for reduction(+ : hubCount)
        for (NodeID_ n=0; n < g.num_nodes(); n++) { 
          degree_id_pairs[n] = std::make_pair(g.out_degree(n), n);
          if (g.out_degree(n) > degreeCutoff) {
            ++hubCount;
          }
        }

        /* Step II: Sort based on degree order */
        __gnu_parallel::sort(degree_id_pairs.begin(), degree_id_pairs.end(),
                             std::greater<degree_node_p>()); //TODO:Use parallel sort

        /* Step III: assigned remap for the hub vertices */
        #pragma omp parallel for 
        for (NodeID_ n=0; n < hubCount; n++) {
          new_ids[degree_id_pairs[n].second] = n;
        }

        /* Step IV: assigning remap for (easy) non hub vertices */
        NodeID_ numHubs = hubCount;
        SlidingQueue<NodeID_> queue(numHubs);
        int numThreads = omp_get_max_threads();
        #pragma omp parallel 
        {
            QueueBuffer<NodeID_> lqueue(queue, numHubs / numThreads);
            #pragma omp for
            for (NodeID_ n = numHubs; n < g.num_nodes(); ++n) {
                if (new_ids[n] == -1) {
                    new_ids[n] = n;
                }
                else {
                    NodeID_ remappedTo  = new_ids[n];
                    if (new_ids[remappedTo] == -1) {
                        new_ids[remappedTo] = n; //Swap ids
                    }
                    else {
                        lqueue.push_back(n);
                    }
                }
            } 
            lqueue.flush();
        }
        queue.slide_window();
        

        /* Step V: assigning remaps for remaining non hubs */
        NodeID_ unassignedCtr = 0;
        auto q_iter = queue.begin();
        #pragma omp parallel for
        for (NodeID_ n = 0; n < numHubs; ++n) {
            if (new_ids[n] == -1) {
                NodeID_ u = *(q_iter + __sync_fetch_and_add(&unassignedCtr, 1));
                new_ids[n] = u; 
            }
        }
        
        /* Step VI: generate degree to build a new graph */
        pvector<NodeID_> degrees(g.num_nodes());
        #pragma omp parallel for 
        for (NodeID_ n=0; n < g.num_nodes(); n++) {
          degrees[new_ids[n]] = g.out_degree(n); 
        }
        
        /* Graph building phase */
        pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
        DestID_* neighs = new DestID_[offsets[g.num_nodes()]];
        DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
        #pragma omp parallel for schedule (dynamic, 1024)
        for (NodeID_ u=0; u < g.num_nodes(); u++) {
          for (NodeID_ v : g.out_neigh(u))
            neighs[offsets[new_ids[u]]++] = new_ids[v];
        }
        t.Stop();
        PrintTime("HubSort time", t.Seconds());
        return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs);
      }
    }
  
  // Similar to the previous function but handles directed graphs,
  // weighted graphs, and also does relabeling by user specified method
  static
  CSRGraph<NodeID_, DestID_, invert> degreeCluster_weighted(
      const CSRGraph<NodeID_, DestID_, invert> &g, bool outDegree, pvector<NodeID_> &new_ids, bool createOnlyDegList, bool createBothCSRs, int cutoffFactor) {
      Timer t;
      t.Start();
       
      typedef std::pair<int64_t, NodeID_> degree_node_p;
      pvector<degree_node_p> degree_id_pairs(g.num_nodes());
      NodeID_ avgDegree    = g.num_edges_directed() / g.num_nodes();
      NodeID_ degreeCutoff = avgDegree * cutoffFactor;
      NodeID_ hubCount(0);
      if (g.directed() == true) {
        /* Step I: Create a list of degrees */
        #pragma omp parallel for reduction(+ : hubCount)
        for (NodeID_ n=0; n < g.num_nodes(); n++) { 
          degree_id_pairs[n] = std::make_pair(g.out_degree(n), n);
          if (outDegree == true) {
            if (g.out_degree(n) > degreeCutoff) {
              ++hubCount;
            }
          }
          else {
            if (g.in_degree(n) > degreeCutoff) {
              ++hubCount;
            }
          }
        }
        
        /* Step II: Sort based on degree order */
        __gnu_parallel::sort(degree_id_pairs.begin(), degree_id_pairs.end(),
                             std::greater<degree_node_p>()); //TODO:Use parallel sort

        /* Step III: assigned remap for the hub vertices */
        #pragma omp parallel for 
        for (NodeID_ n=0; n < hubCount; n++) {
          new_ids[degree_id_pairs[n].second] = n;
        }

        /* Step IV: assigning remap for (easy) non hub vertices */
        NodeID_ numHubs = hubCount;
        SlidingQueue<NodeID_> queue(numHubs);
        int numThreads = omp_get_max_threads();
        #pragma omp parallel 
        {
            QueueBuffer<NodeID_> lqueue(queue, numHubs / numThreads); 
            #pragma omp for
            for (NodeID_ n = numHubs; n < g.num_nodes(); ++n) {
                if (new_ids[n] == -1) {
                    new_ids[n] = n;
                }
                else {
                    NodeID_ remappedTo  = new_ids[n];
                    if (new_ids[remappedTo] == -1) {
                        new_ids[remappedTo] = n; //Swap ids
                    }
                    else {  
                        lqueue.push_back(n);
                    }
                }
            } 
            lqueue.flush();
        }
        queue.slide_window();
        

        /* Step V: assigning remaps for remaining non hubs */
        NodeID_ unassignedCtr = 0;
        auto q_iter = queue.begin();
        #pragma omp parallel for
        for (NodeID_ n = 0; n < numHubs; ++n) {
            if (new_ids[n] == -1) {
                NodeID_ u = *(q_iter + __sync_fetch_and_add(&unassignedCtr, 1));
                new_ids[n] = u; 
            }
        }
        
        /* Step VI: generate degree to build a new graph */
        pvector<NodeID_> degrees(g.num_nodes());
        pvector<NodeID_> inv_degrees(g.num_nodes());
        if (outDegree == true) {
          #pragma omp parallel for 
          for (NodeID_ n=0; n < g.num_nodes(); n++) {
            degrees[new_ids[n]]     = g.out_degree(n); 
            inv_degrees[new_ids[n]] = g.in_degree(n);
          }
        }
        else {
          #pragma omp parallel for 
          for (NodeID_ n=0; n < g.num_nodes(); n++) {
            degrees[new_ids[n]]     = g.in_degree(n); 
            inv_degrees[new_ids[n]] = g.out_degree(n);
          }
        }
        
        /* Graph building phase */
        pvector<SGOffset> offsets = ParallelPrefixSum(inv_degrees);
        DestID_* neighs = new DestID_[offsets[g.num_nodes()]];
        DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
        #pragma omp parallel for schedule(dynamic, 1024)
        for (NodeID_ u=0; u < g.num_nodes(); u++) {
          if (outDegree) {
            for (auto v : g.in_neigh(u)) {
              auto oldWeight = v.w;
              auto newID = new_ids[(NodeID_) v.v];
              DestID_ newV (newID, oldWeight);
              neighs[offsets[new_ids[u]]++] = newV;
            }
          }
          else {
            for (auto v : g.out_neigh(u)) {
              auto oldWeight = v.w;
              auto newID = new_ids[(NodeID_) v.v];
              DestID_ newV (newID, oldWeight);
              neighs[offsets[new_ids[u]]++] = newV;
            }
          }
        }
        DestID_* inv_neighs(nullptr);
        DestID_** inv_index(nullptr);
        if (createBothCSRs == true || createOnlyDegList == true) {
          // making the inverse list (in-degrees in this case)
          pvector<SGOffset> inv_offsets = ParallelPrefixSum(degrees);
          inv_neighs = new DestID_[inv_offsets[g.num_nodes()]];
          inv_index = CSRGraph<NodeID_, DestID_>::GenIndex(inv_offsets, inv_neighs);
        
          if (createBothCSRs == true) {
            #pragma omp parallel for  schedule(dynamic, 1024)
            for (NodeID_ u=0; u < g.num_nodes(); u++) {
              if (outDegree == true) {
                for (auto v : g.out_neigh(u)) {
                  auto oldWeight = v.w;
                  auto newID = new_ids[(NodeID_) v.v];
                  DestID_ newV (newID, oldWeight);
                  inv_neighs[inv_offsets[new_ids[u]]++] = newV;
                }
              }
              else {
                for (auto v : g.in_neigh(u)) {
                  auto oldWeight = v.w;
                  auto newID = new_ids[(NodeID_) v.v];
                  DestID_ newV (newID, oldWeight);
                  inv_neighs[inv_offsets[new_ids[u]]++] = newV;
                }
              }
            }
          }
        }
        t.Stop();
        PrintTime("HubSort time", t.Seconds());
        if (outDegree == true) {
          return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), inv_index, inv_neighs, index, neighs);
        }
        else {
          return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs, inv_index, inv_neighs);
        }
      }
      else {
        /* Undirected graphs - no need to make separate lists for in and out degree */
        /* Step I: Create a list of degrees */
        #pragma omp parallel for reduction(+ : hubCount)
        for (NodeID_ n=0; n < g.num_nodes(); n++) { 
          degree_id_pairs[n] = std::make_pair(g.out_degree(n), n);
          if (g.out_degree(n) > degreeCutoff) {
            ++hubCount;
          }
        }

        /* Step II: Sort based on degree order */
        __gnu_parallel::sort(degree_id_pairs.begin(), degree_id_pairs.end(),
                             std::greater<degree_node_p>()); //TODO:Use parallel sort

        /* Step III: assigned remap for the hub vertices */
        #pragma omp parallel for 
        for (NodeID_ n=0; n < hubCount; n++) {
          new_ids[degree_id_pairs[n].second] = n;
        }

        /* Step IV: assigning remap for (easy) non hub vertices */
        NodeID_ numHubs = hubCount;
        SlidingQueue<NodeID_> queue(numHubs);
        int numThreads = omp_get_max_threads();
        #pragma omp parallel 
        {
            QueueBuffer<NodeID_> lqueue(queue, numHubs / numThreads);
            #pragma omp for
            for (NodeID_ n = numHubs; n < g.num_nodes(); ++n) {
                if (new_ids[n] == -1) {
                    new_ids[n] = n;
                }
                else {
                    NodeID_ remappedTo  = new_ids[n];
                    if (new_ids[remappedTo] == -1) {
                        new_ids[remappedTo] = n; //Swap ids
                    }
                    else {
                        lqueue.push_back(n);
                    }
                }
            } 
            lqueue.flush();
        }
        queue.slide_window();
        



        /* Step V: assigning remaps for remaining non hubs */
        NodeID_ unassignedCtr = 0;
        auto q_iter = queue.begin();
        #pragma omp parallel for
        for (NodeID_ n = 0; n < numHubs; ++n) {
            if (new_ids[n] == -1) {
                NodeID_ u = *(q_iter + __sync_fetch_and_add(&unassignedCtr, 1));
                new_ids[n] = u; 
            }
        }
        
        /* Step VI: generate degree to build a new graph */
        pvector<NodeID_> degrees(g.num_nodes());
        #pragma omp parallel for 
        for (NodeID_ n=0; n < g.num_nodes(); n++) {
          degrees[new_ids[n]] = g.out_degree(n); 
        }
        
        /* Graph building phase */
        pvector<SGOffset> offsets = ParallelPrefixSum(degrees);
        DestID_* neighs = new DestID_[offsets[g.num_nodes()]];
        DestID_** index = CSRGraph<NodeID_, DestID_>::GenIndex(offsets, neighs);
        #pragma omp parallel for schedule (dynamic, 1024)
        for (NodeID_ u=0; u < g.num_nodes(); u++) {
          for (auto v : g.out_neigh(u)) {
            auto oldWeight = v.w;
            auto newID = new_ids[(NodeID_) v.v];
            DestID_ newV (newID, oldWeight);
            neighs[offsets[new_ids[u]]++] = newV;
          }
        }
        t.Stop();
        PrintTime("HubSort time", t.Seconds());
        return CSRGraph<NodeID_, DestID_, invert>(g.num_nodes(), index, neighs);
      }
    }
};

#endif  // BUILDER_H_
