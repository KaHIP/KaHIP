/******************************************************************************
 * cut_flow_problem_solver.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <unordered_map>
#include <cstdint>
#include <limits>

#include "algorithms/bk_maxflow.h"
#include "data_structure/graph_access.h"
#include "cut_flow_problem_solver.h"
#include "most_balanced_minimum_cuts/most_balanced_minimum_cuts.h"
#include "data_structure/flow_graph.h"
#include "io/graph_io.h"



cut_flow_problem_solver::cut_flow_problem_solver() {
}

cut_flow_problem_solver::~cut_flow_problem_solver() {
}

EdgeID cut_flow_problem_solver::regions_no_edges( graph_access & G,
                                               std::vector<NodeID> & lhs_boundary_stripe,
                                               std::vector<NodeID> & rhs_boundary_stripe,
                                               PartitionID & lhs, 
                                               PartitionID & rhs,
                                               std::vector<NodeID> & outer_lhs_boundary_nodes,
                                               std::vector<NodeID> & outer_rhs_boundary_nodes ) {

        EdgeID no_of_edges = 0;
        unsigned idx = 0;
        for( unsigned i = 0; i < lhs_boundary_stripe.size(); i++, idx++) {
                NodeID node = lhs_boundary_stripe[i];
                bool is_outer_boundary = false;
                forall_out_edges(G, e, node) {
                        if(G.getPartitionIndex(G.getEdgeTarget(e)) == BOUNDARY_STRIPE_NODE) no_of_edges++;
                        else is_outer_boundary = true;
                } endfor
                if(is_outer_boundary) {
                        outer_lhs_boundary_nodes.push_back(idx);
                }
        }

        for( unsigned i = 0; i < rhs_boundary_stripe.size(); i++, idx++) {
                NodeID node = rhs_boundary_stripe[i];
                bool is_outer_boundary = false;
                forall_out_edges(G, e, node) {
                        if(G.getPartitionIndex(G.getEdgeTarget(e)) == BOUNDARY_STRIPE_NODE) no_of_edges++;
                        else is_outer_boundary = true;
                } endfor
                if(is_outer_boundary) {
                        outer_rhs_boundary_nodes.push_back(idx);
                }
        }

        return no_of_edges; 
}

EdgeWeight cut_flow_problem_solver::convert_ds( const PartitionConfig & config, 
                                             graph_access & G, 
                                             PartitionID & lhs, 
                                             PartitionID & rhs, 
                                             std::vector<NodeID> & lhs_boundary_stripe,
                                             std::vector<NodeID> & rhs_boundary_stripe,
                                             std::vector<NodeID> & new_to_old_ids,              
                                             flow_graph & fG,
                                             FlowType & lhs_rhs_capacity) {

        //building up the graph as in parse.h of hi_pr code
        NodeID idx = 0;
        new_to_old_ids.resize(lhs_boundary_stripe.size() + rhs_boundary_stripe.size());
        // flow-graph id per stripe node, reused across calls: only stripe nodes (marked
        // BOUNDARY_STRIPE_NODE below, before any lookup) are read, so stale entries are never used
        static std::vector<NodeID> old_to_new;
        if(old_to_new.size() < G.number_of_nodes()) old_to_new.resize(G.number_of_nodes());
        // one bit per node, set exactly for the stripe nodes (the nodes marked BOUNDARY_STRIPE_NODE):
        // the scans below test their neighbors here instead of in the partition index array, which on
        // large graphs is not in the cache (hubs in the stripe have many neighbors); cleared after the scans
        static std::vector<uint64_t> in_stripe;
        if(in_stripe.size() < G.number_of_nodes()/64 + 1) in_stripe.resize(G.number_of_nodes()/64 + 1, 0);
        uint64_t * stripe_bits = in_stripe.data();
        auto set_stripe_bit = [stripe_bits](NodeID v) { stripe_bits[v >> 6] |= (uint64_t)1 << (v & 63); };
        auto is_stripe      = [stripe_bits](NodeID v) { return (stripe_bits[v >> 6] >> (v & 63)) & 1; };
        for( unsigned i = 0; i < lhs_boundary_stripe.size(); i++) {
                G.setPartitionIndex(lhs_boundary_stripe[i], BOUNDARY_STRIPE_NODE);
                set_stripe_bit(lhs_boundary_stripe[i]);
                new_to_old_ids[idx]                = lhs_boundary_stripe[i];
                old_to_new[lhs_boundary_stripe[i]] = idx++ ;
        }
        for( unsigned i = 0; i < rhs_boundary_stripe.size(); i++) {
                G.setPartitionIndex(rhs_boundary_stripe[i], BOUNDARY_STRIPE_NODE); 
                set_stripe_bit(rhs_boundary_stripe[i]);
                new_to_old_ids[idx]                = rhs_boundary_stripe[i];
                old_to_new[rhs_boundary_stripe[i]] = idx++;
        }

        std::vector<NodeID>  outer_lhs_boundary;
        std::vector<NodeID>  outer_rhs_boundary;

        // An edge {u, v} of the stripes becomes one edge pair with capacity w(u, v) in both directions
        // (added at the endpoint with the smaller id): the same flow network as an edge with a
        // reverse edge of capacity 0 per direction, with half the edges.
        // One pass per stripe node finds both its flow edges and whether it is an outer boundary node
        // (formerly a separate counting pass, regions_no_edges). The edges are added in the same order;
        // if a side has no outer boundary node the partially built graph is dropped (the caller returns).
        NodeID n = lhs_boundary_stripe.size() + rhs_boundary_stripe.size() + 2; //+source and target
        fG.start_construction(n);

        NodeID source = n-2;
        NodeID sink   = n-1;
        idx = 0;

        //add LHS stripe to flow problem
        // lhs_rhs_capacity: capacity of the edges from the LHS stripe to the RHS stripe, i.e. of the
        // cut {source} + LHS stripe, an upper bound on the maximum flow
        lhs_rhs_capacity = 0;
        NodeID no_lhs_stripe = lhs_boundary_stripe.size();
        for( unsigned i = 0; i < lhs_boundary_stripe.size(); i++, idx++) {
                NodeID node = lhs_boundary_stripe[i];
                NodeID sourceID = idx;
                bool is_outer_boundary = false;
                forall_out_edges(G, e, node) {
                        if(is_stripe(G.getEdgeTarget(e)))  {
                                NodeID targetID     = old_to_new[G.getEdgeTarget(e)];
                                if(sourceID < targetID) fG.new_edge(sourceID, targetID, G.getEdgeWeight(e), G.getEdgeWeight(e));
                                if(targetID >= no_lhs_stripe) lhs_rhs_capacity += (NodeID) G.getEdgeWeight(e);
                        } else is_outer_boundary = true;
                } endfor
                if(is_outer_boundary) {
                        outer_lhs_boundary.push_back(idx);
                }
        }

        //add RHS stripe to flow problem
        for( unsigned i = 0; i < rhs_boundary_stripe.size(); i++, idx++) {
                NodeID node = rhs_boundary_stripe[i];
                NodeID sourceID = idx;
                bool is_outer_boundary = false;
                forall_out_edges(G, e, node) {
                        if(is_stripe(G.getEdgeTarget(e)))  {
                                NodeID targetID     = old_to_new[G.getEdgeTarget(e)];
                                if(sourceID < targetID) fG.new_edge(sourceID, targetID, G.getEdgeWeight(e), G.getEdgeWeight(e));
                        } else is_outer_boundary = true;
                } endfor
                if(is_outer_boundary) {
                        outer_rhs_boundary.push_back(idx);
                }
        }

        for( unsigned i = 0; i < lhs_boundary_stripe.size(); i++) stripe_bits[lhs_boundary_stripe[i] >> 6] = 0;
        for( unsigned i = 0; i < rhs_boundary_stripe.size(); i++) stripe_bits[rhs_boundary_stripe[i] >> 6] = 0;

        if(outer_lhs_boundary.size() == 0 || outer_rhs_boundary.size() == 0) return false;

        ////connect source and target with outer boundary nodes 
        FlowType max_capacity = std::numeric_limits<FlowType>::max();
        for(unsigned i = 0; i < outer_lhs_boundary.size(); i++) {
                NodeID targetID = outer_lhs_boundary[i];
                fG.new_edge(source, targetID, max_capacity);
        }

        for(unsigned i = 0; i < outer_rhs_boundary.size(); i++) {
                NodeID sourceID = outer_rhs_boundary[i];
                fG.new_edge(sourceID, sink, max_capacity);
        }

        fG.finish_construction();
        return true;
}

EdgeWeight cut_flow_problem_solver::get_min_flow_max_cut(const PartitionConfig & config, 
                                                      graph_access & G, 
                                                      PartitionID & lhs, 
                                                      PartitionID & rhs, 
                                                      std::vector<NodeID> & lhs_boundary_stripe,
                                                      std::vector<NodeID> & rhs_boundary_stripe,
                                                      std::vector<NodeID> & new_to_old_ids,
                                                      EdgeWeight & initial_cut,
                                                      NodeWeight & rhs_part_weight,
                                                      NodeWeight & rhs_stripe_weight,
                                                      std::vector<NodeID> & new_rhs_nodes) {

        flow_graph fG;
        FlowType lhs_rhs_capacity = 0;
        bool do_sth = convert_ds(config, G, lhs, rhs, lhs_boundary_stripe, rhs_boundary_stripe, new_to_old_ids, fG, lhs_rhs_capacity );

        if(!do_sth) return initial_cut;

        NodeID source = fG.number_of_nodes()-2;
        NodeID sink   = fG.number_of_nodes()-1;

        // Without most balanced cuts the caller accepts a result only if its cut is below
        // initial_cut; a flow of value initial_cut ends the refinement of this pair with every
        // change undone, whatever source set comes with it. So if the cut {source} + LHS stripe
        // has exactly that capacity, the solver can stop as soon as the sink has received it.
        FlowType stop_flow = std::numeric_limits<FlowType>::max();
        if(!config.most_balanced_minimum_cuts && lhs_rhs_capacity == (FlowType) initial_cut) {
                stop_flow = lhs_rhs_capacity;
        }
        static bk_maxflow mf; // keeps its arrays between flow problems
        FlowType flowvalue = mf.solve( fG, source, sink, stop_flow);
        if(flowvalue >= stop_flow) return initial_cut;

        if(!config.most_balanced_minimum_cuts) {
                // the source tree is the set of nodes reachable from the source in the residual
                // graph, which is the same for every maximum flow
                forall_nodes(fG, node) {
                        if(!mf.in_source_tree(node) && node < fG.number_of_nodes()-2) {
                                new_rhs_nodes.push_back(node);
                        }
                } endfor
        } else {
                graph_access residualGraph;
                residualGraph.start_construction(fG.number_of_nodes(),fG.number_of_edges());

                // An edge (node, target) of positive capacity gets a residual edge if it is not
                // saturated. For a stripe edge pair with flow f this is f < w, the same relation as
                // "(node, target) unsaturated or (target, node) carries flow" in a network with one
                // edge per direction; the reverse edges of the source and sink edges (capacity 0)
                // are left out as before.
                forall_nodes(fG, node) {
                        NodeID node = residualGraph.new_node(); // for each node here create a new node 

                        if( node < fG.number_of_nodes() -2 ) {
                                residualGraph.setNodeWeight( node, G.getNodeWeight(new_to_old_ids[node]));
                        }

                        forall_out_edges(fG, e, node) {
                                if( fG.getEdgeCapacity(node, e) > 0 && fG.getEdgeFlow(node, e) < (FlowType) fG.getEdgeCapacity(node, e)) {
                                        residualGraph.new_edge(node, fG.getEdgeTarget(node, e));
                                }
                        } endfor
                } endfor
                
                residualGraph.setNodeWeight(source, 0);
                residualGraph.setNodeWeight(sink, 0);
                residualGraph.finish_construction();
                NodeWeight average_partition_weight = ceil(config.work_load / config.k);
                NodeWeight perfect_rhs_stripe_weight = abs((int)average_partition_weight - (int)rhs_part_weight+(int) rhs_stripe_weight);
                
                most_balanced_minimum_cuts mbmc;
                mbmc.compute_good_balanced_min_cut(residualGraph, config, perfect_rhs_stripe_weight, new_rhs_nodes);
        }
        
        return flowvalue;
}

