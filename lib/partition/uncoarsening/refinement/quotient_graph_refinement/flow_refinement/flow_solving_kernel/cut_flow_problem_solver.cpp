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
        for( unsigned i = 0; i < lhs_boundary_stripe.size(); i++) {
                G.setPartitionIndex(lhs_boundary_stripe[i], BOUNDARY_STRIPE_NODE);
                new_to_old_ids[idx]                = lhs_boundary_stripe[i];
                old_to_new[lhs_boundary_stripe[i]] = idx++ ;
        }
        for( unsigned i = 0; i < rhs_boundary_stripe.size(); i++) {
                G.setPartitionIndex(rhs_boundary_stripe[i], BOUNDARY_STRIPE_NODE); 
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
                        if(G.getPartitionIndex(G.getEdgeTarget(e)) == BOUNDARY_STRIPE_NODE)  {
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
                        if(G.getPartitionIndex(G.getEdgeTarget(e)) == BOUNDARY_STRIPE_NODE)  {
                                NodeID targetID     = old_to_new[G.getEdgeTarget(e)];
                                if(sourceID < targetID) fG.new_edge(sourceID, targetID, G.getEdgeWeight(e), G.getEdgeWeight(e));
                        } else is_outer_boundary = true;
                } endfor
                if(is_outer_boundary) {
                        outer_rhs_boundary.push_back(idx);
                }
        }

        if(outer_lhs_boundary.size() == 0 || outer_rhs_boundary.size() == 0) return false;

        // without a source or a sink edge in the unlocked model, use the locked model for this problem
        if(config.flow_unlocked_model && connect_terminals_unlocked(G, lhs, rhs, lhs_boundary_stripe, rhs_boundary_stripe,
                                                                    fG, lhs_rhs_capacity)) {
                return true;
        }

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

// Unlocked flow model (model F_H of Heuer, Sanders, Schlag, "Network Flow-Based Refinement for Multilevel
// Hypergraph Partitioning", SEA 2018, restricted to graphs): a stripe node gets an edge from the source with
// the weight of its edges to lhs nodes outside the stripes and an edge to the sink with the weight of its
// edges to rhs nodes outside the stripes; edges to third blocks are left out, since they are cut on either
// side. The value of an s-t cut is then the lhs-rhs cut of the partition it induces, and stripe nodes at the
// border of the region may change sides (the original model ties them to source or sink).
// The stripe edges are already in fG; the terminal edges are appended in stripe order, source edges first.
bool cut_flow_problem_solver::connect_terminals_unlocked(graph_access & G,
                                                         PartitionID & lhs,
                                                         PartitionID & rhs,
                                                         std::vector<NodeID> & lhs_boundary_stripe,
                                                         std::vector<NodeID> & rhs_boundary_stripe,
                                                         flow_graph & fG,
                                                         FlowType & lhs_rhs_capacity) {
        NodeID no_lhs_stripe = lhs_boundary_stripe.size();
        NodeID no_stripe     = no_lhs_stripe + rhs_boundary_stripe.size();
        NodeID source        = no_stripe;
        NodeID sink          = no_stripe + 1;

        std::vector<FlowType> to_source(no_stripe, 0);
        std::vector<FlowType> to_sink(no_stripe, 0);
        bool has_source_edge = false;
        bool has_sink_edge   = false;
        for( NodeID idx = 0; idx < no_stripe; idx++) {
                NodeID node = idx < no_lhs_stripe ? lhs_boundary_stripe[idx] : rhs_boundary_stripe[idx - no_lhs_stripe];
                forall_out_edges(G, e, node) {
                        PartitionID target_block = G.getPartitionIndex(G.getEdgeTarget(e));
                        if(target_block == lhs) {
                                to_source[idx] += G.getEdgeWeight(e);
                        } else if(target_block == rhs) {
                                to_sink[idx]   += G.getEdgeWeight(e);
                        }
                } endfor
                has_source_edge |= to_source[idx] > 0;
                has_sink_edge   |= to_sink[idx] > 0;
        }
        // without a source (sink) edge, one side of the region would be free to take all of it; the
        // caller then uses the locked model (no edge has been added yet)
        if(!has_source_edge || !has_sink_edge) return false;

        for( NodeID idx = 0; idx < no_stripe; idx++) {
                if(to_source[idx] > 0) fG.new_edge(source, idx, to_source[idx]);
                // lhs_rhs_capacity is the capacity of the cut {source} + LHS stripe
                if(idx >= no_lhs_stripe) lhs_rhs_capacity += to_source[idx];
        }
        for( NodeID idx = 0; idx < no_stripe; idx++) {
                if(to_sink[idx] > 0) fG.new_edge(idx, sink, to_sink[idx]);
                if(idx < no_lhs_stripe) lhs_rhs_capacity += to_sink[idx];
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

        if(!do_sth) {
                // no flow problem could be built: report "no change" (the rhs stripe stays rhs). With an empty
                // new_rhs_nodes the caller would move both stripes to lhs while keeping the old cut value.
                NodeID no_lhs_stripe = lhs_boundary_stripe.size();
                for( NodeID i = 0; i < rhs_boundary_stripe.size(); i++) new_rhs_nodes.push_back(no_lhs_stripe + i);
                return initial_cut;
        }

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
                // The residual graph as arrays (formerly a graph_access): an edge (node, target) of
                // positive capacity gets a residual edge if it is not saturated. For a stripe edge
                // pair with flow f this is f < w, the same relation as "(node, target) unsaturated
                // or (target, node) carries flow" in a network with one edge per direction; the
                // reverse edges of the source and sink edges (capacity 0) are left out as before.
                static std::vector<EdgeID> res_first;
                static std::vector<NodeID> res_targets;
                static std::vector<NodeWeight> res_weights;
                res_first.resize(fG.number_of_nodes()+1);
                res_weights.resize(fG.number_of_nodes());
                res_targets.clear();
                forall_nodes(fG, node) {
                        res_first[node]   = res_targets.size();
                        res_weights[node] = node < fG.number_of_nodes() - 2 ? G.getNodeWeight(new_to_old_ids[node]) : 0;
                        forall_out_edges(fG, e, node) {
                                if( fG.getEdgeCapacity(node, e) > 0 && fG.getEdgeFlow(node, e) < (FlowType) fG.getEdgeCapacity(node, e)) {
                                        res_targets.push_back(fG.getEdgeTarget(node, e));
                                }
                        } endfor
                } endfor
                res_first[fG.number_of_nodes()] = res_targets.size();
                NodeWeight average_partition_weight = ceil(config.work_load / config.k);
                NodeWeight perfect_rhs_stripe_weight = abs((int)average_partition_weight - (int)rhs_part_weight+(int) rhs_stripe_weight);
                
                most_balanced_minimum_cuts mbmc;
                mbmc.compute_good_balanced_min_cut(res_first, res_targets, res_weights, config, perfect_rhs_stripe_weight, new_rhs_nodes);
        }
        
        return flowvalue;
}

