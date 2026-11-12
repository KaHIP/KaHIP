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

#include "algorithms/push_relabel.h"
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
                                fG.new_edge(sourceID, targetID, G.getEdgeWeight(e));
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
                                fG.new_edge(sourceID, targetID, G.getEdgeWeight(e));
                        } else is_outer_boundary = true;
                } endfor
                if(is_outer_boundary) {
                        outer_rhs_boundary.push_back(idx);
                }
        }

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
        // Without most balanced cuts only the set of nodes reachable from the source in the
        // residual graph of a maximum flow is used; it is the same for every maximum flow. It is
        // computed on the reversed network: there it is the set of nodes that can reach the
        // (original) source, which a maximum preflow from the original sink already gives, so the
        // excess return of the second phase is skipped.
        const bool reversed = !config.most_balanced_minimum_cuts;
        fG.set_reversed(reversed);
        FlowType lhs_rhs_capacity = 0;
        bool do_sth = convert_ds(config, G, lhs, rhs, lhs_boundary_stripe, rhs_boundary_stripe, new_to_old_ids, fG, lhs_rhs_capacity );

        if(!do_sth) return initial_cut;

        push_relabel pr;
        NodeID source = fG.number_of_nodes()-2;
        NodeID sink   = fG.number_of_nodes()-1;
        std::vector< NodeID > source_set;

        // Without most balanced cuts the caller accepts a result only if its cut is below
        // initial_cut; a flow of value initial_cut ends the refinement of this pair with every
        // change undone, whatever source set comes with it. So if the cut {source} + LHS stripe
        // has exactly that capacity, the solver can stop as soon as the sink has received it.
        FlowType stop_flow = std::numeric_limits<FlowType>::max();
        if(!config.most_balanced_minimum_cuts && lhs_rhs_capacity == (FlowType) initial_cut) {
                stop_flow = lhs_rhs_capacity;
        }
        FlowType flowvalue = reversed
                ? pr.solve_max_flow_min_cut( fG, sink, source, true, source_set, stop_flow, true)
                : pr.solve_max_flow_min_cut( fG, source, sink, true, source_set, stop_flow);
        if(flowvalue >= stop_flow) return initial_cut;

        std::vector< bool > new_rhs_flag(fG.number_of_nodes(), true);
        for( unsigned int i = 0; i < source_set.size(); i++) {
                new_rhs_flag[source_set[i]] = false;
        }

        if(!config.most_balanced_minimum_cuts) {
                forall_nodes(fG, node) {
                        if(new_rhs_flag[node] && node < fG.number_of_nodes()-2) {
                                new_rhs_nodes.push_back(node);
                        }
                } endfor
        } else {
                graph_access residualGraph;
                residualGraph.start_construction(fG.number_of_nodes(),fG.number_of_edges());

                // For a saturated edge (node, target) the residual graph gets one edge per edge
                // (target, node) with positive flow. The edges (target, node) are exactly the reverse
                // partners of the edges (node, target), so their number is counted in one pass over
                // the edges of node instead of scanning the edges of target per saturated edge.
                std::vector<NodeID> pos_back(fG.number_of_nodes(), 0);
                std::vector<NodeID> pos_stamp(fG.number_of_nodes(), fG.number_of_nodes());
                forall_nodes(fG, node) {
                        NodeID node = residualGraph.new_node(); // for each node here create a new node 

                        if( node < fG.number_of_nodes() -2 ) {
                                residualGraph.setNodeWeight( node, G.getNodeWeight(new_to_old_ids[node]));
                        }

                        bool counted = false;
                        forall_out_edges(fG, e, node) {
                                NodeID target = fG.getEdgeTarget(node, e);
                                if( fG.getEdgeCapacity(node, e) > 0 ) {
                                        if( fG.getEdgeFlow(node, e) < (FlowType) fG.getEdgeCapacity(node, e)) {
                                                residualGraph.new_edge(node, target);
                                        } else {
                                                if(!counted) {
                                                        forall_out_edges(fG, e2, node) {
                                                                NodeID t2 = fG.getEdgeTarget(node, e2);
                                                                if(pos_stamp[t2] != node) { pos_stamp[t2] = node; pos_back[t2] = 0; }
                                                                if(fG.getEdgeFlow(t2, fG.getReverseEdge(node, e2)) > 0) pos_back[t2]++;
                                                        } endfor
                                                        counted = true;
                                                }
                                                for( NodeID c = 0; c < pos_back[target]; c++) {
                                                        residualGraph.new_edge(node, target);
                                                }
                                        }
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

