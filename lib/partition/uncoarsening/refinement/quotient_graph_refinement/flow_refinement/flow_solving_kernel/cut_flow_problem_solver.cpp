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
                                             flow_graph & fG) {

        //building up the graph as in parse.h of hi_pr code
        NodeID idx = 0;
        new_to_old_ids.resize(lhs_boundary_stripe.size() + rhs_boundary_stripe.size());
        std::unordered_map<NodeID, NodeID> old_to_new;
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

        regions_no_edges(G, lhs_boundary_stripe, rhs_boundary_stripe, 
                         lhs, rhs, outer_lhs_boundary, outer_rhs_boundary);
        
        if(outer_lhs_boundary.size() == 0 || outer_rhs_boundary.size() == 0) return false;
        NodeID n = lhs_boundary_stripe.size() + rhs_boundary_stripe.size() + 2; //+source and target
        fG.start_construction(n);

        NodeID source = n-2;
        NodeID sink   = n-1;
        idx = 0;

        //add LHS stripe to flow problem
        for( unsigned i = 0; i < lhs_boundary_stripe.size(); i++, idx++) {
                NodeID node = lhs_boundary_stripe[i];
                NodeID sourceID = idx;
                forall_out_edges(G, e, node) {
                        if(G.getPartitionIndex(G.getEdgeTarget(e)) == BOUNDARY_STRIPE_NODE)  {
                                NodeID targetID     = old_to_new[G.getEdgeTarget(e)];
                                fG.new_edge(sourceID, targetID, G.getEdgeWeight(e));
                        }
                } endfor
        }

        //add RHS stripe to flow problem
        for( unsigned i = 0; i < rhs_boundary_stripe.size(); i++, idx++) {
                NodeID node = rhs_boundary_stripe[i];
                NodeID sourceID = idx;
                forall_out_edges(G, e, node) {
                        if(G.getPartitionIndex(G.getEdgeTarget(e)) == BOUNDARY_STRIPE_NODE)  {
                                NodeID targetID     = old_to_new[G.getEdgeTarget(e)];
                                fG.new_edge(sourceID, targetID, G.getEdgeWeight(e));
                        }
                } endfor
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
        bool do_sth = convert_ds(config, G, lhs, rhs, lhs_boundary_stripe, rhs_boundary_stripe, new_to_old_ids, fG );

        if(!do_sth) return initial_cut;

        push_relabel pr;
        NodeID source = fG.number_of_nodes()-2;
        NodeID sink   = fG.number_of_nodes()-1;
        std::vector< NodeID > source_set;
        FlowType flowvalue = pr.solve_max_flow_min_cut( fG, source, sink, true, source_set);

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

                forall_nodes(fG, node) {
                        NodeID node = residualGraph.new_node(); // for each node here create a new node 

                        if( node < fG.number_of_nodes() -2 ) {
                                residualGraph.setNodeWeight( node, G.getNodeWeight(new_to_old_ids[node]));
                        }

                        forall_out_edges(fG, e, node) {
                                NodeID target = fG.getEdgeTarget(node, e);
                                if( fG.getEdgeCapacity(node, e) > 0 ) {
                                        if( fG.getEdgeFlow(node, e) < (FlowType) fG.getEdgeCapacity(node, e)) {
                                                residualGraph.new_edge(node, target);
                                        } else {
                                                forall_out_edges(fG, e_bar, target) {
                                                        NodeID target_prime = fG.getEdgeTarget(target, e_bar);
                                                        if( target_prime == node && fG.getEdgeFlow(target, e_bar) > 0) {
                                                                residualGraph.new_edge(node, target);
                                                        }
                                                } endfor
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

