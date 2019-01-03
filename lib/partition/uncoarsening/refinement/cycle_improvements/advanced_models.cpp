/******************************************************************************
 * advanced_models.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <algorithm>

#include "advanced_models.h"
#include "algorithms/cycle_search.h"
#include "data_structure/priority_queues/bucket_pq.h"
#include "graph_io.h"
#include "partition_snapshooter.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_core.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_stop_rule.h"
#include "uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/vertex_moved_hashtable.h"

unsigned long advanced_models::conflicts = 0;


advanced_models::advanced_models() {

}

advanced_models::~advanced_models() {

}

bool advanced_models::compute_vertex_movements_rebalance_ultra( PartitionConfig & config, 
                                                                graph_access & G, 
                                                                complete_boundary & boundary, 
                                                                augmented_Qgraph & aqg, 
                                                                unsigned & steps) {


        graph_access G_bar; 
        boundary.getUnderlyingQuotientGraph(G_bar); 

        aqg.prepare(config, G, G_bar, steps);

        std::vector<bool> feasable_edge;
        std::vector<NodeID> id_mapping;
        NodeID s; NodeID t; // start vertex 
        do {
                graph_access cycle_problem;
                build_rebalance_model( config, G, G_bar, boundary, aqg, 
                                       feasable_edge, steps, cycle_problem, s, t, id_mapping);

                //********************************************************************
                //solve the problem
                //********************************************************************
                cycle_search cs;
                std::vector<NodeID> path;
                cs.find_shortest_path(cycle_problem, s, t, path);

                //detect conflict -- a block should be at most one time in a cycle
                bool conflict_detected = handle_ultra_model_conflicts(config, cycle_problem, 
                                                                      boundary, id_mapping, 
                                                                      feasable_edge, path, 
                                                                      s, aqg, true);
                if(!conflict_detected) {
                        perform_augmented_move(config, G, boundary, path, s, t, aqg);
                        return true;
                } 
        } while( true );  // at some point the model will become feasable! (otherwise the method would not have been entered and the fall back algorithm would have been applied


        return false;
}

bool advanced_models::compute_vertex_movements_rebalance( PartitionConfig & config, 
                graph_access & G, 
                complete_boundary & boundary, 
                augmented_Qgraph & aqg, 
                unsigned & steps) {

        graph_access cycle_problem;
        graph_access G_bar; 
        boundary.getUnderlyingQuotientGraph(G_bar); 

        aqg.prepare(config, G, G_bar, steps);

        //***********************************************************************
        //build the model 
        //***********************************************************************
        NodeID number_of_nodes  = G_bar.number_of_nodes()*steps + 2; 
        EdgeID number_of_edges  = G_bar.number_of_edges()*steps + 2*number_of_nodes;

        cycle_problem.start_construction(number_of_nodes, number_of_edges);
        NodeID s = number_of_nodes - 2;
        NodeID t = number_of_nodes - 1 ;

        for( unsigned s_idx = 0; s_idx < steps; s_idx++) {
                //create a new layer
                forall_nodes(G_bar, lhs) {
                        NodeID cur_node = cycle_problem.new_node();
                        forall_out_edges(G_bar, e, lhs) {
                                EdgeID rhs   = G_bar.getEdgeTarget(e);

                                //find the right edge in the augmented quotient graph
                                boundary_pair bp;
                                bp.k   = config.k;
                                bp.lhs = lhs;
                                bp.rhs = rhs;


                                unsigned load_difference = s_idx + 1;
                                if( aqg.exists_vmovements_of_diff(bp, load_difference) ) {
                                        EdgeID e_bar = cycle_problem.new_edge(cur_node, s_idx*config.k+rhs);
                                        cycle_problem.setEdgeWeight(e_bar, -aqg.get_gain_of_vmovements(bp, load_difference)); 
                                }

                        } endfor

                        //create a backward edge if it can take s_idx+1 vertices
                        if( boundary.getBlockWeight(lhs) + s_idx < config.upper_bound_partition) {
                                EdgeID e_bar = cycle_problem.new_edge(cur_node, t);
                                cycle_problem.setEdgeWeight(e_bar, 0);
                        }
                } endfor
        }

        s = cycle_problem.new_node();

        //now connect s to all vertices of the layer graph with weight zero
        for( unsigned s_idx = 0; s_idx < steps; s_idx++) {
                forall_nodes(G_bar, node) {
                        if( boundary.getBlockWeight(node) > config.upper_bound_partition ) {
                                EdgeID e = cycle_problem.new_edge(s, s_idx*config.k + node);
                                cycle_problem.setEdgeWeight(e, 0);
                        } 
                } endfor
        }

        t = cycle_problem.new_node();
        cycle_problem.finish_construction();

        //*************************************************************************************
        //solve shortest path problem in model 
        //*************************************************************************************
        //check wether t is reachable from s by performing a bfs
	std::deque<NodeID>* bfsqueue = new std::deque<NodeID>;
        std::vector<bool> touched(cycle_problem.number_of_nodes(), false);
        bfsqueue->push_back(s);
        touched[s] = true;

        cycle_search cs;
        std::vector<NodeID> path;
        cs.find_shortest_path(cycle_problem, s, t, path);

        //perform the found movements
        perform_augmented_move(config, G, boundary, path, s, t, aqg);

        return true;
}

bool advanced_models::compute_vertex_movements_ultra_model( PartitionConfig & config, 
                                                            graph_access & G, 
                                                            complete_boundary & boundary, 
                                                            augmented_Qgraph & aqg,
                                                            unsigned & steps, bool zero_weight_cycle) { 
        graph_access G_bar; 
        boundary.getUnderlyingQuotientGraph(G_bar); 

        if(!zero_weight_cycle) {
                aqg.prepare(config, G, G_bar, steps);
        }

        bool found_some;
        std::vector<bool> feasable_edge;
        std::vector<NodeID> id_mapping;
        NodeID s; // start vertex 
        do {
                graph_access cycle_problem;
                build_ultra_model( config, G, G_bar, boundary, aqg, feasable_edge, steps, cycle_problem, s, id_mapping);

                //********************************************************************
                //solve the problem
                //********************************************************************
                cycle_search cs; std::vector<NodeID> cycle;                 
                if( zero_weight_cycle ) {
                        found_some = cs.find_zero_weight_cycle(cycle_problem, s, cycle);
                } else {
                        found_some = cs.find_negative_cycle(cycle_problem, s, cycle);
                }

                if(found_some) {
                        //detect conflict -- a block should be at most one time in a cycle
                        bool conflict_detected = handle_ultra_model_conflicts(config, cycle_problem, 
                                                                              boundary, id_mapping, 
                                                                              feasable_edge, cycle, s, aqg);
                        if(!conflict_detected) {
                                perform_augmented_move(config, G, boundary, cycle, s, s, aqg);
                                return true;
                        }
                }
        } while( found_some ); 

        return false;
}



