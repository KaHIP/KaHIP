/******************************************************************************
 * problem_factory.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PROBLEM_FACTORY_KHGQXT9H
#define PROBLEM_FACTORY_KHGQXT9H

#include "definitions.h"
#include "partition_config.h"
#include "uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"

class problem_factory {
        public:
                problem_factory();
                virtual ~problem_factory();

                void build_cycle_problem( PartitionConfig & partition_config,
                                          complete_boundary & boundary, 
                                          graph_access & G_bar, 
                                          graph_access & cycle_problem, 
                                          NodeID & s); 

                void build_cycle_problem_with_reverse( PartitionConfig & partition_config, 
                                                       complete_boundary & boundary, 
                                                       graph_access & G_bar, 
                                                       graph_access & cycle_problem, 
                                                       NodeID & s);

                void build_shortest_path_problem( PartitionConfig & partition_config, 
                                                  complete_boundary & boundary, 
                                                  graph_access & G_bar, 
                                                  graph_access & shortest_path_problem, 
                                                  NodeID & s,
                                                  NodeID & t); 


};

inline
void problem_factory::build_cycle_problem( PartitionConfig & partition_config, 
                                           complete_boundary & boundary, 
                                           graph_access & G_bar, 
                                           graph_access & cycle_problem, 
                                           NodeID & s) {
        //build new graph
        cycle_problem.start_construction(G_bar.number_of_nodes()+1, G_bar.number_of_nodes() + G_bar.number_of_edges());

        forall_nodes(G_bar, node) {
                NodeID cur_node = cycle_problem.new_node();
                forall_out_edges(G_bar, e, node) {
                        EdgeID target = G_bar.getEdgeTarget(e);
                        EdgeID e_bar  = cycle_problem.new_edge(cur_node, target);
                        cycle_problem.setEdgeWeight(e_bar, -G_bar.getEdgeWeight(e));
                } endfor
        } endfor

        s = cycle_problem.new_node();
        forall_nodes(G_bar, node) {
                EdgeID e = cycle_problem.new_edge(s, node);
                cycle_problem.setEdgeWeight(e, 0);
        } endfor

        cycle_problem.finish_construction();
}

inline
void problem_factory::build_cycle_problem_with_reverse( PartitionConfig & partition_config, 
                                                        complete_boundary & boundary, 
                                                        graph_access & G_bar, 
                                                        graph_access & cycle_problem, 
                                                        NodeID & s) {
        //build new graph
        cycle_problem.start_construction(G_bar.number_of_nodes()+1, 2*G_bar.number_of_nodes() + G_bar.number_of_edges());

        s = G_bar.number_of_nodes();
        forall_nodes(G_bar, node) {
                NodeID cur_node = cycle_problem.new_node();
                forall_out_edges(G_bar, e, node) {
                        EdgeID target = G_bar.getEdgeTarget(e);
                        EdgeID e_bar  = cycle_problem.new_edge(cur_node, target);
                        cycle_problem.setEdgeWeight(e_bar, -G_bar.getEdgeWeight(e));
                } endfor
                if( boundary.getBlockWeight(cur_node) < partition_config.upper_bound_partition ) {
                        EdgeID e_bar = cycle_problem.new_edge(cur_node, s);
                        cycle_problem.setEdgeWeight(e_bar, 0);
                }
        } endfor

        NodeID real_s = cycle_problem.new_node();
        ASSERT_EQ(s, real_s);

        forall_nodes(G_bar, node) {
                EdgeID e = cycle_problem.new_edge(real_s, node);
                cycle_problem.setEdgeWeight(e, 0);
        } endfor

        cycle_problem.finish_construction();
}

inline
void problem_factory::build_shortest_path_problem( PartitionConfig & partition_config, 
                                                   complete_boundary & boundary, 
                                                   graph_access & G_bar, 
                                                   graph_access & shortest_path_problem, 
                                                   NodeID & s,
                                                   NodeID & t) {
        //build new graph
        shortest_path_problem.start_construction(G_bar.number_of_nodes()+2, 2*G_bar.number_of_nodes() + G_bar.number_of_edges());

        s = G_bar.number_of_nodes();
        t = G_bar.number_of_nodes()+1;

        forall_nodes(G_bar, node) {
                NodeID cur_node = shortest_path_problem.new_node();
                forall_out_edges(G_bar, e, node) {
                        EdgeID target = G_bar.getEdgeTarget(e);
                        EdgeID e_bar  = shortest_path_problem.new_edge(cur_node, target);
                        shortest_path_problem.setEdgeWeight(e_bar, -G_bar.getEdgeWeight(e));
                } endfor
                if( boundary.getBlockWeight(cur_node) < partition_config.upper_bound_partition ) {
                        EdgeID e_bar = shortest_path_problem.new_edge(cur_node, t);
                        shortest_path_problem.setEdgeWeight(e_bar, 0);
                }
        } endfor

        NodeID tmp = shortest_path_problem.new_node();
        ASSERT_EQ(s, tmp);
        forall_nodes(G_bar, node) {
                if( boundary.getBlockWeight(node) > partition_config.upper_bound_partition ) {
                        EdgeID e = shortest_path_problem.new_edge(tmp, node);
                        shortest_path_problem.setEdgeWeight(e, 0);
                }
        } endfor
        tmp = shortest_path_problem.new_node();
        ASSERT_EQ(t, tmp);

        shortest_path_problem.finish_construction();
}

#endif /* end of include guard: PROBLEM_FACTORY_KHGQXT9H */
