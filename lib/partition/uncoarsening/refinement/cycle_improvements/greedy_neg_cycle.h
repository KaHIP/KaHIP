/******************************************************************************
 * greedy_neg_cycle.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef GREEDY_NEG_CYCLE_IVBKH6WD
#define GREEDY_NEG_CYCLE_IVBKH6WD

#include "algorithms/cycle_search.h"
#include "cycle_definitions.h"
#include "data_structure/graph_access.h"
#include "definitions.h"
#include "partition_config.h"
#include "problem_factory.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.h"
#include "uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"

class greedy_neg_cycle {
public:
        greedy_neg_cycle( PartitionConfig & partition_config );
        virtual ~greedy_neg_cycle();

        EdgeWeight shortest_path_rebalance(PartitionConfig & partition_config, 
                                           graph_access & G, 
                                           complete_boundary & boundary);

        EdgeWeight negative_cycle_test(PartitionConfig & partition_config, 
                                       graph_access & G, 
                                       complete_boundary & boundary, 
                                       bool zero_weight_cycle = false);

private:
        void init_gains( PartitionConfig & partition_config, 
                         graph_access & G, 
                         graph_access & G_bar, 
                         complete_boundary & boundary, 
                         edge_movements & em); 

        void init_gains_new( PartitionConfig & partition_config, 
                             graph_access & G, 
                             graph_access & G_bar, 
                             complete_boundary & boundary, 
                             edge_movements & em); 

        problem_factory m_pf;
        kway_graph_refinement_commons* commons;

};

inline
EdgeWeight greedy_neg_cycle::shortest_path_rebalance(PartitionConfig & partition_config, 
                                                     graph_access & G, 
                                                     complete_boundary & boundary) {
        graph_access G_bar;
        boundary.getUnderlyingQuotientGraph(G_bar); 

        edge_movements em;
        init_gains_new(partition_config, G, G_bar, boundary, em); 

        NodeID s,t;
        graph_access shortest_path_graph;
        m_pf.build_shortest_path_problem(partition_config, boundary, G_bar, shortest_path_graph, s,t);

        //now we constructed the graph where we can find all negative cycles
        cycle_search cs;
        std::vector<NodeID> path;

        cs.find_shortest_path(shortest_path_graph, s,t, path);

        Gain overall_gain = 0;
        for( unsigned i = 0; i < path.size()-1; i++) {
                PartitionID lhs = path[i];
                PartitionID rhs = path[i+1];

                if(lhs == s || rhs == s || lhs == t || rhs == t) {
                        continue; // in this case it is really a path
                }

                boundary_pair bp;
                bp.k   = partition_config.k;
                bp.lhs = lhs;
                bp.rhs = rhs;

                G.setPartitionIndex(em[bp].to_move, rhs);

                boundary.postMovedBoundaryNodeUpdates(em[bp].to_move, &bp, true, true);
                boundary.setBlockNoNodes(lhs, boundary.getBlockNoNodes(lhs)-1);
                boundary.setBlockNoNodes(rhs, boundary.getBlockNoNodes(rhs)+1);
                boundary.setBlockWeight(lhs, boundary.getBlockWeight(lhs)-1);
                boundary.setBlockWeight(rhs, boundary.getBlockWeight(rhs)+1);


                ASSERT_TRUE(boundary.assert_bnodes_in_boundaries());
                ASSERT_TRUE(boundary.assert_boundaries_are_bnodes());
                overall_gain += em[bp].gain;
        } 

        path.clear();
        return overall_gain;
}


inline
EdgeWeight greedy_neg_cycle::negative_cycle_test(PartitionConfig & partition_config, 
                                                 graph_access & G, 
                                                 complete_boundary & boundary, 
                                                 bool zero_weight_cycle) {
        graph_access G_bar;
        boundary.getUnderlyingQuotientGraph(G_bar); 

        edge_movements em;
        init_gains_new(partition_config, G, G_bar, boundary, em); 
        
        NodeID s;
        graph_access cycle_graph;

        if(partition_config.kaba_include_removal_of_paths) 
                m_pf.build_cycle_problem_with_reverse(partition_config, boundary, G_bar, cycle_graph, s);
        else 
                m_pf.build_cycle_problem(partition_config, boundary, G_bar, cycle_graph, s);

        //now we constructed the graph where we can find all negative cycles
        cycle_search cs;
        std::vector<NodeID> cycle;

        bool found_some;
        if(zero_weight_cycle) {
                //this option assumes that there are no negative cycles in the graph
                found_some = cs.find_zero_weight_cycle(cycle_graph, s, cycle);
        } else {
                found_some = cs.find_negative_cycle(cycle_graph, s, cycle);
        }
        Gain overall_gain = 0;
        if(found_some) {
                for( unsigned i = 0; i < cycle.size()-1; i++) {
                        PartitionID lhs = cycle[i];
                        PartitionID rhs = cycle[i+1];

                        if(lhs == s || rhs == s) {
                                continue; // in this case it is really a path
                        }

                        boundary_pair bp;
                        bp.k   = partition_config.k;
                        bp.lhs = lhs;
                        bp.rhs = rhs;

                        G.setPartitionIndex(em[bp].to_move, rhs);
                        boundary.postMovedBoundaryNodeUpdates(em[bp].to_move, &bp, true, true);
                        boundary.setBlockNoNodes(lhs, boundary.getBlockNoNodes(lhs)-1);
                        boundary.setBlockNoNodes(rhs, boundary.getBlockNoNodes(rhs)+1);
                        boundary.setBlockWeight( lhs, boundary.getBlockWeight(lhs)-1);
                        boundary.setBlockWeight( rhs, boundary.getBlockWeight(rhs)+1);

                        ASSERT_TRUE(boundary.assert_bnodes_in_boundaries());
                        ASSERT_TRUE(boundary.assert_boundaries_are_bnodes());

                        overall_gain += em[bp].gain;
                } 
        }

        cycle.clear();
        return overall_gain;
}

inline void greedy_neg_cycle::init_gains_new( PartitionConfig & partition_config, 
                                              graph_access & G, 
                                              graph_access & G_bar, 
                                              complete_boundary & boundary, 
                                              edge_movements & em) {

        std::vector<bool>   blocked(G.number_of_nodes(), false);
        std::vector<NodeID> movements(G.number_of_edges(), 0);

        NodeID max_gainer = std::numeric_limits<NodeID>::max();
        forall_nodes(G_bar, block) {
                forall_out_edges(G_bar, e, block) {
                        NodeID target_block = G_bar.getEdgeTarget(e);
                        EdgeWeight max_gain = std::numeric_limits<EdgeWeight>::min()/100;
              
                        PartialBoundary & lhs_b = boundary.getDirectedBoundary(block, block, target_block);
                        
                        std::vector<NodeID> lhs_boundary;
                        forall_boundary_nodes(lhs_b, node) {
                                lhs_boundary.push_back(node);
                        } endfor
                        random_functions::permutate_vector_good(lhs_boundary, false);
                        for( unsigned i = 0; i < lhs_boundary.size(); i++) {
                                NodeID node = lhs_boundary[i];
                                ASSERT_EQ(G.getPartitionIndex(node), block);
                                bool valid = true;
                                forall_out_edges(G, e_bar, node) {
                                        NodeID target = G.getEdgeTarget(e_bar);
                                        if(blocked[target] == true) {
                                                valid = false;
                                                break;
                                        }

                                } endfor
                                if(blocked[node] == false && valid) {
                                        EdgeWeight int_degree = 0, ext_degree = 0;
                                        commons->int_ext_degree(G, node, block, target_block, int_degree, ext_degree); 
                                        Gain cur_gain = ext_degree - int_degree;
                                        if(cur_gain > max_gain) {
                                                max_gain   = cur_gain;
                                                max_gainer = node;
                                        }
                                }
                        } 

                        blocked[max_gainer] = true;
                        G_bar.setEdgeWeight(e,max_gain);

                        //store the information on the qgraph edge (using hashing)
                        boundary_pair bp;
                        bp.k   = partition_config.k;
                        bp.lhs = block;
                        bp.rhs = target_block;

                        em[bp].to_move = max_gainer;
                        em[bp].gain    = max_gain;

                } endfor
        } endfor



}

inline void greedy_neg_cycle::init_gains( PartitionConfig & partition_config, 
                                          graph_access & G, 
                                          graph_access & G_bar, 
                                          complete_boundary & boundary, 
                                          edge_movements & em) {

        std::vector<bool>   blocked(G.number_of_nodes(), false);
        std::vector<NodeID> movements(G.number_of_edges(), 0);
        NodeID max_gainer;

        forall_nodes(G_bar, block) {
                forall_out_edges(G_bar, e, block) {
                        NodeID target_block = G_bar.getEdgeTarget(e);
                        EdgeWeight max_gain = std::numeric_limits<EdgeWeight>::min()/100;

                        forall_nodes(G, node) {
                                bool valid = true;
                                forall_out_edges(G, e_bar, node) {
                                        NodeID target = G.getEdgeTarget(e_bar);
                                        if(blocked[target] == true) {
                                                valid = false;
                                                break;
                                        }
                                } endfor
                                if(G.getPartitionIndex(node) == block && blocked[node] == false && valid) {
                                        EdgeWeight int_degree = 0, ext_degree = 0;
                                        commons->int_ext_degree(G, node, block, target_block, int_degree, ext_degree); 
                                        Gain cur_gain = ext_degree - int_degree;
                                        if(cur_gain > max_gain) {
                                                max_gain   = cur_gain;
                                                max_gainer = node;
                                        }
                                }
                        } endfor

                        blocked[max_gainer] = true;
                        G_bar.setEdgeWeight(e,max_gain);

                        //store the information on the qgraph edge (using hashing)
                        boundary_pair bp;
                        bp.k   = partition_config.k;
                        bp.lhs = block;
                        bp.rhs = target_block;

                        em[bp].to_move = max_gainer;
                        em[bp].gain    = max_gain;

                } endfor
        } endfor



}

#endif /* end of include guard: GREEDY_NEG_CYCLE_IVBKH6WD */
