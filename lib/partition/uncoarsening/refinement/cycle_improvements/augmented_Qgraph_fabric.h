/******************************************************************************
 * augmented_Qgraph_fabric.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef AUGMENTED_QGRAPH_FABRIC_MULTITRY_FM_PVGY97EW
#define AUGMENTED_QGRAPH_FABRIC_MULTITRY_FM_PVGY97EW

#include <algorithm>
#include <vector>

#include "augmented_Qgraph.h"
#include "definitions.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.h"
#include "uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/two_way_fm.h"
#include "uncoarsening/refinement/refinement.h"

class augmented_Qgraph_fabric {
        public:
                augmented_Qgraph_fabric( );
                virtual ~augmented_Qgraph_fabric();

                //return false if the network will be feasable for the desired model
                //returns true iff rebalance = true and the fall back solution has been applied
                bool build_augmented_quotient_graph( PartitionConfig & config, 
                                                     graph_access & G, 
                                                     complete_boundary & boundary, 
                                                     augmented_Qgraph & aqg,
                                                     unsigned & s, bool rebalance, bool plus = false);

                void cleanup_eligible();

        private:
                bool construct_local_searches_on_qgraph_edge( PartitionConfig & config, 
                                                              graph_access & G, 
                                                              complete_boundary & boundary,
                                                              augmented_Qgraph & aqg,
                                                              boundary_pair & pair,
                                                              unsigned s,
                                                              bool plus);

                bool local_search(PartitionConfig & config, 
                                  bool  plus,
                                  graph_access & G, 
                                  complete_boundary & boundary,
                                  augmented_Qgraph & aqg,
                                  boundary_pair & bp,
                                  unsigned s);


                void directed_more_locallized_search(PartitionConfig & config, graph_access & G, 
                                complete_boundary & boundary, 
                                PartitionID & lhs, PartitionID & rhs,
                                NodeID start_node, unsigned & number_of_swaps, pairwise_local_search & pls);

        public:
                void more_locallized_search(PartitionConfig & config, graph_access & G, 
                                complete_boundary & boundary, 
                                PartitionID & lhs, PartitionID & rhs,
                                NodeID start_node, unsigned & number_of_swaps, pairwise_local_search & pls);
        private:
                void directed_more_locallized_search_all_bnd(PartitionConfig & config, graph_access & G, 
                                complete_boundary & boundary,  
                                PartitionID & lhs, PartitionID & rhs,
                                unsigned & number_of_swaps, pairwise_local_search & pls);

                void move_node(PartitionConfig & config, 
                               graph_access & G, 
                               NodeID & node, 
                               refinement_pq * queue, 
                               complete_boundary & boundary, 
                               PartitionID & from, 
                               PartitionID & to);

                void move_node(PartitionConfig & config, 
                               graph_access & G, 
                               NodeID & node, 
                               refinement_pq * queue, 
                               refinement_pq * to_queue, 
                               complete_boundary & boundary, 
                               PartitionID & from, 
                               PartitionID & to);

                Gain find_eligible_start_node( graph_access  & G, 
                                               PartitionID & lhs, 
                                               PartitionID & rhs, 
                                               std::vector<NodeID> & lhs_boundary, 
                                               std::vector<bool> & eligible,
                                               NodeID & start_node, bool rebalance = false); 

                void rebalance_fall_back(PartitionConfig & config, 
                                         graph_access & G,     
                                         graph_access & G_bar, 
                                         complete_boundary & boundary, 
                                         std::vector< NodeID > & candidates, 
                                         std::vector< int > & parent,
                                         augmented_Qgraph & aqg); 

                void perform_simple_move( PartitionConfig & config, 
                                          graph_access & G, 
                                          complete_boundary & boundary, 
                                          NodeID & node,
                                          PartitionID & from, 
                                          PartitionID & to);

                kway_graph_refinement_commons* commons;
                two_way_fm          m_twfm;
                std::vector<bool>   m_eligible;
                std::vector<NodeID> m_tomake_eligible;
};


inline 
Gain augmented_Qgraph_fabric::find_eligible_start_node( graph_access  & G, 
                                                         PartitionID & lhs, 
                                                         PartitionID & rhs, 
                                                         std::vector<NodeID> & lhs_boundary,  
                                                         std::vector<bool> & eligible_,
                                                         NodeID & start_node, bool rebalance) {
        //select start node
        unsigned max_idx = lhs_boundary.size();
        int random_idx   = 0;
        start_node       = lhs_boundary[0];

        Gain max_gain = std::numeric_limits<Gain>::min();
        do {
                max_gain   = std::numeric_limits<Gain>::min();
                random_idx = random_functions::nextInt(0, max_idx-1);
                for( unsigned i = 0; i < lhs_boundary.size(); i++) {
                        NodeID node = lhs_boundary[i];
                        if(eligible_[node]) {
                                EdgeWeight int_degree = 0;
                                EdgeWeight ext_degree = 0;
                                m_twfm.int_ext_degree(G, node, lhs, rhs, int_degree, ext_degree);
                                if( ext_degree - int_degree > max_gain) { //todo tiebreaking
                                        max_gain = ext_degree - int_degree;
                                } 
                        } 
                }

                if(  max_gain == std::numeric_limits<Gain>::min() ) {
                        // no node is eligible
                        break;
                }

                std::vector<NodeID> eligibles;
                for( unsigned i = 0; i < lhs_boundary.size(); i++) {
                        NodeID node = lhs_boundary[i];
                        if(eligible_[node]) {
                                EdgeWeight int_degree = 0;
                                EdgeWeight ext_degree = 0;
                                m_twfm.int_ext_degree(G, node, lhs, rhs, int_degree, ext_degree);
                                if( ext_degree - int_degree == max_gain) { 
                                        eligibles.push_back(node); 
                                } 
                        } 
                }

                random_idx = random_functions::nextInt(0, eligibles.size()-1);
                start_node = eligibles[random_idx]; 

                for( unsigned i = 0; i < lhs_boundary.size(); i++) {
                        NodeID node = lhs_boundary[i];
                        if( node == start_node ) {
                                random_idx = i;
                        }
                }

                std::swap(lhs_boundary[random_idx], lhs_boundary[max_idx-1]); lhs_boundary.pop_back();
                max_idx--;
        } while(!m_eligible[start_node] && max_idx != 0);
        return max_gain;
}

inline
void augmented_Qgraph_fabric::rebalance_fall_back(PartitionConfig & config, 
                                                  graph_access  & G, 
                                                  graph_access & G_bar, 
                                                  complete_boundary & boundary, 
                                                  std::vector< NodeID > & candidates, 
                                                  std::vector< int > & parent,
                                                  augmented_Qgraph & aqg) {

        std::vector<bool> eligible_(G.number_of_nodes(), true); 
        random_functions::permutate_vector_good_small(candidates);

        std::vector<simple_move> best_path;
        Gain best_path_gain = std::numeric_limits< Gain >::min();
        for( unsigned i = 0; i < candidates.size(); i++) {
                //for each start_vertice compute a path to a root and perform moves that can depend on each other
                int cur_block      = candidates[i];
                Gain cur_path_gain = 0;

                std::vector<simple_move> cur_path;
                if(parent[cur_block] != -1) {
                        std::vector<PartitionID> block_path;
                        while( boundary.getBlockWeight( cur_block ) <= config.upper_bound_partition ) {
                                block_path.push_back(cur_block);
                                cur_block = parent[cur_block];
                        }
                        block_path.push_back(cur_block);

                        std::reverse(block_path.begin(), block_path.end());

                        //push down from high to low
                        for( unsigned i = 0; i < block_path.size()-1; i++) {
                                PartitionID lhs = block_path[i];
                                PartitionID rhs = block_path[i+1];

                                PartialBoundary & lhs_b = boundary.getDirectedBoundary(lhs, lhs, rhs);

                                std::vector<NodeID> lhs_boundary; // todo list
                                forall_boundary_nodes(lhs_b, node) {
                                        lhs_boundary.push_back(node);
                                } endfor

                                NodeID node;
                                cur_path_gain += find_eligible_start_node( G, lhs, rhs, lhs_boundary, eligible_, node);

                                simple_move sm;
                                sm.from = lhs;
                                sm.to   = rhs;
                                sm.node = node;

                                cur_path.push_back(sm);
                                perform_simple_move(config, G, boundary, node, lhs, rhs);
                        }

               }

               //undo these changes
               for( unsigned i = 0; i < cur_path.size(); i++) {
                       perform_simple_move(config, G, boundary, cur_path[i].node, cur_path[i].to, cur_path[i].from);
               }

               if(cur_path_gain > best_path_gain) {
                       best_path = cur_path;
                       best_path_gain = cur_path_gain;
               }
 

        }

        //apply the best found movements
        for( unsigned i = 0; i < best_path.size(); i++) {
                perform_simple_move(config, G, boundary, best_path[i].node, best_path[i].from, best_path[i].to);
        }
}

inline
void augmented_Qgraph_fabric::perform_simple_move( PartitionConfig & config, 
                                                    graph_access & G, 
                                                    complete_boundary & boundary, 
                                                    NodeID & node,
                                                    PartitionID & from, 
                                                    PartitionID & to) {
        boundary_pair pair;
        pair.k   = config.k;
        pair.lhs = from;
        pair.rhs = to;

        // perform the move
        G.setPartitionIndex(node, to);

        boundary.postMovedBoundaryNodeUpdates(node, &pair, true, true);

        NodeWeight this_nodes_weight = G.getNodeWeight(node);
        boundary.setBlockNoNodes(from, boundary.getBlockNoNodes(from)-1);
        boundary.setBlockNoNodes(to,   boundary.getBlockNoNodes(to)+1);
        boundary.setBlockWeight( from, boundary.getBlockWeight(from)-this_nodes_weight);
        boundary.setBlockWeight( to,   boundary.getBlockWeight(to)+this_nodes_weight);

}

inline 
void augmented_Qgraph_fabric::move_node(PartitionConfig & config, 
                                        graph_access &  G, 
                                        NodeID & node, 
                                        refinement_pq * queue, 
                                        refinement_pq * to_queue, 
                                        complete_boundary & boundary, 
                                        PartitionID & from, 
                                        PartitionID & to) {

        G.setPartitionIndex(node, to);        
        m_eligible[node] = false;

        boundary_pair pair;
        pair.k   = config.k;
        pair.lhs = from;
        pair.rhs = to;

        boundary.postMovedBoundaryNodeUpdates(node, &pair, true, true);

        NodeWeight this_nodes_weight = G.getNodeWeight(node);
        boundary.setBlockNoNodes(from, boundary.getBlockNoNodes(from)-1);
        boundary.setBlockNoNodes(to,   boundary.getBlockNoNodes(to)+1);
        boundary.setBlockWeight( from, boundary.getBlockWeight(from)-this_nodes_weight);
        boundary.setBlockWeight( to,   boundary.getBlockWeight(to)+this_nodes_weight);

        ASSERT_TRUE(boundary.assert_bnodes_in_boundaries());
        ASSERT_TRUE(boundary.assert_boundaries_are_bnodes());

        //update gain of neighbors / the boundaries have allready been updated
        refinement_pq* cur_queue = NULL;

        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                if(G.getPartitionIndex(target) == from) {
                        cur_queue = queue;         
                } else if (G.getPartitionIndex(target) == to) {
                        cur_queue = to_queue;
                } else { continue; }
        
                EdgeWeight int_degree = 0;
                EdgeWeight ext_degree = 0;

                PartitionID target_pid = G.getPartitionIndex(target);
                PartitionID other_pid  = target_pid == from ? to : from;

                m_twfm.int_ext_degree(G, target, target_pid, other_pid, int_degree, ext_degree); 
                Gain gain = ext_degree - int_degree;

                if(cur_queue->contains(target)) {
                        if(ext_degree > 0) {
                                cur_queue->changeKey(target, gain);
                        } else {
                                cur_queue->deleteNode(target);
                        }
                } else {
                        if(ext_degree > 0) {
                                if(m_eligible[target]) {
                                        cur_queue->insert(target, gain);
                                } 
                        } 
                }

        } endfor
}

inline 
void augmented_Qgraph_fabric::move_node(PartitionConfig & config, 
                                        graph_access &  G, 
                                        NodeID & node, 
                                        refinement_pq * queue, 
                                        complete_boundary & boundary, 
                                        PartitionID & from, 
                                        PartitionID & to) {

        G.setPartitionIndex(node, to);        
        m_eligible[node] = false;

        boundary_pair pair;
        pair.k   = config.k;
        pair.lhs = from;
        pair.rhs = to;

        boundary.postMovedBoundaryNodeUpdates(node, &pair, true, true);

        NodeWeight this_nodes_weight = G.getNodeWeight(node);
        boundary.setBlockNoNodes(from, boundary.getBlockNoNodes(from)-1);
        boundary.setBlockNoNodes(to,   boundary.getBlockNoNodes(to)+1);
        boundary.setBlockWeight( from, boundary.getBlockWeight(from)-this_nodes_weight);
        boundary.setBlockWeight( to,   boundary.getBlockWeight(to)+this_nodes_weight);

        ASSERT_TRUE(boundary.assert_bnodes_in_boundaries());
        ASSERT_TRUE(boundary.assert_boundaries_are_bnodes());

        //update gain of neighbors / the boundaries have allready been updated
        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                if(G.getPartitionIndex(target) != from) continue;
        
                EdgeWeight int_degree = 0;
                EdgeWeight ext_degree = 0;

                m_twfm.int_ext_degree(G, target, from, to, int_degree, ext_degree); 
                Gain gain = ext_degree - int_degree;

                if(queue->contains(target)) {
                        if(ext_degree > 0) {
                                queue->changeKey(target, gain);
                        } else {
                                queue->deleteNode(target);
                        }
                } else {
                        if(ext_degree > 0) {
                                if(m_eligible[target]) {
                                        queue->insert(target, gain);
                                } 
                        } 
                }

        } endfor
}


inline
bool augmented_Qgraph_fabric::local_search(PartitionConfig & config, 
                                        bool plus,
                                        graph_access & G, 
                                        complete_boundary & boundary,
                                        augmented_Qgraph & aqg,
                                        boundary_pair & pair,
                                        unsigned s ) {
        return construct_local_searches_on_qgraph_edge( config, G, boundary, aqg, pair, s, plus);
}

#endif /* end of include guard: AUGMENTED_QGRAPH_FABRIC_MULTITRY_FM_PVGY97EW*/


