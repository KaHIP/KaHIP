//
// Author: Christian Schulz <christian.schulz.phone@gmail.com>
// 

#ifndef FM_NS_LOCAL_SEARCH_P621XWW8
#define FM_NS_LOCAL_SEARCH_P621XWW8

#include "uncoarsening/refinement/quotient_graph_refinement/partial_boundary.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "definitions.h"
#include "partition_config.h"
#include "data_structure/graph_access.h"
#include "uncoarsening/refinement/quotient_graph_refinement/partial_boundary.h"

// information to perform undo 
struct change_set {
        NodeID node;
        PartitionID block;
};

class fm_ns_local_search {
public:
        fm_ns_local_search();
        virtual ~fm_ns_local_search();

        EdgeWeight perform_refinement(const PartitionConfig & config, graph_access & G, bool balance = false, PartitionID to = 4);
        EdgeWeight perform_refinement(const PartitionConfig & config, graph_access & G, std::vector< NodeWeight > & block_weight, 
                                      std::vector< bool > & moved_out_of_separator,
                                      PartialBoundary & separator, bool balance = false, PartitionID to = 4);

private: 
        void compute_gain( graph_access & G, NodeID node, Gain & toLHS, Gain & toRHS);
        void move_node( graph_access & G,  NodeID & node, PartitionID & to_block, PartitionID & other_block, 
                        std::vector< NodeWeight > & block_weights,
                        std::vector< bool > & moved_out_of_S,
                        std::vector< maxNodeHeap > & heaps,
                        std::vector< change_set > & rollback_info);

        void move_node( graph_access & G,  NodeID & node, PartitionID & to_block, PartitionID & other_block, 
                        std::vector< NodeWeight > & block_weights,
                        std::vector< bool > & moved_out_of_S,
                        std::vector< maxNodeHeap > & heaps,
                        std::vector< change_set > & rollback_info,
                        PartialBoundary & separator);

        std::vector< NodeID > moved_nodes;
};


inline
void fm_ns_local_search::compute_gain( graph_access & G, NodeID node, Gain & toLHS, Gain & toRHS) {
        toLHS = G.getNodeWeight(node);
        toRHS = G.getNodeWeight(node);

        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                if( G.getPartitionIndex(target) == 0) {
                        toRHS -= G.getNodeWeight(target);
                } else if( G.getPartitionIndex(target) == 1 ) {
                        toLHS -= G.getNodeWeight(target);
                }
        } endfor

}

inline
void fm_ns_local_search::move_node( graph_access & G, NodeID & node, PartitionID & to_block, PartitionID & other_block, 
                                    std::vector< NodeWeight > & block_weights,
                                    std::vector< bool > & moved_out_of_S, 
                                    std::vector< maxNodeHeap > & queues,
                                    std::vector< change_set > & rollback_info) {

        change_set cur_move;
        cur_move.node = node;
        cur_move.block = G.getPartitionIndex(node);
        rollback_info.push_back(cur_move);

        G.setPartitionIndex(node, to_block);
        block_weights[to_block] += G.getNodeWeight(node);
        block_weights[2] -= G.getNodeWeight(node);
        moved_out_of_S[node] = true;

        std::vector< NodeID > to_be_added;
        std::vector< NodeID > to_be_updated; // replace by hashmap?
        Gain gain_achieved = G.getNodeWeight(node);
        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);

                if( G.getPartitionIndex( target ) == other_block ) {
                        change_set cur_move;
                        cur_move.node = target;
                        cur_move.block = G.getPartitionIndex(target);
                        rollback_info.push_back(cur_move);

                        G.setPartitionIndex(target, 2);
                        block_weights[other_block] -= G.getNodeWeight(target);
                        block_weights[2]           += G.getNodeWeight(target);
                        gain_achieved              -= G.getNodeWeight(target);

                        if( !moved_out_of_S[target] ) {
                                to_be_added.push_back(target);
                        }

                        forall_out_edges(G, e_bar, target) {
                                NodeID v = G.getEdgeTarget(e_bar);
                                if( queues[0].contains(v) ) {
                                        to_be_updated.push_back(v);
                                } 
                        } endfor
                } else if(  G.getPartitionIndex( target ) == 2 ) {
                        to_be_updated.push_back(target);
                }
        } endfor

        Gain toLHS = 0;
        Gain toRHS = 0;

        for( NodeID node : to_be_added ) {
                compute_gain( G, node, toLHS, toRHS);
                queues[0].insert(node, toLHS);
                queues[1].insert(node, toRHS);
        }

        for( NodeID node : to_be_updated) {
                compute_gain( G, node, toLHS, toRHS);
                queues[0].changeKey(node, toLHS);
                queues[1].changeKey(node, toRHS);
        }
}

inline
void fm_ns_local_search::move_node( graph_access & G, NodeID & node, PartitionID & to_block, PartitionID & other_block, 
                                    std::vector< NodeWeight > & block_weights,
                                    std::vector< bool > & moved_out_of_S, 
                                    std::vector< maxNodeHeap > & queues,
                                    std::vector< change_set > & rollback_info,
                                    PartialBoundary & separator) {

        change_set cur_move;
        cur_move.node = node;
        cur_move.block = G.getPartitionIndex(node);
        rollback_info.push_back(cur_move);
        separator.deleteNode(node);

        G.setPartitionIndex(node, to_block);
        block_weights[to_block] += G.getNodeWeight(node);
        block_weights[2] -= G.getNodeWeight(node);
        moved_out_of_S[node] = true;
        moved_nodes.push_back(node);

        std::vector< NodeID > to_be_added;
        std::vector< NodeID > to_be_updated; // replace by hashmap?
        Gain gain_achieved = G.getNodeWeight(node);
        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);

                if( G.getPartitionIndex( target ) == other_block ) {
                        change_set cur_move;
                        cur_move.node = target;
                        cur_move.block = G.getPartitionIndex(target);
                        rollback_info.push_back(cur_move);

                        G.setPartitionIndex(target, 2);
                        separator.insert(target);

                        block_weights[other_block] -= G.getNodeWeight(target);
                        block_weights[2]           += G.getNodeWeight(target);
                        gain_achieved              -= G.getNodeWeight(target);

                        if( !moved_out_of_S[target] ) {
                                to_be_added.push_back(target);
                        }

                        forall_out_edges(G, e_bar, target) {
                                NodeID v = G.getEdgeTarget(e_bar);
                                if( queues[0].contains(v) ) {
                                        to_be_updated.push_back(v);
                                } 
                        } endfor
                } else if(  G.getPartitionIndex( target ) == 2 ) {
                        to_be_updated.push_back(target);
                }
        } endfor

        Gain toLHS = 0;
        Gain toRHS = 0;

        for( NodeID node : to_be_added ) {
                compute_gain( G, node, toLHS, toRHS);
                queues[0].insert(node, toLHS);
                queues[1].insert(node, toRHS);
        }

        for( NodeID node : to_be_updated) {
                compute_gain( G, node, toLHS, toRHS);
                queues[0].changeKey(node, toLHS);
                queues[1].changeKey(node, toRHS);
        }
}

#endif /* end of include guard: FM_NS_LOCAL_SEARCH_P621XWW8 */
