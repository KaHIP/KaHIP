//
// Author: Christian Schulz <christian.schulz.phone@gmail.com>
// 

#ifndef GREEDY_NS_LOCAL_SEARCH_P9KLE4NH
#define GREEDY_NS_LOCAL_SEARCH_P9KLE4NH

#include "definitions.h"
#include "partition_config.h"
#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/maxNodeHeap.h"

class greedy_ns_local_search {
public:
        greedy_ns_local_search();
        virtual ~greedy_ns_local_search();

        EdgeWeight perform_refinement(const PartitionConfig & config, graph_access & G);

private: 
        void compute_gain( graph_access & G, NodeID node, Gain & toLHS, Gain & toRHS);
        void move_node( graph_access & G,  NodeID & node, PartitionID & to_block, PartitionID & other_block, 
                        std::vector< NodeWeight > & block_weights,
                        std::vector< bool > & moved_out_of_S,
                        std::vector< maxNodeHeap > & heaps);
};

inline
void greedy_ns_local_search::compute_gain( graph_access & G, NodeID node, Gain & toLHS, Gain & toRHS) {
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
void greedy_ns_local_search::move_node( graph_access & G, NodeID & node, PartitionID & to_block, PartitionID & other_block, 
                                    std::vector< NodeWeight > & block_weights,
                                    std::vector< bool > & moved_out_of_S, 
                                    std::vector< maxNodeHeap > & queues) {
        G.setPartitionIndex(node, to_block);
        block_weights[to_block] += G.getNodeWeight(node);
        block_weights[2] -= G.getNodeWeight(node);
        moved_out_of_S[node] = true;

        std::vector< NodeID > to_be_added;
        std::vector< NodeID > to_be_updated;
        Gain gain_achieved = G.getNodeWeight(node);
        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);

                if( G.getPartitionIndex( target ) == other_block ) {
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


#endif /* end of include guard: GREEDY_NS_LOCAL_SEARCH_P9KLE4NH */
