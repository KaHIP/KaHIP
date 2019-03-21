//
// Author: Christian Schulz <christian.schulz.phone@gmail.com>
// 

#include <algorithm>
#include "fm_ns_local_search.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "tools/random_functions.h"

fm_ns_local_search::fm_ns_local_search() {
                
}

fm_ns_local_search::~fm_ns_local_search() {
                
}

EdgeWeight fm_ns_local_search::perform_refinement(const PartitionConfig & config, graph_access & G, bool balance, PartitionID to) {

        std::vector< maxNodeHeap > queues; queues.resize(2);
        std::vector< bool > moved_out_of_separator(G.number_of_nodes(), false);
        std::vector< change_set > rollback_info;

        std::vector< NodeID > start_nodes;
        forall_nodes(G, node) {
                if( G.getPartitionIndex(node) == 2 ) {
                        start_nodes.push_back(node);
                }
        } endfor
        if(start_nodes.empty()) return 0;

        random_functions::permutate_vector_good(start_nodes, false);
        
        for( NodeID node : start_nodes ) {
                Gain toLHS = 0;
                Gain toRHS = 0;
                compute_gain( G, node, toLHS, toRHS);

                queues[0].insert(node, toLHS);
                queues[1].insert(node, toRHS);
        }
        
        std::vector< NodeWeight > block_weights(3,0);
        forall_nodes(G, node) {
                if( G.getPartitionIndex(node) == 0) {
                        block_weights[0] += G.getNodeWeight(node);
                } else if( G.getPartitionIndex(node) == 1 ) {
                        block_weights[1] += G.getNodeWeight(node);
                } else {
                        block_weights[2] += G.getNodeWeight(node);
                }
        } endfor

        NodeWeight best_separator  = block_weights[2];
        NodeWeight input_separator = block_weights[2];
        int best_diff              = abs((int)block_weights[1]-(int)block_weights[0]);
        int undo_idx = 0;

        int steps_till_last_improvement = 0;
        //roll forwards
        while( steps_till_last_improvement < config.sep_fm_unsucc_steps) {
                Gain gainToA = queues[0].maxValue();
                Gain gainToB = queues[1].maxValue();

                Gain top_gain        = 0;
                PartitionID to_block = 0;
               
                if(balance) {
                        top_gain = queues[to].maxValue();
                        to_block = to;
                } else {
                        if( gainToA == gainToB ) {
                                top_gain = gainToA;
                                to_block = random_functions::nextInt(0,1); 
                        } else {
                                top_gain = gainToA > gainToB ? gainToA : gainToB;
                                to_block = top_gain == gainToA ? 0 : 1;
                        }
                }

                Gain other_gain = gainToA > gainToB ? gainToB : gainToA;
                PartitionID other_block = to_block == 0 ? 1 : 0;

                NodeID nodeToBlock = queues[to_block].maxElement();
                if( block_weights[to_block] + G.getNodeWeight(nodeToBlock) < config.upper_bound_partition ) {
                        queues[to_block].deleteMax();
                        queues[other_block].deleteNode(nodeToBlock);
                        move_node(G, nodeToBlock, to_block, other_block, block_weights, moved_out_of_separator, queues, rollback_info);
                } else {
                        NodeID nodeOtherBlock = queues[other_block].maxElement();
                        if( other_gain >= 0 && block_weights[other_block] + G.getNodeWeight(nodeOtherBlock) < config.upper_bound_partition) {
                                queues[other_block].deleteMax();
                                queues[to_block].deleteNode(nodeOtherBlock);
                                move_node(G, nodeOtherBlock, other_block, to_block, block_weights, moved_out_of_separator, queues, rollback_info);
                        } else {
                                // need to make progress (remove a random node from the queues)
                                if( nodeOtherBlock == nodeToBlock ) {
                                        queues[0].deleteMax();
                                        queues[1].deleteMax();
                                } else {
                                        int block = random_functions::nextInt(0,1);
                                        queues[block].deleteMax();
                                }
                        }
                }

                int cur_diff = abs((int)block_weights[1]-(int)block_weights[0]);
                if( block_weights[2] < best_separator || (block_weights[2] == best_separator && cur_diff < best_diff)  ) {
                        best_separator = block_weights[2];
                        undo_idx = rollback_info.size();
                        steps_till_last_improvement = 0;
                }  else {
                        steps_till_last_improvement++;
                }

                if( queues[0].empty() || queues[1].empty() ) {
                        break;
                }
        }

        // roll back 
        for( int i = rollback_info.size()-1; i >= undo_idx; i--) {
                G.setPartitionIndex(rollback_info[i].node, rollback_info[i].block);
        }
                  
        return input_separator - best_separator;

}

EdgeWeight fm_ns_local_search::perform_refinement(const PartitionConfig & config, graph_access & G, 
                                                  std::vector< NodeWeight > & block_weights, 
                                                  std::vector< bool > & moved_out_of_separator,
                                                  PartialBoundary & separator, bool balance, PartitionID to) {

        std::vector< maxNodeHeap > queues; queues.resize(2);
        std::vector< change_set > rollback_info;

        std::vector< NodeID > start_nodes;
        forall_boundary_nodes( separator, node ) {
                 start_nodes.push_back(node);
        } endfor 

        random_functions::permutate_vector_good(start_nodes, false);
        
        for( NodeID node : start_nodes ) {
                Gain toLHS = 0;
                Gain toRHS = 0;
                compute_gain( G, node, toLHS, toRHS);

                queues[0].insert(node, toLHS);
                queues[1].insert(node, toRHS);
        }
        
        std::vector< NodeWeight > best_block_weights(3,0);
        best_block_weights = block_weights;
        NodeWeight best_separator  = block_weights[2];
        NodeWeight input_separator = block_weights[2];
        int best_diff              = abs((int)block_weights[1]-(int)block_weights[0]);
        int undo_idx = 0;

        int steps_till_last_improvement = 0;
        //roll forwards
        while( steps_till_last_improvement < config.sep_fm_unsucc_steps) {
                Gain gainToA = queues[0].maxValue();
                Gain gainToB = queues[1].maxValue();

                Gain top_gain        = 0;
                PartitionID to_block = 0;
               
                if(balance) {
                        top_gain = queues[to].maxValue();
                        to_block = to;
                } else {
                        if( gainToA == gainToB ) {
                                top_gain = gainToA;
                                to_block = random_functions::nextInt(0,1); 
                        } else {
                                top_gain = gainToA > gainToB ? gainToA : gainToB;
                                to_block = top_gain == gainToA ? 0 : 1;
                        }
                }
                
                Gain other_gain = gainToA > gainToB ? gainToB : gainToA;
                PartitionID other_block = to_block == 0 ? 1 : 0;

                NodeID nodeToBlock = queues[to_block].maxElement();
                if( block_weights[to_block] + G.getNodeWeight(nodeToBlock) < config.upper_bound_partition ) {
                        queues[to_block].deleteMax();
                        queues[other_block].deleteNode(nodeToBlock);
                        move_node(G, nodeToBlock, to_block, other_block, block_weights, moved_out_of_separator, queues, rollback_info, separator);
                } else {
                        NodeID nodeOtherBlock = queues[other_block].maxElement();
                        if( other_gain >= 0 && block_weights[other_block] + G.getNodeWeight(nodeOtherBlock) < config.upper_bound_partition) {
                                queues[other_block].deleteMax();
                                queues[to_block].deleteNode(nodeOtherBlock);
                                move_node(G, nodeOtherBlock, other_block, to_block, block_weights, moved_out_of_separator, queues, rollback_info, separator);
                        } else {
                                // need to make progress (remove a random node from the queues)
                                if( nodeOtherBlock == nodeToBlock ) {
                                        queues[0].deleteMax();
                                        queues[1].deleteMax();
                                } else {
                                        int block = random_functions::nextInt(0,1);
                                        queues[block].deleteMax();
                                }
                        }
                }

                int cur_diff = abs((int)block_weights[1]-(int)block_weights[0]);
                if( block_weights[2] < best_separator || (block_weights[2] == best_separator && cur_diff < best_diff)  ) {
                        best_separator              = block_weights[2];
                        undo_idx                    = rollback_info.size();
                        steps_till_last_improvement = 0;
                        best_block_weights          = block_weights;
                }  else {
                        steps_till_last_improvement++;
                }

                if( queues[0].empty() || queues[1].empty() ) {
                        break;
                }
        }

        // roll back 
        for( int i = rollback_info.size()-1; i >= undo_idx; i--) {
                if( G.getPartitionIndex( rollback_info[i].node ) == 2 ) separator.deleteNode( rollback_info[i].node );
                G.setPartitionIndex(rollback_info[i].node, rollback_info[i].block);
                if( G.getPartitionIndex( rollback_info[i].node ) == 2 ) separator.insert( rollback_info[i].node );
        }
        block_weights = best_block_weights;
                  
        for( NodeID node : moved_nodes ) {
                moved_out_of_separator[node] = false;
        }
        moved_nodes.clear();

        return input_separator - best_separator;
}

