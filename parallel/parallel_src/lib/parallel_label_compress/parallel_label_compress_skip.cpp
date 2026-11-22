/******************************************************************************
 * parallel_label_compress_skip.cpp
 * Refinement LP that skips provably clean nodes (see sweeps_skip). In its own translation unit, so
 * that the code of the plain LP (parallel_label_compress.h, compiled in distributed_partitioner.cpp)
 * stays as it is without the skip.
 *****************************************************************************/

#include "parallel_label_compress/parallel_label_compress.h"

// by_work: rounds of equal work (refinement), otherwise of equal node counts (coarsening)
// skip: skip clean nodes (refinement): a node is clean if it has no ghost neighbour, its last scan in
// this call saw only its own label and drew no random number, and no local neighbour moved since.
// Its scan would return stay without a random draw (every candidate is the own block, so no size or
// balance test applies, and the value grows with every edge), so it is skipped; its round work counts.
template < bool by_work, typename T >
PARHIP_ALIGN_CODE_64
static void sweeps_skip( PPartitionConfig & config, parallel_graph_access & G, bool balance,
         NodeWeight cluster_upperbound, std::vector< NodeID > & permutation,
         hmap_wrapper< T > & hash_map, std::vector< unsigned char > & m_clean ) {
        const bool skip = true;
        for( ULONG i = 0; i < config.label_iterations; i++) {
                NodeID prev_node = 0;
                forall_local_nodes(G, rnode) {
                        NodeID node = permutation[rnode]; // use the current random node
                        // the nodes come in permutation order, so the next nodes' records, first edges,
                        // data and first target labels miss; hint them a few nodes ahead (no effect on results)
                        if( rnode + 4 < end ) G.prefetchNodeRecord(permutation[rnode + 4]);
                        if( rnode + 2 < end ) G.prefetchEdge(G.get_first_edge(permutation[rnode + 2]));
                        if( rnode + 1 < end ) {
                                NodeID next = permutation[rnode + 1];
                                G.prefetchNodeData(next);
                                EdgeID b = G.get_first_edge(next), b_end = G.get_first_invalid_edge(next);
                                if( b_end > b + 8 ) b_end = b + 8;
                                for( EdgeID e = b; e < b_end; e++) G.prefetchTargetData(e);
                        }

                        if( skip && m_clean[node] ) {
                                prev_node = node;
                                if( by_work ) G.update_ghost_node_data_by_work( G.getNodeDegree(node) + 1 );
                                else G.update_ghost_node_data();
                                continue;
                        }

                        //move the node to the cluster that is most common in the neighborhood
                        //second sweep for finding max and resetting array
                        PartitionID max_block   = G.getNodeLabel(node);
                        PartitionID old_block   = G.getNodeLabel(node);
                        PartitionID max_value   = 0;
                        bool        mixed       = false; // skip: a foreign label or a random draw seen
                        NodeWeight  node_weight = G.getNodeWeight(node);
                        bool own_block_balanced = G.getBlockSize(old_block) <= cluster_upperbound || !balance;

                        if( G.getNodeDegree(node) == 0) {
                                // find a block to assign it to
                                if(config.vcycle) {
                                        NodeWeight prev_block_size = G.getBlockSize( G.getNodeLabel( prev_node ) );
                                        bool same_block = G.getSecondPartitionIndex(prev_node)==G.getSecondPartitionIndex(node);
                                        if( prev_block_size  + node_weight <= cluster_upperbound && same_block ) {
                                                max_block = G.getNodeLabel( prev_node );
                                        }
                                } else {
                                        NodeWeight prev_block_size = G.getBlockSize( G.getNodeLabel( prev_node ) );
                                        if( prev_block_size  + node_weight <= cluster_upperbound) {
                                                max_block = G.getNodeLabel( prev_node );
                                        }
                                }

                        } else {
                                hash_map.prepare( G.getNodeDegree(node) );
                                forall_out_edges(G, e, node) {
                                        if( e + 8 < end ) G.prefetchTargetData(e + 8);
                                        NodeID target             = G.getEdgeTarget(e);
                                        PartitionID cur_block     = G.getNodeLabel(target);
                                        NodeWeight & cur_weight   = hash_map[cur_block]; // one lookup per edge
                                        cur_weight               += G.getEdgeWeight(e);
                                        PartitionID cur_value     = cur_weight;

                                        if( skip ) mixed |= (cur_block != old_block) | (cur_value == max_value);
                                        bool improvement = cur_value > max_value;
                                        improvement |= cur_value == max_value && random_functions::nextBool();

                                        if( !improvement ) continue; // nothing below has side effects

                                        bool cycle = !config.vcycle;
                                        cycle |= G.getSecondPartitionIndex( node ) == G.getSecondPartitionIndex(target);

                                        bool balancing = own_block_balanced || cur_block != old_block;
                                        if( !(cycle && balancing) ) continue;

                                        bool sizeconstraint = cur_block == old_block;
                                        sizeconstraint = sizeconstraint || G.getBlockSize(cur_block) + node_weight <= cluster_upperbound;
                                        if( sizeconstraint ) {
                                                max_value = cur_value;
                                                max_block = cur_block;
                                        }
                                } endfor
                        }

                        if( old_block != max_block ) {
                                G.setNodeLabel(node, max_block);

                                G.setBlockSize(old_block, G.getBlockSize(old_block) - node_weight);
                                G.setBlockSize(max_block, G.getBlockSize(max_block) + node_weight);
                                if( skip ) { // its local neighbours see a new label
                                        const NodeID num_local = G.number_of_local_nodes();
                                        forall_out_edges(G, e, node) {
                                                NodeID target = G.getEdgeTarget(e);
                                                if( target < num_local ) m_clean[target] = 0;
                                        } endfor
                                }
                        }
                        if( skip ) m_clean[node] = !mixed && G.getNodeDegree(node) > 0 && old_block == max_block && !G.is_interface_node(node);

                        prev_node = node;
                        if( by_work ) G.update_ghost_node_data_by_work( G.getNodeDegree(node) + 1 ); 
                        else G.update_ghost_node_data(); 
                        hash_map.clear();

                } endfor
                G.update_ghost_node_data_finish(); 
        }
}

template <typename T>
bool parallel_label_compress< T >::skipping_sweeps( PPartitionConfig & config, parallel_graph_access & G, bool balance,
                                                    NodeWeight cluster_upperbound, std::vector< NodeID > & permutation,
                                                    hmap_wrapper< T > & hash_map, int size ) {
        const NodeID n = G.number_of_local_nodes();
        if( n == 0 ) return false;
        if( size > 1 ) {
                NodeID non_interface = 0;
                for( NodeID v = 0; v < n; v++) non_interface += !G.is_interface_node(v);
                if( 2 * (ULONG)non_interface < (ULONG)n ) return false;
        }
        std::vector< unsigned char > clean( n, 0 ); // clean flag per local node, released with the call
        if( size == 1 ) sweeps_skip< false >( config, G, balance, cluster_upperbound, permutation, hash_map, clean );
        else sweeps_skip< true >( config, G, balance, cluster_upperbound, permutation, hash_map, clean );
        return true;
}

// only this member: instantiating the class here would emit a second copy of the header's members
template bool parallel_label_compress< std::vector< NodeWeight > >::skipping_sweeps( PPartitionConfig &, parallel_graph_access &, bool,
        NodeWeight, std::vector< NodeID > &, hmap_wrapper< std::vector< NodeWeight > > &, int );
template bool parallel_label_compress< linear_probing_hashmap >::skipping_sweeps( PPartitionConfig &, parallel_graph_access &, bool,
        NodeWeight, std::vector< NodeID > &, hmap_wrapper< linear_probing_hashmap > &, int );
