/******************************************************************************
 * parallel_label_compress.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PARALLEL_LABEL_COMPRESS_9ME4H8DK
#define PARALLEL_LABEL_COMPRESS_9ME4H8DK

#include "data_structure/parallel_graph_access.h"
#include "partition_config.h"
#include "tools/random_functions.h"

// 64-byte code alignment of a function on GCC and Clang, nothing elsewhere
#if defined(__GNUC__) || defined(__clang__)
#define PARHIP_ALIGN_CODE_64 __attribute__((aligned(64)))
#else
#define PARHIP_ALIGN_CODE_64
#endif
#include "hmap_wrapper.h"
#include "node_ordering.h"


template <typename T> 
class parallel_label_compress {
        public:
                parallel_label_compress() {};
                virtual ~parallel_label_compress() {};

                // fixed code alignment: the refinement LP was up to 33% slower at some offsets of this function
                PARHIP_ALIGN_CODE_64
                void perform_parallel_label_compression( PPartitionConfig & config, 
                                parallel_graph_access & G, bool balance, bool for_coarsening = true) {

                        if( config.label_iterations == 0) return;
                        NodeWeight cluster_upperbound = config.upper_bound_cluster;

                        std::vector< NodeID > permutation( G.number_of_local_nodes() );
                        if( for_coarsening ) {
                                node_ordering no; 
                                no.order_nodes( config, G, permutation); 
                        } else {
                                random_functions::permutate_vector_fast( permutation, true);
                        }

                        //std::unordered_map<NodeID, NodeWeight> hash_map;
                        hmap_wrapper< T > hash_map(config);
                        hash_map.init( G.get_max_degree() );
                        unsigned long long prev_moved = 0;
                        for( ULONG i = 0; i < config.label_iterations; i++) {
                                unsigned long long moved = 0;
                                NodeID prev_node = 0;
                                forall_local_nodes(G, rnode) {
                                        NodeID node = permutation[rnode]; // use the current random node

                                        //move the node to the cluster that is most common in the neighborhood
                                        //second sweep for finding max and resetting array
                                        PartitionID max_block   = G.getNodeLabel(node);
                                        PartitionID old_block   = G.getNodeLabel(node);
                                        PartitionID max_value   = 0;
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
                                                moved++;
                                                G.setNodeLabel(node, max_block);

                                                G.setBlockSize(old_block, G.getBlockSize(old_block) - node_weight);
                                                G.setBlockSize(max_block, G.getBlockSize(max_block) + node_weight);
                                        }

                                        prev_node = node;
                                        G.update_ghost_node_data(); 
                                        hash_map.clear();

                                } endfor
                                G.update_ghost_node_data_finish(); 

                                // refinement: once the global number of moves stops falling by at least
                                // 5% per iteration, the remaining moves are mostly ties flipping
                                // between blocks; stop (the same decision on every PE)
                                if( !for_coarsening ) {
                                        unsigned long long global_moved = 0;
                                        MPI_Allreduce(&moved, &global_moved, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, G.getCommunicator());
                                        if( global_moved == 0 || (i > 0 && 20 * global_moved > 19 * prev_moved) ) break;
                                        prev_moved = global_moved;
                                }
                        }
                }

};


#endif /* end of include guard: PARALLEL_LABEL_COMPRESS_9ME4H8DK */
