/******************************************************************************
 * parallel_label_compress.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PARALLEL_LABEL_COMPRESS_9ME4H8DK
#define PARALLEL_LABEL_COMPRESS_9ME4H8DK

#include <math.h>

#include "data_structure/parallel_graph_access.h"
#include "data_structure/processor_tree.h"
#include "ppartition_config.h"
#include "tools/random_functions.h"
#include "hmap_wrapper.h"
#include "node_ordering.h"


template <typename T> 
class parallel_label_compress {
        public:
                parallel_label_compress() {};
                virtual ~parallel_label_compress() {};

                void perform_parallel_label_compression( PPartitionConfig & config, 
                                parallel_graph_access & G, bool balance, bool for_coarsening = true,
                                processor_tree PEtree=processor_tree() ) {

                        if( config.label_iterations == 0) return;
                        NodeWeight cluster_upperbound = config.upper_bound_cluster;

//TODO: when coarsening, do we need to check all nodes? maybe we can limit to only boundary nodes
                        std::vector< NodeID > permutation( G.number_of_local_nodes() );
                        if( for_coarsening ) {
                                node_ordering no; 
                                no.order_nodes( config, G, permutation); 
                        } else {
                                random_functions::permutate_vector_fast( permutation, true);
                        }

                        //use distance if integrated mapping is activated and we do uncoarsening
                        bool usePEdistances = !for_coarsening && config.integrated_mapping ? true : false ;

                        //checks, keep?
                        if( usePEdistances ){
                                //tree should not be empty
                                assert( PEtree.get_numPUs()>1 );
                        }
                        if( PEtree.get_numPUs()==1 ){
                                //if tree is empty, PE distances should not be used
                                assert( !usePEdistances );
                        }

                        int clz = __builtin_clzll(G.number_of_global_nodes()); // index of highest bit
                        const int label_size = 8*sizeof(unsigned long long int) - clz;

                        std::cout << __FILE__ <<": "<< G.number_of_global_nodes() << " bit label size = " << label_size
                                << ", usePEdistances " << usePEdistances << std::endl;
                        
                        //std::unordered_map<NodeID, NodeWeight> hash_map;
                        hmap_wrapper< T > hash_map(config);
                        hash_map.init( G.get_max_degree() );
                        for( ULONG i = 0; i < config.label_iterations; i++) {
                                NodeID prev_node = 0;
                                forall_local_nodes(G, rnode) {
                                        const NodeID node = permutation[rnode]; // use the current random node

                                        //move the node to the cluster that is most common in the neighborhood
                                        //second sweep for finding max and resetting array
                                        const PartitionID old_block   = G.getNodeLabel(node);
                                        const NodeWeight  node_weight = G.getNodeWeight(node);
                                        PartitionID max_value   = 0;
                                        PartitionID max_block   = G.getNodeLabel(node);

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
                                                forall_out_edges(G, e, node) {
                                                        const NodeID target             = G.getEdgeTarget(e);
                                                        const PartitionID cur_block     = G.getNodeLabel(target);
                                                        //PartitionID cur_value;
                                                        long long cur_value;

                                                        if( usePEdistances ){
//TODO: verify that these are the correct variables
//TODO: take minimum if PU distances are used or calculate negative
//TODO: when coarsening, which are the blocks? is there a partition?
//TODO: change when for_coarsening
//TODO: does the hash map accepts negative values?
                                                                cur_value += (-1) * G.getEdgeWeight(e) * PEtree.getDistance_PxPy( old_block, cur_block) ;
printf("%lld",cur_value);
                                                        }else{
                                                                cur_value += G.getEdgeWeight(e);
                                                        }
                                                        hash_map[cur_block] = cur_value;

                                                        bool improvement = cur_value > max_value;
                                                        improvement |= cur_value == max_value && random_functions::nextBool();

                                                        bool sizeconstraint = G.getBlockSize(cur_block) + node_weight <= cluster_upperbound;
                                                        sizeconstraint |= cur_block == old_block;

                                                        bool cycle = !config.vcycle;
                                                        cycle |= G.getSecondPartitionIndex( node ) == G.getSecondPartitionIndex(target);

                                                        bool balancing = own_block_balanced || cur_block != old_block;
                                                        if( improvement  && sizeconstraint && cycle && balancing) {
                                                                max_value = cur_value;
                                                                max_block = cur_block;
                                                        }
                                                } endfor
                                        }

                                        if( old_block != max_block ) {
                                                G.setNodeLabel(node, max_block);

                                                G.setBlockSize(old_block, G.getBlockSize(old_block) - node_weight);
                                                G.setBlockSize(max_block, G.getBlockSize(max_block) + node_weight);
                                        }

                                        prev_node = node;
                                        G.update_ghost_node_data(); 
                                        hash_map.clear();

                                } endfor
                                G.update_ghost_node_data_finish(); 
                        }
                }

};


#endif /* end of include guard: PARALLEL_LABEL_COMPRESS_9ME4H8DK */
