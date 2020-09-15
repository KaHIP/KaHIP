/******************************************************************************
 * parallel_label_compress.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PARALLEL_LABEL_COMPRESS_9ME4H8DK
#define PARALLEL_LABEL_COMPRESS_9ME4H8DK

#include <math.h>
#include <map>

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
                                 const processor_tree & PEtree = processor_tree()) {

                        if( config.label_iterations == 0) return;
                        NodeWeight cluster_upperbound = config.upper_bound_cluster;

                        std::vector< NodeID > permutation( G.number_of_local_nodes() );
                        if( for_coarsening ) {
                                //std::cout <<    ": will coarsen" << std::endl;
                                node_ordering no; 
                                no.order_nodes( config, G, permutation);
                        } else {
                                //std::cout << __LINE__ << ": will refine" << std::endl;
                                random_functions::permutate_vector_fast( permutation, true);
                        }

                        //use distance if integrated mapping is activated and we do uncoarsening
                        const bool usePEdistances = !for_coarsening && config.integrated_mapping ? true : false ;

                        //checks, keep?
                        if( usePEdistances ){
                                //tree should not be empty
                                assert( PEtree.get_numPUs()>1 );
                        }
                        if( PEtree.get_numPUs()==1 ){
                                //if tree is empty, PE distances should not be used
                                assert( !usePEdistances );
                        }

                        // std::cout << "TEST print: "<< G.number_of_global_nodes()
                        //         << ", usePEdistances " << usePEdistances << ", only_boundary " << config.only_boundary << std::endl;

                        //std::unordered_map<NodeID, NodeWeight> hash_map;
                        hmap_wrapper< T > hash_map(config);
                        hash_map.init( G.get_max_degree() );

                        for( ULONG i = 0; i < config.label_iterations; i++) {
                                ULONG numChanges = 0;
                                NodeID prev_node = 0;
                                forall_local_nodes(G, rnode) {
                                        const NodeID node = permutation[rnode]; // use the current random node

                                        //TODO: evaluate how this affects quality and running time.
                                        //when NOT coarsening and we check only boundary nodes and this node is NOT a boundary node, skip it
                                        if( !for_coarsening && config.only_boundary && !is_boundary(node,G) ){
                                            continue;
                                        }

                                        //move the node to the cluster that is most common in the neighborhood
                                        //second sweep for finding max and resetting array
                                        const PartitionID old_block   = G.getNodeLabel(node);
                                        const NodeWeight  node_weight = G.getNodeWeight(node);
                                        
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
                                                if( usePEdistances){
                                                        assert( !config.vcycle );
                                                        assert( config.only_boundary ? is_boundary(node,G):true );
                                                        max_block = refine_with_PU_distances(node, G, PEtree, cluster_upperbound);
                                                }else{
                                                        max_block = coarse_or_refine(node, G, hash_map, cluster_upperbound, own_block_balanced, config.vcycle );
                                                }
                                        }

                                        if( old_block != max_block ) {
                                                G.setNodeLabel(node, max_block);

                                                G.setBlockSize(old_block, G.getBlockSize(old_block) - node_weight);
                                                G.setBlockSize(max_block, G.getBlockSize(max_block) + node_weight);
                                                numChanges++;
                                        }

                                        prev_node = node;
                                        G.update_ghost_node_data(); 
                                        hash_map.clear();

                                } endfor //for G nodes
                                G.update_ghost_node_data_finish(); 
                                //std::cout << "in iteration round " << i << ", we moved " << numChanges << " vertices" <<std::endl;
                        }//for( ULONG i = 0; i < config.label_iterations; i++)
                }


        private:
                bool is_boundary(NodeID node, parallel_graph_access & G) {

                        PartitionID my_part = G.getNodeLabel(node);

                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                if ( G.getNodeLabel(target)!=my_part ) {
                                        return true;
                                }
                        } endfor
                        return false;
                }


                PartitionID refine_with_PU_distances(
                        const NodeID& node,
                        parallel_graph_access& G,
                        const processor_tree& PEtree,
                        const NodeWeight& cluster_upperbound
                ){

                        const PartitionID old_block   = G.getNodeLabel(node);
                        const NodeWeight  node_weight = G.getNodeWeight(node);

                        //create the map with all neighboring blocks and the sum of the edge weights
                        std::map<PartitionID, EdgeWeight> ngbr_blocks;

                        forall_out_edges(G, e, node) {
                                const NodeID target             = G.getEdgeTarget(e);
                                const PartitionID ngbr_block    = G.getNodeLabel(target);
                                ngbr_blocks[ngbr_block]        += G.getEdgeWeight(e);
                        }endfor

                        //go over the map, assign this vertex to all neighboring blocks and each time
                        //calculate the communication costs 
                        EdgeWeight min_comm_cost = std::numeric_limits<EdgeWeight>::max();
                        PartitionID best_block = old_block;

                        for( auto const& x : ngbr_blocks ){
                                //the candidate new block for the node
                                const PartitionID new_block = x.first;
                                //check: if moving to new block is gonna violate the weight bound, do not consider this move
                                // if the new block is the old, then the size constraint is not violated
                                bool sizeconstraint = new_block==old_block || G.getBlockSize(new_block) + node_weight <= cluster_upperbound;
                                if( !sizeconstraint ){
                                       continue;
                                }

                                EdgeWeight comm_cost = 0;
                                //go over the map again to calculate the communication cost for the case
                                //that the node is moved to new_block
                                for( auto const& xx : ngbr_blocks ){
                                        const PartitionID ngbr_block = xx.first;
                                        const EdgeWeight comm_vol = xx.second;
                                        //skip costs for the same block
                                        if( ngbr_block==new_block ){
                                            continue;
                                        }
                                        comm_cost += PEtree.getDistance_PxPy(new_block, ngbr_block) * comm_vol;
                                }

                                //if cost is equal the min, improve or not at random
                                bool improvement = comm_cost < min_comm_cost;
                                improvement |= comm_cost == min_comm_cost && random_functions::nextBool();

                                if( improvement ){
                                        min_comm_cost = comm_cost;
                                        best_block = new_block;
                                }
                        }

                        // if(best_block!=old_block){
                        //     std::cout<< "node " << node << " will move from block " << old_block << " to " << best_block << " with min cost " << min_comm_cost  << " num ngbr blocks = " << ngbr_blocks.size() << endl;
                        // }
                        return best_block;
                }


                PartitionID coarse_or_refine(
                        const NodeID& node,
                        parallel_graph_access& G,
                        hmap_wrapper<T>& hash_map,
                        const NodeWeight& cluster_upperbound,
                        const bool& own_block_balanced,
                        const bool& vcycle
                ){
                        PartitionID max_block           = G.getNodeLabel(node);
                        PartitionID max_value           = 0;
                        const NodeWeight node_weight    = G.getNodeWeight(node);
                        const PartitionID old_block     = G.getNodeLabel(node);

                        forall_out_edges(G, e, node) {
                                NodeID target             = G.getEdgeTarget(e);
                                PartitionID cur_block     = G.getNodeLabel(target);
                                hash_map[cur_block] += G.getEdgeWeight(e);
                                PartitionID cur_value     = hash_map[cur_block];

                                bool improvement = cur_value > max_value;
                                improvement |= cur_value == max_value && random_functions::nextBool();

                                bool sizeconstraint = G.getBlockSize(cur_block) + node_weight <= cluster_upperbound;
                                sizeconstraint |= cur_block == old_block;

                                bool cycle = !vcycle;
                                cycle |= G.getSecondPartitionIndex( node ) == G.getSecondPartitionIndex(target);

                                bool balancing = own_block_balanced || cur_block != old_block;
                                if( improvement  && sizeconstraint && cycle && balancing) {
                                        max_value = cur_value;
                                        max_block = cur_block;
                                }
                        } endfor

                        return max_block;
                }

};


#endif /* end of include guard: PARALLEL_LABEL_COMPRESS_9ME4H8DK */
