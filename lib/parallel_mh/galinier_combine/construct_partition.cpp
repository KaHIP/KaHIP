/******************************************************************************
 * construct_partition.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "construct_partition.h"
#include "graph_partitioner.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "random_functions.h"
#include "uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "uncoarsening/refinement/mixed_refinement.h"
#include "uncoarsening/refinement/refinement.h"
#include "uncoarsening/refinement/tabu_search/tabu_search.h"


construct_partition::construct_partition() {
                
}

construct_partition::~construct_partition() {
                
}

void construct_partition::construct_starting_from_partition( PartitionConfig & config, graph_access & G) {
	std::vector< std::queue< NodeID > >     queues(config.k);
	std::vector< NodeID >                   unassigned_vertices;
	std::vector< std::vector< NodeID > >    blocks(config.k);

	forall_nodes(G, node) {
		if( G.getPartitionIndex(node) == config.k ) {
			unassigned_vertices.push_back(node);
		} else {
			blocks[G.getPartitionIndex(node)].push_back(node);
		}
	} endfor


	//shuffle the blocks
	std::vector< NodeWeight > block_weights(config.k, 0);
	for( unsigned block = 0; block < config.k; block++) {
		random_functions::permutate_vector_good(blocks[block], false);
		block_weights[block] = blocks[block].size();
	}

	for( unsigned block = 0; block < config.k; block++) {
		for( unsigned j = 0; j < blocks[block].size(); j++) {
			NodeID node = blocks[block][j];
			forall_out_edges(G, e, node) {
				NodeID target = G.getEdgeTarget(e);
				if(G.getPartitionIndex(target) == config.k) {//unassigned
					queues[block].push(target);
				}
			} endfor
		}
	}

	unsigned no_unassigned = unassigned_vertices.size();
	while(no_unassigned > 0) {
		PartitionID block = 0;
		NodeWeight min_weight = G.number_of_nodes();
		for( unsigned cur_block = 0; cur_block < config.k; cur_block++) {
			if( block_weights[cur_block] < min_weight) {
				block = cur_block;
				min_weight = block_weights[cur_block];
			} 
		}
		
		if( queues[block].size() != 0) {
			NodeID front_node = queues[block].front();
			queues[block].pop();
			if(G.getPartitionIndex(front_node) == config.k) { //safe to assign
				G.setPartitionIndex(front_node, block);
				block_weights[block]++;
				no_unassigned--;
				forall_out_edges(G, e, front_node) {
					NodeID target = G.getEdgeTarget(e);
					if(G.getPartitionIndex(target) == config.k) {
						queues[block].push(target);
					}
				} endfor
			}
		} else {
			if( unassigned_vertices.size() > 0) {
				NodeID node = unassigned_vertices[0];
				do {
					unsigned idx = random_functions::nextInt(0, unassigned_vertices.size()-1);
					node = unassigned_vertices[idx];
					if( G.getPartitionIndex(node) != config.k ) {
						std::swap(unassigned_vertices[idx], 
                                                          unassigned_vertices[unassigned_vertices.size() -1]);	

						unassigned_vertices.pop_back();
					} else {
						queues[block].push(node);
						break;
					}
				} while( unassigned_vertices.size() != 0 );
			}
		}

	
	}
}

void construct_partition::createIndividuum( PartitionConfig & config, graph_access & G, Individuum & ind, bool output) {
        std::cout <<  "creating individuum "  << std::endl;
        forall_nodes(G, node) {
                G.setPartitionIndex(node, config.k);
        } endfor

        for( unsigned block = 0; block < config.k; block++) {
                NodeID node = 0;
                do {
                        node = random_functions::nextInt(0, G.number_of_nodes() - 1);
                } while( G.getPartitionIndex(node) != config.k );

                G.setPartitionIndex(node, block);
        }

        construct_starting_from_partition( config, G);

	complete_boundary boundary(&G);
	boundary.build();

	tabu_search ts;
	PartitionConfig copy = config;
	copy.maxIter         = G.number_of_nodes();

	ts.perform_refinement( copy, G, boundary);

	int* partition_map = new int[G.number_of_nodes()];
	forall_nodes(G, node) {
		partition_map[node] = G.getPartitionIndex(node);
	} endfor

        quality_metrics qm; 
        ind.objective     = qm.objective(config, G, partition_map);
        ind.partition_map = partition_map;
        ind.cut_edges     = new std::vector<EdgeID>();

        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if(partition_map[node] != partition_map[target]) {
                                ind.cut_edges->push_back(e);
                        }
                } endfor
        } endfor
}
