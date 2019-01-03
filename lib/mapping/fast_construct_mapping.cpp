/******************************************************************************
 * fast_construct_mapping.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/


#include <fstream>

#include "balance_configuration.h"
#include "graph_partitioner.h"
#include "configuration.h"
#include "partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "fast_construct_mapping.h"
#include "tools/graph_extractor.h"

fast_construct_mapping::fast_construct_mapping() {

}

fast_construct_mapping::~fast_construct_mapping() {

}

void fast_construct_mapping::construct_initial_mapping_bottomup( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank) {
        m_tmp_num_nodes = C.number_of_nodes();
        construct_initial_mapping_bottomup_internal( config, C, D, 0, perm_rank);
}

void fast_construct_mapping::construct_initial_mapping_bottomup_internal( PartitionConfig & config, graph_access & C, matrix & D, int idx,  std::vector< NodeID > & perm_rank) {

        PartitionID num_parts = C.number_of_nodes()/config.group_sizes[idx];
        partition_C_perfectly_balanced( config, C, num_parts);

        if( idx ==(int)(config.group_sizes.size() - 1) ) {
                // build initial offsets 
                int nodes_per_block = m_tmp_num_nodes / config.group_sizes[idx];
                perm_rank[0] = 0;
                for( unsigned int block = 1; block < perm_rank.size(); block++) {
                        perm_rank[block] = perm_rank[block-1]+nodes_per_block;
                }
        } else {
                //contract partitioned graph
                graph_access Q; complete_boundary bnd(&C);
                bnd.build();
                bnd.getUnderlyingQuotientGraph(Q);
               
                std::vector< NodeID > rec_ranks( num_parts, 0);
                construct_initial_mapping_bottomup_internal( config, Q, D, idx+1, rec_ranks);

                //recompute offsets 
                forall_nodes(C, node) {
                        PartitionID block = C.getPartitionIndex(node);
                        perm_rank[node]   = rec_ranks[block];
                        rec_ranks[block] += C.getNodeWeight(node);
                } endfor
        }
}

void fast_construct_mapping::construct_initial_mapping_topdown( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank) {

        std::vector< NodeID > m_mapping(C.number_of_nodes());
        forall_nodes(C, node) {
                m_mapping[node] = node;
        } endfor

        construct_initial_mapping_topdown_internal( config, C, config.group_sizes, 0, m_mapping, perm_rank);
}

void fast_construct_mapping::construct_initial_mapping_topdown_internal( PartitionConfig & config, 
                graph_access & C, 
                std::vector< int > group_sizes, 
                int start_id, 
                std::vector< NodeID > & map_to_original, 
                std::vector< NodeID > & perm_rank) {

        PartitionID num_parts = group_sizes[group_sizes.size()-1];
        if( num_parts == 1 ) {
                if( group_sizes.size() == 1 ) return;
                group_sizes.pop_back();
                return construct_initial_mapping_topdown_internal( config, C, group_sizes, 0, map_to_original, perm_rank);
        } 

        partition_C_perfectly_balanced( config, C, num_parts);

        int nodes_per_block = C.number_of_nodes() / num_parts;
        std::vector< int > count( num_parts, start_id);
        for( PartitionID block = 1; block < num_parts; block++) {
                count[block] = count[block-1]+nodes_per_block;
        }

        int stop_number = std::max(2, (int)config.group_sizes.size() - config.max_recursion_levels_construction);
        if( (int)group_sizes.size() == stop_number ) {
                forall_nodes(C, node) {
                        PartitionID block = C.getPartitionIndex(node);
                        perm_rank[map_to_original[node]] = count[block];
                        count[block]++;
                } endfor
        } else {
                // extract subgraphs and recurse on them
                group_sizes.pop_back();
                for( PartitionID block = 0; block < num_parts; block++) {
                        graph_extractor ge; graph_access Q;
                        std::vector<NodeID> mapping;
                        ge.extract_block( C, Q, block, mapping);

                        forall_nodes(Q, node) {
                                mapping[node] = map_to_original[mapping[node]];
                        } endfor

                        construct_initial_mapping_topdown_internal( config, Q, group_sizes, count[block], mapping, perm_rank);
                }
        }
}

void fast_construct_mapping::partition_C_perfectly_balanced( PartitionConfig & config, graph_access & C, PartitionID blocks) {
        std::streambuf* backup = std::cout.rdbuf();
        std::ofstream ofs;
        ofs.open("/dev/null");
        std::cout.rdbuf(ofs.rdbuf()); 

        PartitionConfig partition_config = config;
        configuration cfg; 
        switch(partition_config.preconfiguration_mapping) {
                case PRE_CONFIG_MAPPING_FAST:
                        cfg.fast( partition_config );
                        break;
                case PRE_CONFIG_MAPPING_ECO:
                        cfg.eco( partition_config );
                        break;
                case PRE_CONFIG_MAPPING_STRONG:
                        cfg.strong( partition_config );
                        break;
                default:
                        cfg.fast( partition_config );
        }

        partition_config.k = blocks;
        partition_config.imbalance = 0;
        partition_config.epsilon = 0;

        std::vector< NodeWeight > weights(C.number_of_nodes());
        forall_nodes(C, node) {
                weights[node] = C.getNodeWeight(node);
                C.setNodeWeight(node, 1);
        } endfor

        graph_partitioner partitioner;
        balance_configuration bc;
        bc.configurate_balance( partition_config, C);

        partitioner.perform_partitioning(partition_config, C);

        complete_boundary boundary(&C);
        boundary.build();

        cycle_refinement cr;

        partition_config.upper_bound_partition    = ceil(C.number_of_nodes()/(double)partition_config.k);
        cr.perform_refinement(partition_config, C, boundary);

        forall_nodes(C, node) {
                C.setNodeWeight(node, weights[node]);
        } endfor

        ofs.close();
        std::cout.rdbuf(backup);
}
