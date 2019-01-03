/******************************************************************************
 * mapping_algorithms.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <algorithm>
#include "communication_graph_search_space.h"
#include "construct_distance_matrix.h"
#include "construct_mapping.h"
#include "full_search_space.h"
#include "full_search_space_pruned.h"
#include "local_search_mapping.h"
#include "mapping_algorithms.h"
#include "partition/partition_config.h"
#include "tools/random_functions.h"

mapping_algorithms::mapping_algorithms() {

}

mapping_algorithms::~mapping_algorithms() {

}

void mapping_algorithms::construct_a_mapping( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank) {
        PRINT(std::cout <<  "computing distance matrix "  << std::endl;)
        construct_distance_matrix cdm;
        cdm.construct_matrix( config, D );

        t.restart();
        construct_mapping cm;
        cm.construct_initial_mapping( config, C, D, perm_rank);
        PRINT(std::cout <<  "construction took " <<  t.elapsed() << std::endl;)
        
        t.restart();
        local_search_mapping lsm;
        switch( config.ls_neighborhood ) {
                case NSQUARE:
                        lsm.perform_local_search< full_search_space > ( config, C, D, perm_rank);
                        break;
                case NSQUAREPRUNED:
                        lsm.perform_local_search< full_search_space_pruned > ( config, C, D, perm_rank);
                        break;
                case COMMUNICATIONGRAPH:
                        lsm.perform_local_search< communication_graph_search_space > ( config, C, D, perm_rank);
                        break;
        }

        PRINT(std::cout <<  "local search took " <<  t.elapsed()  << std::endl;)
}

void mapping_algorithms::graph_to_matrix( graph_access & C, matrix & C_bar) {
        for( unsigned int i = 0; i < C.number_of_nodes(); i++) {
                for( unsigned int j = 0; j < C.number_of_nodes(); j++) {
                        C_bar.set_xy(i,j,0);
                }
        }

        forall_nodes(C, node) {
                forall_out_edges(C, e, node) {
                        NodeID target = C.getEdgeTarget(e);
                        C_bar.set_xy(node, target, C.getEdgeWeight(e));
                } endfor
        } endfor
}
