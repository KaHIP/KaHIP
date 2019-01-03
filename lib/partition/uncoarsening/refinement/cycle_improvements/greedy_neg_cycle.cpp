/******************************************************************************
 * greedy_neg_cycle.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "greedy_neg_cycle.h"

greedy_neg_cycle::greedy_neg_cycle( PartitionConfig & partition_config ) {
        commons = new kway_graph_refinement_commons( partition_config );
}

greedy_neg_cycle::~greedy_neg_cycle() {
        if( commons != NULL ) delete commons;
}

