/******************************************************************************
 * full_search_space_pruned.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "full_search_space_pruned.h"

full_search_space_pruned::full_search_space_pruned(PartitionConfig & config, NodeID number_of_nodes) {
        m_ub         = config.search_space_s*(config.search_space_s-1);
        m_ub /= 2;
        this->config   = config;
        m_internal_k   = 0;
        m_unsucc_tries = 0;
        m_number_of_nodes = number_of_nodes;

        for( unsigned int k = 0; k < ceil(number_of_nodes/config.search_space_s); k++) {
                NodeID lb = k*config.search_space_s;
                m_search_space_pointers.push_back( std::pair< NodeID, NodeID>(lb, lb+1));
        }
}

full_search_space_pruned::~full_search_space_pruned() {
                
}

