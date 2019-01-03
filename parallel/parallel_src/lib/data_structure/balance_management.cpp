/******************************************************************************
 * balance_management.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "balance_management.h"
#include "data_structure/parallel_graph_access.h"

balance_management::balance_management( parallel_graph_access * G, NodeID total_num_labels ) : m_G ( G ), m_total_num_labels ( total_num_labels ) {
                
 
}

balance_management::~balance_management() {
                
}


