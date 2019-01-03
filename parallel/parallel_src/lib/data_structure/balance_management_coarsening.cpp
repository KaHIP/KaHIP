/******************************************************************************
 * balance_management_coarsening.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <mpi.h>
#include "balance_management_coarsening.h"
#include "data_structure/parallel_graph_access.h"

balance_management_coarsening::balance_management_coarsening(parallel_graph_access * G, PartitionID total_num_labels) 
: balance_management( G, total_num_labels)
{
        init();
}

balance_management_coarsening::~balance_management_coarsening() {
                
}

void balance_management_coarsening::init(  ) {
        forall_local_nodes((*m_G), node) {
                PartitionID label = m_G->getNodeLabel(node);
                if( m_fuzzy_block_weights.find(label) == m_fuzzy_block_weights.end() ) {
                        m_fuzzy_block_weights[label] = 0;
                }

                m_fuzzy_block_weights[label] += m_G->getNodeWeight(node);
        } endfor

        forall_ghost_nodes((*m_G),node) {
                PartitionID label = m_G->getNodeLabel(node);
                if( m_fuzzy_block_weights.find(label) == m_fuzzy_block_weights.end() ) {
                        m_fuzzy_block_weights[label] = 0;
                }
                m_fuzzy_block_weights[label] += m_G->getNodeWeight(node);
        } endfor
}

void balance_management_coarsening::update( ) {
}
