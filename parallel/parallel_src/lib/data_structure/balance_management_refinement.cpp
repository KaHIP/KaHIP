/******************************************************************************
 * balance_management_refinement.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "balance_management_refinement.h"
#include "parallel_graph_access.h"

balance_management_refinement::balance_management_refinement(parallel_graph_access * G, PartitionID total_num_labels)
: balance_management( G, total_num_labels) {
        m_total_block_weights.resize( total_num_labels );
        m_local_block_weights.resize( total_num_labels );

        for( long block = 0; block < (long) total_num_labels; block++) {
                m_local_block_weights[block] = 0;
                m_total_block_weights[block] = 0;
        }
        
        init();
}

balance_management_refinement::~balance_management_refinement() {
                
}

// init local and total block sizes
void balance_management_refinement::init() {
        forall_local_nodes((*m_G), node) {
                PartitionID label = m_G->getNodeLabel(node);
                m_local_block_weights[label] += m_G->getNodeWeight(node);
        } endfor
        update();
}

void balance_management_refinement::update() {
        MPI_Allreduce(&m_local_block_weights[0], &m_total_block_weights[0], 
                       m_total_num_labels, MPI_UNSIGNED_LONG_LONG, MPI_SUM, m_G->getCommunicator());
}
