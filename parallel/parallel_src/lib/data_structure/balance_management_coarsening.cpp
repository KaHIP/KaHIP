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
        m_local_from = m_G->get_from_range();
        m_local_weights.assign(m_G->number_of_local_nodes(), 0);
        m_fuzzy_block_weights.reserve(m_G->number_of_ghost_nodes());
        forall_local_nodes((*m_G), node) {
                PartitionID label = m_G->getNodeLabel(node);
                local_or_map_ref(label) += m_G->getNodeWeight(node);
        } endfor

        const NodeID ghost_end = m_G->number_of_local_nodes() + 1 + m_G->number_of_ghost_nodes();
        forall_ghost_nodes((*m_G),node) {
                if( node + 8 < ghost_end ) prefetch_block(m_G->getNodeLabel(node + 8)); // random keys: hide the miss
                PartitionID label = m_G->getNodeLabel(node);
                local_or_map_ref(label) += m_G->getNodeWeight(node);
        } endfor
}

void balance_management_coarsening::update( ) {
}
