/******************************************************************************
 * balance_management_refinement.cpp
 *
 * Source of KaHIP -- Karlsruhe High Quality Graph Partitioning 
 ******************************************************************************
 * Copyright (C) 2017 Christian Schulz 
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
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
