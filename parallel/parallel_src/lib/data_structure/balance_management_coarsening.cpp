/******************************************************************************
 * balance_management_coarsening.cpp
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
