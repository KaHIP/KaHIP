/******************************************************************************
 * active_block_quotient_graph_scheduler.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
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

#include "active_block_quotient_graph_scheduler.h"

active_block_quotient_graph_scheduler::active_block_quotient_graph_scheduler( const PartitionConfig & config, 
                                                                              QuotientGraphEdges & qgraph_edges, 
                                                                              unsigned int bank_account) :  
                                                                              m_quotient_graph_edges(qgraph_edges) {

        m_is_block_active.resize(config.k);
        for( unsigned int i = 0; i < m_is_block_active.size(); i++) {
                m_is_block_active[i] = true;
        }
         
        m_no_of_active_blocks = config.k;     
        init(); 
}

active_block_quotient_graph_scheduler::~active_block_quotient_graph_scheduler() {
                
}

