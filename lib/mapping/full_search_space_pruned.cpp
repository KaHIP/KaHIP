/******************************************************************************
 * full_search_space_pruned.cpp
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

