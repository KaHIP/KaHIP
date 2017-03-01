/******************************************************************************
 * full_search_space.cpp
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

#include "full_search_space.h"

full_search_space::full_search_space(PartitionConfig & config, NodeID number_of_nodes) {
        m_unsucc_tries = 0;
        m_swap_lhs = 0;
        m_swap_rhs = 1;
        m_ub  = number_of_nodes*(number_of_nodes-1);
        m_ub /= 2;
        m_number_of_nodes = number_of_nodes;
}

full_search_space::~full_search_space() {
                
}

