/******************************************************************************
 * full_search_space.h
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

#ifndef FULL_SEARCH_SPACE_H49ZQ8A4
#define FULL_SEARCH_SPACE_H49ZQ8A4

#include <utility>
#include "data_structure/graph_access.h"
#include "partition_config.h"

class full_search_space {
public:
        full_search_space(PartitionConfig & config, NodeID number_of_nodes);
        virtual ~full_search_space();

        void set_graph_ref( graph_access * C ) {}

        bool done() {
             return !(m_unsucc_tries < m_ub);
        }; // are we done?

        void commit_status( bool success ) {
                if(success) m_unsucc_tries = 0;
                else m_unsucc_tries++;
	}

        std::pair< NodeID, NodeID > nextPair() {
                std::pair< NodeID, NodeID > ret_value(m_swap_lhs, m_swap_rhs);
                if(  m_swap_rhs+1 < m_number_of_nodes )
                        m_swap_rhs++;
                else {
                        if( m_swap_lhs + 2 < m_number_of_nodes ) {
                                m_swap_lhs += 1;
                                m_swap_rhs = m_swap_lhs + 1;
                        } else {
                                m_swap_lhs = 0;
                                m_swap_rhs = 1;
                        }
                }

                return ret_value;
        }
private:
        NodeID m_ub;
        NodeID m_swap_lhs;
        NodeID m_swap_rhs;
        NodeID m_unsucc_tries;
        NodeID m_number_of_nodes;
};


#endif /* end of include guard: FULL_SEARCH_SPACE_H49ZQ8A4 */
