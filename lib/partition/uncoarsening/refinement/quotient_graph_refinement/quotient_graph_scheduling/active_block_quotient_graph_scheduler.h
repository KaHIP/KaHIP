/******************************************************************************
 * active_block_quotient_graph_scheduler.h 
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

#ifndef ACTIVE_BLOCK_QUOTIENT_GRAPH_SCHEDULER_2QATIGSY
#define ACTIVE_BLOCK_QUOTIENT_GRAPH_SCHEDULER_2QATIGSY

#include <unordered_map>

#include "partition_config.h"
#include "quotient_graph_scheduling.h"
#include "random_functions.h"

class active_block_quotient_graph_scheduler : public quotient_graph_scheduling {
        public:
                active_block_quotient_graph_scheduler( const PartitionConfig & config,
                                QuotientGraphEdges & qgraph_edges, 
                                unsigned int bank_account);

                virtual ~active_block_quotient_graph_scheduler();

                virtual bool hasFinished();
                virtual boundary_pair & getNext();
                virtual void pushStatistics(qgraph_edge_statistics & statistic);
                virtual void init();

                void activate_blocks(std::unordered_map<PartitionID, PartitionID> & blocks);

        private: 
                QuotientGraphEdges & m_quotient_graph_edges;
                QuotientGraphEdges   m_active_quotient_graph_edges;
                PartitionID          m_no_of_active_blocks;
                std::vector<bool>    m_is_block_active;
};

inline void active_block_quotient_graph_scheduler::init() {
        m_no_of_active_blocks = 0;
        m_active_quotient_graph_edges.clear();

        for( unsigned int i = 0; i < m_quotient_graph_edges.size(); i++) {
                PartitionID lhs = m_quotient_graph_edges[i].lhs;                      
                PartitionID rhs = m_quotient_graph_edges[i].rhs;  

                if(m_is_block_active[lhs]) m_no_of_active_blocks++;
                if(m_is_block_active[rhs]) m_no_of_active_blocks++;

                if(m_is_block_active[lhs] || m_is_block_active[rhs]) {
                        m_active_quotient_graph_edges.push_back(m_quotient_graph_edges[i]);
                }
        }

        random_functions::permutate_vector_good_small(m_active_quotient_graph_edges);

        for( unsigned int i = 0; i < m_is_block_active.size(); i++) {
                m_is_block_active[i] = false;
        }
}

inline bool active_block_quotient_graph_scheduler::hasFinished( ) {
        if(m_active_quotient_graph_edges.empty()) {
                init();
        }

        return m_no_of_active_blocks == 0;        
}

inline boundary_pair & active_block_quotient_graph_scheduler::getNext( ) {
        boundary_pair & ret_value = m_active_quotient_graph_edges.back();
        m_active_quotient_graph_edges.pop_back();

        return ret_value; 
}

inline void active_block_quotient_graph_scheduler::pushStatistics(qgraph_edge_statistics & statistic) {
        if(statistic.something_changed) {
                m_is_block_active[statistic.pair->lhs] = true;
                m_is_block_active[statistic.pair->rhs] = true;
        }
}

inline void active_block_quotient_graph_scheduler::activate_blocks(std::unordered_map<PartitionID, PartitionID> & blocks) {
        std::unordered_map<PartitionID, PartitionID>::iterator it;
        for(it = blocks.begin(); it != blocks.end(); ++it) {
             m_is_block_active[it->first] = true;
        }
}

#endif /* end of include guard: ACTIVE_BLOCK_QUOTIENT_GRAPH_SCHEDULER_2QATIGSY */
