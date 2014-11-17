/******************************************************************************
 * simple_quotient_graph_scheduler.h 
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

#ifndef SIMPLE_QUOTIENT_GRAPH_SCHEDULER_YG9BEBH0
#define SIMPLE_QUOTIENT_GRAPH_SCHEDULER_YG9BEBH0

#include "partition_config.h"
#include "quotient_graph_scheduling.h"

class simple_quotient_graph_scheduler : public quotient_graph_scheduling  { 
public:
        simple_quotient_graph_scheduler(PartitionConfig & config, 
                                        QuotientGraphEdges & qgraph_edges,  
                                        unsigned int account);

        virtual ~simple_quotient_graph_scheduler();

        virtual bool hasFinished();
        virtual boundary_pair & getNext();
        virtual void pushStatistics(qgraph_edge_statistics & statistic) {};

private: 
        QuotientGraphEdges m_quotient_graph_edges_pool;
};

inline bool simple_quotient_graph_scheduler::hasFinished( ) {
        return m_quotient_graph_edges_pool.empty();                
}

inline boundary_pair & simple_quotient_graph_scheduler::getNext( ) {
        boundary_pair & ret_value = m_quotient_graph_edges_pool.back();
        m_quotient_graph_edges_pool.pop_back();
        return ret_value; 
}

#endif /* end of include guard: SIMPLE_QUOTIENT_GRAPH_SCHEDULER_YG9BEBH0 */
