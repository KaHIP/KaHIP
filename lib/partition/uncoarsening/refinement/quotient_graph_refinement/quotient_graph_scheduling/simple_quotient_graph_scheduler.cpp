/******************************************************************************
 * simple_quotient_graph_scheduler.cpp 
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

#include "random_functions.h"
#include "simple_quotient_graph_scheduler.h"

simple_quotient_graph_scheduler::simple_quotient_graph_scheduler(PartitionConfig & config, 
                                                                 QuotientGraphEdges & qgraph_edges,  
                                                                 unsigned int account) {
        unsigned added_edges = 0;
        for( unsigned i = 0; i < (unsigned)ceil(config.bank_account_factor) && added_edges <= account; i++) {
                random_functions::permutate_vector_good_small(qgraph_edges);               
                 for( unsigned i = 0; i < qgraph_edges.size() && added_edges <= account; i++) {
                        m_quotient_graph_edges_pool.push_back(qgraph_edges[i]);
                        added_edges++;
                 }
        }

}

simple_quotient_graph_scheduler::~simple_quotient_graph_scheduler() {
                
}

