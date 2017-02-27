/******************************************************************************
 * vertex_separator_algorithm.h
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

#ifndef VERTEX_SEPARTATOR_ALGORITHM_XUDNZZM8
#define VERTEX_SEPARTATOR_ALGORITHM_XUDNZZM8

#include <unordered_map>

#include "data_structure/graph_access.h"
#include "partition_config.h"
#include "uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"

class vertex_separator_algorithm {
        public:
                vertex_separator_algorithm();
                virtual ~vertex_separator_algorithm();

                void compute_vertex_separator(const PartitionConfig & config, 
                                              graph_access & G, 
                                              complete_boundary & boundary, 
                                              std::vector<NodeID> & overall_separator);

                void compute_vertex_separator_simple(const PartitionConfig & config, 
                                                     graph_access & G, 
                                                     complete_boundary & boundary, 
                                                     std::vector<NodeID> & overall_separator);

                void compute_vertex_separator(const PartitionConfig & config, 
                                              graph_access & G, 
                                              complete_boundary & boundary);

                //ASSERTIONS
                bool is_vertex_separator(graph_access & G, std::unordered_map<NodeID, bool> & separator);

};


#endif /* end of include guard: VERTEX_SEPARTATOR_ALGORITHM_XUDNZZM8 */
