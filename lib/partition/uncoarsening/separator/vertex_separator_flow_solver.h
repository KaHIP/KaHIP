/******************************************************************************
 * vertex_separator_flow_solver.h 
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

#ifndef VERTEX_SEPARATOR_FLOW_SOLVER_FLA4518Q
#define VERTEX_SEPARATOR_FLOW_SOLVER_FLA4518Q

#include "definitions.h"
#include "partition_config.h"
#include "data_structure/flow_graph.h"

class vertex_separator_flow_solver  {

public:
        vertex_separator_flow_solver();
        virtual ~vertex_separator_flow_solver();

        bool build_flow_pb( const PartitionConfig & config, 
                               graph_access & G, 
                               PartitionID & lhs, 
                               PartitionID & rhs, 
                               std::vector<NodeID> & lhs_nodes,
                               std::vector<NodeID> & rhs_nodes,
                               std::vector<NodeID> & new_to_old_ids,              
                               flow_graph & fG); 



        void find_separator(const PartitionConfig & config, 
                            graph_access & G, 
                            PartitionID lhs, 
                            PartitionID rhs,  
                            boundary_starting_nodes start_nodes_lhs,
                            boundary_starting_nodes start_nodes_rhs,
                            std::vector<NodeID> & separator);
};


#endif /* end of include guard: VERTEX_SEPARATOR_FLOW_SOLVER_FLA4518Q */
