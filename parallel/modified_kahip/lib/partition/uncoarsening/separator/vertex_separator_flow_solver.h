/******************************************************************************
 * vertex_separator_flow_solver.h
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

#ifndef VERTEX_SEPARATOR_FLOW_SOLVER_FLA4518Q
#define VERTEX_SEPARATOR_FLOW_SOLVER_FLA4518Q

#include "flow_solving_kernel/flow_solver.h"

class vertex_separator_flow_solver : public flow_solver {

public:
        vertex_separator_flow_solver();
        virtual ~vertex_separator_flow_solver();

        bool construct_flow_pb( const PartitionConfig & config, 
                                graph_access & G, 
                                PartitionID & lhs, 
                                PartitionID & rhs, 
                                std::vector<NodeID> & lhs_boundary_stripe,
                                std::vector<NodeID> & rhs_boundary_stripe,
                                std::vector<NodeID> & new_to_old_ids,              
                                long *n_ad, 
                                long* m_ad, 
                                node** nodes_ad, 
                                arc** arcs_ad, 
                                long ** cap_ad,
                                node** source_ad, 
                                node** sink_ad, 
                                long* node_min_ad,
                                EdgeID & no_edges_in_flow_graph); 


        void find_separator(const PartitionConfig & config, 
                            graph_access & G, 
                            PartitionID lhs, 
                            PartitionID rhs,  
                            boundary_starting_nodes start_nodes_lhs,
                            boundary_starting_nodes start_nodes_rhs,
                            std::vector<NodeID> & separator);
};


#endif /* end of include guard: VERTEX_SEPARATOR_FLOW_SOLVER_FLA4518Q */
