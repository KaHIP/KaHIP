/******************************************************************************
 * parallel_contraction.h
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

#ifndef PARALLEL_CONTRACTION_64O127GD
#define PARALLEL_CONTRACTION_64O127GD

#include "data_structure/hashed_graph.h"
#include "data_structure/parallel_graph_access.h"
#include "partition_config.h"

class parallel_contraction {
public:
        parallel_contraction();
        virtual ~parallel_contraction();

        void contract_to_distributed_quotient( MPI_Comm communicator, PPartitionConfig & config, 
                                               parallel_graph_access & G, 
                                               parallel_graph_access & Q);
private:
        // compute mapping of labels id into contiguous intervall [0, ...., num_lables)
        void compute_label_mapping( MPI_Comm communicator, parallel_graph_access & G, 
                                    NodeID & global_num_distinct_ids,
                                    std::unordered_map< NodeID, NodeID > & label_mapping);

        void get_nodes_to_cnodes_ghost_nodes( MPI_Comm communicator, parallel_graph_access & G );   

	void build_quotient_graph_locally( parallel_graph_access & G, 
                                           NodeID number_of_distinct_labels, 
					   hashed_graph & hG, 
					   std::unordered_map< NodeID, NodeWeight > & node_weights);

        void redistribute_hased_graph_and_build_graph_locally( MPI_Comm communicator, hashed_graph &  hG, 
                                                               std::unordered_map< NodeID, NodeWeight > & node_weights,
                                                               NodeID number_of_cnodes,
                                                               parallel_graph_access & Q);

        void update_ghost_nodes_weights( MPI_Comm communicator, parallel_graph_access & G ); 

        // some send buffers
        std::vector< std::vector< NodeID > >  m_messages;
        std::vector< std::vector< NodeID > >  m_out_messages;
        std::vector< std::vector< NodeID > >  m_send_buffers; // buffers to send messages
};


#endif /* end of include guard: PARALLEL_CONTRACTION_64O127GD */
