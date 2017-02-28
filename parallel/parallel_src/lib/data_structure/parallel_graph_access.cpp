/******************************************************************************
 * parallel_graph_access.cpp
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

#include "balance_management_coarsening.h"
#include "balance_management_refinement.h"
#include "parallel_graph_access.h"


ULONG parallel_graph_access::m_comm_rounds = 128; 
ULONG parallel_graph_access::m_comm_rounds_up = 128; 

parallel_graph_access::parallel_graph_access( MPI_Comm communicator ) : m_num_local_nodes(0), 
                                                 from(0), 
                                                 to(0),
                                                 m_num_ghost_nodes(0), m_max_node_degree(0), m_bm(NULL) {


                m_communicator = communicator;
                MPI_Comm_rank( m_communicator, &rank);
                MPI_Comm_size( m_communicator, &size);
                
                m_gnc = new ghost_node_communication(m_communicator);
                m_gnc->setGraphReference(this);
}

parallel_graph_access::~parallel_graph_access() {
        m_comm_rounds = std::min(m_comm_rounds, m_comm_rounds_up); 
        delete m_gnc;
        if ( m_bm ) delete m_bm;
}

void parallel_graph_access::init_balance_management( PPartitionConfig & config ) {
        if( m_bm != NULL ) {
                delete m_bm;
        }

        if( config.total_num_labels != config.k ) {
                m_bm = new balance_management_coarsening( this, config.total_num_labels );
        } else {
                m_bm = new balance_management_refinement( this, config.total_num_labels );
        }
}

void parallel_graph_access::update_non_contained_block_balance( PartitionID from, PartitionID to, NodeWeight node_weight) {
        m_bm->update_non_contained_block_balance( from, to, node_weight);
}
void parallel_graph_access::update_block_weights() {
        m_bm->update();
}

void parallel_graph_access::update_ghost_node_data( bool check_iteration_counter ) {
        m_gnc->update_ghost_node_data( check_iteration_counter );
}
void parallel_graph_access::update_ghost_node_data_global() {
        m_gnc->update_ghost_node_data_global();
}

void parallel_graph_access::update_ghost_node_data_finish() {
        m_gnc->update_ghost_node_data_finish();
}

void parallel_graph_access::set_comm_rounds(ULONG comm_rounds) {
        m_comm_rounds = comm_rounds;
        set_comm_rounds_up(comm_rounds);
}

void parallel_graph_access::set_comm_rounds_up(ULONG comm_rounds) {
        m_comm_rounds_up = comm_rounds;
}

