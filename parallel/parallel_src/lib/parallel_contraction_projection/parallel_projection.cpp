/******************************************************************************
 * parallel_projection.cpp
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

#include "parallel_projection.h"

parallel_projection::parallel_projection() {
                
}

parallel_projection::~parallel_projection() {
                
}

//issue recv before send
void parallel_projection::parallel_project( MPI_Comm communicator, parallel_graph_access & finer, parallel_graph_access & coarser ) {
        PEID rank, size;
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);
        
        NodeID divisor = ceil(coarser.number_of_global_nodes() / (double)size);

        m_messages.resize(size);

        std::unordered_map< NodeID, std::vector< NodeID > > cnode_to_nodes;
        forall_local_nodes(finer, node) {
                NodeID cnode = finer.getCNode(node);
                //std::cout <<  "cnode " <<  cnode  << std::endl;
                if( coarser.is_local_node_from_global_id(cnode) ) {
                        NodeID new_label = coarser.getNodeLabel(coarser.getLocalID(cnode));
                        finer.setNodeLabel(node, new_label);
                } else {
                        //we have to request it from another PE
                        PEID peID = cnode / divisor; // cnode is 

                        if( cnode_to_nodes.find( cnode ) == cnode_to_nodes.end()) {
                                m_messages[peID].push_back(cnode); // we are requesting the label of this node 
                        }

                        cnode_to_nodes[cnode].push_back(node);
                }
        } endfor

        for( PEID peID = 0; peID < size; peID++) {
                if( peID != rank ) {
                        if( m_messages[peID].size() == 0 ){
                                m_messages[peID].push_back(std::numeric_limits<NodeID>::max());
                        }

                        MPI_Request rq;
                        MPI_Isend( &m_messages[peID][0], 
                                                m_messages[peID].size(), 
                                                MPI_UNSIGNED_LONG_LONG, 
                                                peID, peID+size, communicator, &rq);
                }
        }

        std::vector< std::vector< NodeID > > out_messages;
        out_messages.resize(size);

        PEID counter = 0;
        while( counter < size - 1) {
                // wait for incomming message of an adjacent processor
                MPI_Status st;
                MPI_Probe(MPI_ANY_SOURCE, rank+size, communicator, &st);
                
                int message_length;
                MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);
                std::vector<NodeID> incmessage; incmessage.resize(message_length);

                MPI_Status rst;
                MPI_Recv( &incmessage[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, rank+size, communicator, &rst); 
                counter++;

                PEID peID = st.MPI_SOURCE;
                // now integrate the changes
                if( incmessage[0] == std::numeric_limits< NodeID >::max()) {
                        out_messages[peID].push_back(std::numeric_limits< NodeID >::max());
                        MPI_Request rq; 
                        MPI_Isend( &out_messages[peID][0], 
                                        out_messages[peID].size(), 
                                        MPI_UNSIGNED_LONG_LONG, 
                                        peID, peID+2*size, communicator, &rq);

                        continue; // nothing to do
                }


                for( int i = 0; i < message_length; i++) {
                        NodeID cnode = coarser.getLocalID(incmessage[i]);
                        out_messages[peID].push_back(coarser.getNodeLabel(cnode));
                }

                MPI_Request rq;
                MPI_Isend( &out_messages[peID][0], 
                                out_messages[peID].size(), 
                                MPI_UNSIGNED_LONG_LONG, 
                                peID, peID+2*size, communicator, &rq);

        }

        counter = 0;
        while( counter < size - 1) {
                // wait for incomming message of an adjacent processor
                MPI_Status st; ULONG tag = rank+2*size;
                MPI_Probe(MPI_ANY_SOURCE, tag, communicator, &st);

                int message_length;
                MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);
                std::vector<NodeID> incmessage; incmessage.resize(message_length);

                MPI_Status rst;
                MPI_Recv( &incmessage[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, tag, communicator, &rst); 
                counter++;

                // now integrate the changes
                if( incmessage[0] == std::numeric_limits< NodeID >::max()) {
                        continue; // nothing to do
                }

                PEID peID = st.MPI_SOURCE;
                for( ULONG i = 0; i < (ULONG)incmessage.size(); i++) {
                        std::vector< NodeID > & proj = cnode_to_nodes[m_messages[peID][i]];
                        NodeID label = incmessage[i];

                        for( ULONG j = 0; j < proj.size(); j++) {
                                finer.setNodeLabel(proj[j], label);
                        }
                }
        }

        finer.update_ghost_node_data_global(); // blocking
}

//initial assignment after initial partitioning
void parallel_projection::initial_assignment( parallel_graph_access & G, complete_graph_access & Q) {
        forall_local_nodes(G, node) {
                G.setNodeLabel(node, Q.getNodeLabel(G.getGlobalID(node)));
                if( G.is_interface_node(node) ) {
                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                if( !G.is_local_node( target ) ) {
                                        G.setNodeLabel(target, Q.getNodeLabel(G.getGlobalID(target)));
                                }
                        } endfor
                }
        } endfor
}
