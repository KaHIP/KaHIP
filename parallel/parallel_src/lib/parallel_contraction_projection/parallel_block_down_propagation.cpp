/******************************************************************************
 * parallel_block_down_propagation.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "parallel_block_down_propagation.h"

parallel_block_down_propagation::parallel_block_down_propagation() {
                
}

parallel_block_down_propagation::~parallel_block_down_propagation() {
                
}

void parallel_block_down_propagation::propagate_block_down( MPI_Comm communicator, PPartitionConfig & config, 
                                                            parallel_graph_access & G, 
                                                            parallel_graph_access & Q) {


        std::unordered_map< NodeID, NodeID > coarse_block_ids; 

        forall_local_nodes(G, node) {
                NodeID cur_cnode = G.getCNode( node );
                coarse_block_ids[cur_cnode] = G.getSecondPartitionIndex( node );
        } endfor

        PEID rank, size;
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);
        
        NodeID divisor          = ceil( Q.number_of_global_nodes()/(double)size);

        m_messages.resize(size);

        //now distribute the block idw
        //pack messages
        for( auto it = coarse_block_ids.begin(); it != coarse_block_ids.end(); it++) {
                NodeID node       = it->first;
                NodeID block      = it->second;
                PEID peID         = node / divisor;

                m_messages[ peID ].push_back( node );
                m_messages[ peID ].push_back( block );
        }

        for( PEID peID = 0; peID < size; peID++) {
                if( peID != rank ) {
                        if( m_messages[peID].size() == 0 ){
                                m_messages[peID].push_back(std::numeric_limits<NodeID>::max());
                        }

                        MPI_Request rq; 
                        MPI_Isend( &m_messages[peID][0], 
                                   m_messages[peID].size(), 
                                   MPI_UNSIGNED_LONG_LONG, 
                                   peID, peID+10*size, communicator, &rq );
                }
        }

        if( m_messages[ rank ].size() != 0 ) {
                for( ULONG i = 0; i < (ULONG)m_messages[rank].size()-1; i+=2) {
                        NodeID globalID   = m_messages[rank][i];
                        NodeID node       = Q.getLocalID(globalID);
                        NodeWeight block  = m_messages[rank][i+1];
                        Q.setSecondPartitionIndex(node , block);
                }
        }

        PEID counter = 0;
        while( counter < size - 1) {
                // wait for incomming message of an adjacent processor
                MPI_Status st;
                MPI_Probe(MPI_ANY_SOURCE, rank+10*size, communicator, &st);

                int message_length;
                MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);
                
                std::vector<NodeID> incmessage; incmessage.resize(message_length);

                MPI_Status rst;
                MPI_Recv( &incmessage[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, rank+10*size, communicator, &rst); 
                counter++;

                // now integrate the changes
                if( incmessage[0] == std::numeric_limits< NodeID >::max()) continue; // nothing to do

                for( ULONG i = 0; i < incmessage.size()-1; i+=2) {
                        NodeID globalID   = incmessage[i];
                        NodeWeight block  = incmessage[i+1];
                        NodeID node       = Q.getLocalID(globalID);
                        Q.setSecondPartitionIndex( node , block);
                }
        }

        update_ghost_nodes_blocks( communicator, Q );
}

void parallel_block_down_propagation::update_ghost_nodes_blocks( MPI_Comm communicator, parallel_graph_access & G ) {
        PEID rank, size;
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);
        
        m_send_buffers.resize(size);
        std::vector< bool > PE_packed(size, false);
        forall_local_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if( !G.is_local_node(target)  ) {
                                PEID peID = G.getTargetPE(target);
                                if( !PE_packed[peID] ) { // make sure a node is sent at most once
                                        m_send_buffers[peID].push_back(G.getGlobalID(node));
                                        m_send_buffers[peID].push_back(G.getSecondPartitionIndex(node));
                                        PE_packed[peID] = true;
                                }
                        }
                } endfor
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if( !G.is_local_node(target)  ) {
                                PE_packed[G.getTargetPE(target)] = false;
                        }
                } endfor
        } endfor

        //send all neighbors their packages using Isends
        //a neighbor that does not receive something gets a specific token
        for( PEID peID = 0; peID < (PEID)m_send_buffers.size(); peID++) {
                if( G.is_adjacent_PE(peID) ) {
                        //now we have to send a message
                        if( m_send_buffers[peID].size() == 0 ){
                                // length 1 encode no message
                                m_send_buffers[peID].push_back(0);
                        }

                        MPI_Request rq;
                        MPI_Isend( &m_send_buffers[peID][0], 
                                    m_send_buffers[peID].size(), MPI_UNSIGNED_LONG_LONG, peID, peID+11*size, communicator, &rq);
                }
        }

        //receive incomming
        PEID counter = 0;
        while( counter < G.getNumberOfAdjacentPEs()) {
                // wait for incomming message of an adjacent processor
                MPI_Status st;
                MPI_Probe(MPI_ANY_SOURCE, rank+11*size, communicator, &st);
                
                int message_length;
                MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);
                std::vector<NodeID> message; message.resize(message_length);

                MPI_Status rst;
                MPI_Recv( &message[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, rank+11*size, communicator, &rst); 
                counter++;

                // now integrate the changes
                if(message_length == 1) continue; // nothing to do

                for( int i = 0; i < message_length-1; i+=2) {
                        NodeID global_id   = message[i];
                        NodeWeight  block  = message[i+1];

                        G.setSecondPartitionIndex( G.getLocalID(global_id), block );
                }
        }

}
