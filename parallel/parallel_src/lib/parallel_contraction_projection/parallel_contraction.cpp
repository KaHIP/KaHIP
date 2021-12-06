/******************************************************************************
 * parallel_contraction.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "parallel_contraction.h"
#include "data_structure/hashed_graph.h"
#include "tools/helpers.h"

parallel_contraction::parallel_contraction() {
                
}

parallel_contraction::~parallel_contraction() {
                
}

void parallel_contraction::contract_to_distributed_quotient( MPI_Comm communicator, PPartitionConfig & config, 
                                                             parallel_graph_access & G, 
                                                             parallel_graph_access & Q) {

        NodeID number_of_distinct_labels; // equals global number of coarse nodes

        // maps old ids to new ids in interval [0, ...., num_of_distinct_labels
        // and stores this information only for the local nodes 
        std::unordered_map< NodeID, NodeID > label_mapping;

        compute_label_mapping( communicator, G, number_of_distinct_labels, label_mapping);
        
        // compute the projection table
        G.allocate_node_to_cnode();
        forall_local_nodes(G, node) {
                G.setCNode( node, label_mapping[ G.getNodeLabel( node )]);
        } endfor

        get_nodes_to_cnodes_ghost_nodes( communicator, G );   

        //now we can really build the edges of the quotient graph
        hashed_graph hG;
        std::unordered_map< NodeID, NodeWeight > node_weights;

        build_quotient_graph_locally( G, number_of_distinct_labels, hG, node_weights);
        
        MPI_Barrier(communicator);

        m_messages.resize(0);
        std::vector< std::vector< NodeID > >(m_messages).swap(m_messages);
        m_out_messages.resize(0);
        std::vector< std::vector< NodeID > >(m_out_messages).swap(m_out_messages);
        m_send_buffers.resize(0); 
        std::vector< std::vector< NodeID > >(m_send_buffers).swap(m_send_buffers);

        redistribute_hased_graph_and_build_graph_locally( communicator, hG, node_weights, number_of_distinct_labels, Q );
        update_ghost_nodes_weights( communicator, Q );
}

void parallel_contraction::compute_label_mapping( MPI_Comm communicator, parallel_graph_access & G, 
                                                  NodeID & global_num_distinct_ids,
                                                  std::unordered_map< NodeID, NodeID > & label_mapping ) {
        PEID rank, size;
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);
        
        NodeID divisor  = ceil( G.number_of_global_nodes()/ (double)size);

        helpers helper;
        m_messages.resize(size);

        std::vector< std::unordered_map< NodeID, bool > > filter;
        filter.resize(size);
        forall_local_nodes(G, node) {
                PEID peID = G.getNodeLabel(node) / divisor;
                filter[ peID ][G.getNodeLabel(node)] = true;
        } endfor

        for( PEID peID = 0; peID < (PEID) size; peID++) {
                std::unordered_map< NodeID, bool >::iterator it;
                for( it = filter[peID].begin(); it != filter[peID].end(); it++) {
                        m_messages[peID].push_back(it->first);
                }
        }

        // now flood the network
        for( PEID peID = 0; peID < size; peID++) {
                if( peID != rank ) {
                        if( m_messages[peID].size() == 0 ){
                                m_messages[peID].push_back(std::numeric_limits<NodeID>::max());
                        }

                        MPI_Request rq; 
                        MPI_Isend( &m_messages[peID][0], 
                                    m_messages[peID].size(), 
                                    MPI_UNSIGNED_LONG_LONG, 
                                    peID, peID+4*size, communicator, &rq);
                }
        }
        std::vector< std::vector< NodeID > > local_labels_byPE;
        local_labels_byPE.resize(size);

        for( ULONG i = 0; i < m_messages[rank].size(); i++) {
                local_labels_byPE[rank].push_back(m_messages[rank][i]);
        }


        std::vector< std::vector< NodeID > >  inc_messages;
        inc_messages.resize(size);

        PEID counter = 0;
        while( counter < size - 1) {
                // wait for incomming message of an adjacent processor
                MPI_Status st;
                MPI_Probe(MPI_ANY_SOURCE, rank+4*size, communicator, &st);

                int message_length;
                MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);
                std::vector<NodeID> incmessage; incmessage.resize(message_length);

                MPI_Status rst;
                MPI_Recv( &incmessage[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, rank+4*size, communicator, &rst); 
                counter++;

                PEID peID = st.MPI_SOURCE;
                for( int i = 0; i < message_length; i++) {
                        inc_messages[peID].push_back(incmessage[i]);
                } // store those because we need to send them their mapping back

                // now integrate the changes
                if( incmessage[0] == std::numeric_limits< NodeID >::max()) continue; // nothing to do

                for( int i = 0; i < message_length; i++) {
                        local_labels_byPE[peID].push_back(incmessage[i]);
                }
        }

        std::vector< NodeID > local_labels;
        for( PEID peID = 0; peID < size; peID++) {
                for( ULONG i = 0; i < local_labels_byPE[peID].size(); i++) {
                        local_labels.push_back(local_labels_byPE[peID][i]);
                }
        }


        // filter duplicates locally
        helper.filter_duplicates( local_labels, 
                        [](const NodeID & lhs, const NodeID & rhs) -> bool { 
                        return (lhs <  rhs); 
                        }, 
                        [](const NodeID & lhs, const NodeID & rhs) -> bool { 
                        return (lhs ==  rhs); 
                        });
        //afterwards they are sorted!

        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // %%%%%%%%%%%%%%%%%%%%%%%Labels are unique on all PEs%%%%%%%%%%%%%%%%%%%%
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // now counting

        NodeID local_num_labels  = local_labels.size();
        NodeID prefix_sum        = 0;

        MPI_Scan(&local_num_labels, &prefix_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator); 

        global_num_distinct_ids = prefix_sum;
        // Broadcast global number of ids
        MPI_Bcast(&global_num_distinct_ids, 1, MPI_UNSIGNED_LONG_LONG, size-1, communicator); 

        NodeID num_smaller_ids = prefix_sum - local_num_labels;

        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // %%%%%Now Build the mapping and send information back to PEs%%%%%%%%%%%%
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        // build the mapping locally
        std::unordered_map< NodeID, NodeID > label_mapping_to_cnode;
        NodeID cur_id = num_smaller_ids;
        for( ULONG i = 0; i < local_labels.size(); i++) {
                label_mapping_to_cnode[local_labels[i]] = cur_id++;
        }

        // now send the processes the mapping back
        //std::vector< std::vector< NodeID > >  m_out_messages;
        m_out_messages.resize(size);

        for( PEID peID = 0; peID < (PEID)size; peID++) {
                if( peID == rank ) continue;

                if( inc_messages[peID][0] == std::numeric_limits<NodeID>::max()) {
                        m_out_messages[peID].push_back(std::numeric_limits<NodeID>::max());
                        continue;
                }

                for( ULONG i = 0; i < inc_messages[peID].size(); i++) {
                        m_out_messages[peID].push_back( label_mapping_to_cnode[ inc_messages[peID][i] ] );
                }
        }

        for( PEID peID = 0; peID < size; peID++) {
                if( peID != rank ) {
                        MPI_Request rq;
                        MPI_Isend( &m_out_messages[peID][0], 
                                    m_out_messages[peID].size(), 
                                    MPI_UNSIGNED_LONG_LONG, 
                                    peID, peID+5*size, communicator, &rq);
                }
        }

        // first the local labels 
        for( ULONG i = 0; i < m_messages[rank].size(); i++) {
                label_mapping[ m_messages[rank][i] ] = label_mapping_to_cnode[m_messages[rank][i]];
        }

        counter = 0;
        while( counter < size - 1) {
                // wait for incomming message of an adjacent processor
                MPI_Status st;
                MPI_Probe(MPI_ANY_SOURCE, rank+5*size, communicator, &st);
                
                int message_length;
                MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);
                std::vector<NodeID> incmessage; incmessage.resize(message_length);

                MPI_Status rst;
                MPI_Recv( &incmessage[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, rank+5*size, communicator, &rst); 
                counter++;

                // now integrate the changes
                if( incmessage[0] == std::numeric_limits< NodeID >::max()) continue; // nothing to do

                PEID peID = st.MPI_SOURCE;
                for( int i = 0; i < message_length; i++) {
                        label_mapping[ m_messages[peID][i] ] = incmessage[i];
                }
        }
}


void parallel_contraction::get_nodes_to_cnodes_ghost_nodes( MPI_Comm communicator, parallel_graph_access & G ) {
        PEID rank, size;
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);
        
        std::vector< bool > PE_packed( size, false );
        m_send_buffers.resize( size );

        forall_local_nodes(G, node) {
                if(G.is_interface_node(node)) {
                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                if( !G.is_local_node(target)  ) {
                                        PEID peID = G.getTargetPE(target);
                                        if( !PE_packed[peID] ) { // make sure a node is sent at most once
                                                m_send_buffers[peID].push_back(G.getGlobalID(node));
                                                m_send_buffers[peID].push_back(G.getCNode(node));
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
                }
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
                                    m_send_buffers[peID].size(), MPI_UNSIGNED_LONG_LONG, peID, peID+6*size, communicator, &rq);
                }
        }

        ////receive incomming
        PEID num_adjacent = G.getNumberOfAdjacentPEs();
        PEID counter = 0;
        while( counter < num_adjacent ) {
                // wait for incomming message of an adjacent processor
                MPI_Status st;
                MPI_Probe(MPI_ANY_SOURCE, rank+6*size, communicator, &st);
                
                int message_length;
                MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);
                std::vector<NodeID> message; message.resize(message_length);

                MPI_Status rst;
                MPI_Recv( &message[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, rank+6*size, communicator, &rst); 
                counter++;

                // now integrate the changes
                if(message_length == 1) continue; // nothing to do

                for( int i = 0; i < message_length-1; i+=2) {
                        NodeID global_id = message[i];
                        NodeID cnode     = message[i+1];

                        G.setCNode( G.getLocalID(global_id), cnode);
                }
        }
}


void parallel_contraction::build_quotient_graph_locally( parallel_graph_access & G, 
                                                         NodeID number_of_distinct_labels, 
                                                         hashed_graph & hG, 
                                                         std::unordered_map< NodeID, NodeWeight > & node_weights) {
        forall_local_nodes(G, node) {
                NodeID cur_cnode = G.getCNode( node );
                if( node_weights.find(cur_cnode) == node_weights.end()) {
                        node_weights[cur_cnode] = 0;
                }

                node_weights[cur_cnode] += G.getNodeWeight( node );

                forall_out_edges(G, e, node) {
                        NodeID target       = G.getEdgeTarget(e);
                        NodeID target_cnode = G.getCNode(target);
                        if( cur_cnode != target_cnode ) {
                                // update the edge
                                hashed_edge he;
                                he.k            = number_of_distinct_labels;
                                he.source       = cur_cnode;
                                he.target       = target_cnode;

                                hG[he].weight  += G.getEdgeWeight(e);
                        }
                } endfor
        } endfor
}



void parallel_contraction::redistribute_hased_graph_and_build_graph_locally( MPI_Comm communicator, hashed_graph &  hG, 
                                                                             std::unordered_map< NodeID, NodeWeight > & node_weights,
                                                                             NodeID number_of_cnodes, 
                                                                             parallel_graph_access & Q  ) {
        PEID rank, size;
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);
        
        NodeID divisor          = ceil( number_of_cnodes/(double)size);

        //std::vector< std::vector< NodeID > >  messages;
        m_messages.resize(size);

        //build messages
        hashed_graph::iterator it;
        for( it = hG.begin(); it != hG.end(); it++) {
                data_hashed_edge & e = it->second;
                hashed_edge he       = it->first;

                PEID peID = he.source / divisor;
                m_messages[ peID ].push_back( he.source );
                m_messages[ peID ].push_back( he.target );
                m_messages[ peID ].push_back( e.weight );

                peID = he.target / divisor;
                m_messages[ peID ].push_back( he.target );
                m_messages[ peID ].push_back( he.source );
                m_messages[ peID ].push_back( e.weight );
        }

        // now flood the network
        for( PEID peID = 0; peID < size; peID++) {
                if( peID != rank ) {
                        if( m_messages[peID].size() == 0 ){
                                m_messages[peID].push_back(std::numeric_limits<NodeID>::max());
                        }

                        MPI_Request rq;
                        MPI_Isend( &m_messages[peID][0], 
                                    m_messages[peID].size(), 
                                    MPI_UNSIGNED_LONG_LONG, 
                                    peID, peID+7*size, communicator, &rq);
                }
        }

        // build the local part of the graph
        //
        std::vector< std::vector< NodeID > > local_msg_byPE;
        local_msg_byPE.resize(size);


        if( m_messages[ rank ].size() != 0 ) {
                local_msg_byPE[rank] = m_messages[rank];
        }
       
        PEID counter = 0;
        while( counter < size - 1) {
                // wait for incomming message of an adjacent processor
                MPI_Status st;
                MPI_Probe(MPI_ANY_SOURCE, rank+7*size, communicator, &st);
                
                int message_length;
                MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);
                std::vector<NodeID> incmessage; incmessage.resize(message_length);

                MPI_Status rst;
                MPI_Recv( &incmessage[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, rank+7*size, communicator, &rst); 
                counter++;

                // now integrate the changes
                if( incmessage[0] == std::numeric_limits< NodeID >::max()) continue; // nothing to do


                PEID peID = rst.MPI_SOURCE;
                local_msg_byPE[peID] = incmessage;
        }

        hashed_graph local_graph;
        for( PEID peID = 0; peID < size; peID++) {
                if(local_msg_byPE[peID].size() > 0) {
                for( ULONG i = 0; i < local_msg_byPE[peID].size()-2; i+=3) {
                        hashed_edge he;
                        he.k = number_of_cnodes;
                        he.source = local_msg_byPE[peID][i];
                        he.target = local_msg_byPE[peID][i+1];

                        local_graph[he].weight += local_msg_byPE[peID][i+2];
                }}
        }


        ULONG from = rank     * ceil(number_of_cnodes / (double)size);
        ULONG to   = (rank+1) * ceil(number_of_cnodes / (double)size) - 1;
        // handle the case where we dont have local edges
        from = std::min(from, number_of_cnodes);
        to   = std::min(to, number_of_cnodes - 1);
        ULONG local_num_cnodes = to - from + 1;

        std::vector < std::vector< std::pair<NodeID, NodeWeight > > > sorted_graph;
        sorted_graph.resize( local_num_cnodes );

        EdgeID edge_counter = 0;
        for( it = local_graph.begin(); it != local_graph.end(); it++) {
                data_hashed_edge & e = it->second;
                hashed_edge he       = it->first;

                if( from <= he.target && he.target <= to) {
                        std::pair< NodeID, NodeWeight > edge;
                        edge.first  = he.target;
                        edge.second = e.weight/4;

                        std::pair< NodeID, NodeWeight > e_bar;
                        e_bar.first  = he.source;
                        e_bar.second = e.weight/4;

                        sorted_graph[ he.target - from ].push_back( e_bar);
                        sorted_graph[ he.source - from ].push_back( edge );
                        edge_counter+=2;
                } else {
                        std::pair< NodeID, NodeWeight > edge;
                        edge.first  = he.target;
                        edge.second = e.weight/2;
                        sorted_graph[ he.source - from ].push_back( edge );
                        edge_counter++;
                }
        }
 
        ULONG global_edges = 0;
        MPI_Allreduce(&edge_counter, &global_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);

        Q.start_construction(local_num_cnodes, edge_counter, number_of_cnodes, global_edges);
        Q.set_range(from, to);

        std::vector< NodeID > vertex_dist( size+1, 0 );
        for( PEID peID = 0; peID <= size; peID++) {
                vertex_dist[peID] = std::min(number_of_cnodes, (NodeID) (peID * ceil(number_of_cnodes / (double)size))); // from positions
        }
        //vertex_dist[size] = std::min(to, number_of_cnodes - 1);
        Q.set_range_array(vertex_dist);

        for (NodeID i = 0; i < local_num_cnodes; ++i) {
                NodeID node = Q.new_node();
                NodeID globalID = from+node;
                Q.setNodeWeight(node, 0); 
                Q.setNodeLabel(node, globalID);
               
                for( EdgeID e = 0; e < sorted_graph[node].size(); e++) {
                        NodeID target = sorted_graph[node][e].first;
                        EdgeID e_bar = Q.new_edge(node, target);
                        Q.setEdgeWeight(e_bar, sorted_graph[node][e].second);
                }
        }

        Q.finish_construction();

        for( PEID peID = 0; peID < size; peID++) {
                m_messages[peID].clear();
        }
        //now distribute the node weights
        //pack messages
        std::unordered_map< NodeID, NodeWeight >::iterator wit;
        for( wit = node_weights.begin(); wit != node_weights.end(); wit++) {
                NodeID node       = wit->first;
                NodeWeight weight = wit->second;
                PEID peID         = node / divisor;

                m_messages[ peID ].push_back( node );
                m_messages[ peID ].push_back( weight );
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
                                    peID, peID+8*size, communicator, &rq);
                }
        }

        if( m_messages[ rank ].size() != 0 ) {
                for( ULONG i = 0; i < m_messages[rank].size()-1; i+=2) {
                        NodeID globalID   = m_messages[rank][i];
                        NodeID node       = globalID - from;
                        NodeWeight weight = m_messages[rank][i+1];
                        Q.setNodeWeight( node , Q.getNodeWeight(node) + weight);
                }
        }

        counter = 0;
        while( counter < size - 1) {
                // wait for incomming message of an adjacent processor
                MPI_Status st;
                MPI_Probe(MPI_ANY_SOURCE, rank+8*size, communicator, &st);

                int message_length;
                MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);
                std::vector<NodeID> incmessage; incmessage.resize(message_length);

                MPI_Status rst;
                MPI_Recv( &incmessage[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, rank+8*size, communicator, &rst); 
                counter++;

                // now integrate the changes
                if( incmessage[0] == std::numeric_limits< NodeID >::max()) continue; // nothing to do

                for( ULONG i = 0; i < incmessage.size()-1; i+=2) {
                        NodeID globalID   = incmessage[i];
                        NodeWeight weight = incmessage[i+1];
                        NodeID node       = globalID - from;
                        Q.setNodeWeight( node , Q.getNodeWeight(node) + weight);
                }
        }
}


void parallel_contraction::update_ghost_nodes_weights( MPI_Comm communicator, parallel_graph_access & G ) {
        PEID rank, size;
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);
        
               //std::vector< std::vector< NodeID > > send_buffers; // buffers to send messages
        m_send_buffers.resize(size);
        std::vector< bool > PE_packed(size, false);
        forall_local_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if( !G.is_local_node(target)  ) {
                                PEID peID = G.getTargetPE(target);
                                if( !PE_packed[peID] ) { // make sure a node is sent at most once
                                        m_send_buffers[peID].push_back(G.getGlobalID(node));
                                        m_send_buffers[peID].push_back(G.getNodeWeight(node));
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
                                    m_send_buffers[peID].size(), MPI_UNSIGNED_LONG_LONG, peID, peID+9*size, communicator, &rq);
                }
        }

        //receive incomming
        PEID counter = 0;
        while( counter < G.getNumberOfAdjacentPEs()) {
                // wait for incomming message of an adjacent processor
                MPI_Status st;
                MPI_Probe(MPI_ANY_SOURCE, rank+9*size, communicator, &st);
                
                int message_length;
                MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);
                std::vector<NodeID> message; message.resize(message_length);

                MPI_Status rst;
                MPI_Recv( &message[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, rank+9*size, communicator, &rst); 
                counter++;

                // now integrate the changes
                if(message_length == 1) continue; // nothing to do

                for( int i = 0; i < message_length-1; i+=2) {
                        NodeID global_id   = message[i];
                        NodeWeight  weight = message[i+1];

                        G.setNodeWeight( G.getLocalID(global_id), weight);
                }
        }

}
