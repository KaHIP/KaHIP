/******************************************************************************
 * parallel_contraction.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "parallel_contraction.h"
#include "data_structure/ordered_edge_map.h"
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
        open_hash_map< NodeID > label_mapping;

        compute_label_mapping( communicator, G, number_of_distinct_labels, label_mapping);
        
        // compute the projection table
        G.allocate_node_to_cnode();
        const NodeID ahead = 8; // hash map lookups below prefetch their slot this many entries ahead
        forall_local_nodes(G, node) {
                if( node + ahead < G.number_of_local_nodes() ) label_mapping.prefetch( G.getNodeLabel( node + ahead ));
                G.setCNode( node, label_mapping.ref( G.getNodeLabel( node )));
        } endfor

        get_nodes_to_cnodes_ghost_nodes( communicator, G );   

        //now we can really build the edges of the quotient graph
        ordered_edge_map hG(number_of_distinct_labels);
        open_hash_map< NodeWeight > node_weights;

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
                                                  open_hash_map< NodeID > & label_mapping ) {
        PEID rank, size;
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);
        
        NodeID divisor  = ceil( G.number_of_global_nodes()/ (double)size);

        helpers helper;
        m_messages.resize(size);

        // distinct labels of the local nodes, grouped by the PE that numbers them
        // (the order inside a message is irrelevant: the answers are aligned with it)
        open_hash_map< char > filter;
        filter.reserve(G.number_of_local_nodes());
        forall_local_nodes(G, node) {
                filter.ref(G.getNodeLabel(node)) = 1;
        } endfor
        filter.for_each([&](NodeID label, char) {
                m_messages[ label / divisor ].push_back(label);
        });

        // now flood the network
        std::vector<MPI_Request> requests(size, MPI_REQUEST_NULL);
        for( PEID peID = 0; peID < size; peID++) {
                if( peID != rank ) {
                        if( m_messages[peID].size() == 0 ){
                                m_messages[peID].push_back(std::numeric_limits<NodeID>::max());
                        }

                        MPI_Isend( &m_messages[peID][0], 
                                    m_messages[peID].size(), 
                                    MPI_UNSIGNED_LONG_LONG, 
                                    peID, peID+4*size, communicator, &requests[peID]);
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
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

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
        open_hash_map< NodeID > label_mapping_to_cnode;
        label_mapping_to_cnode.reserve(local_labels.size());
        NodeID cur_id = num_smaller_ids;
        for( ULONG i = 0; i < local_labels.size(); i++) {
                label_mapping_to_cnode.ref(local_labels[i]) = cur_id++;
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
                        if( i + 8 < inc_messages[peID].size() ) label_mapping_to_cnode.prefetch( inc_messages[peID][i + 8] );
                        m_out_messages[peID].push_back( label_mapping_to_cnode.ref( inc_messages[peID][i] ) );
                }
        }

        for( PEID peID = 0; peID < size; peID++) {
                if( peID != rank ) {
                        MPI_Isend( &m_out_messages[peID][0], 
                                    m_out_messages[peID].size(), 
                                    MPI_UNSIGNED_LONG_LONG, 
                                    peID, peID+5*size, communicator, &requests[peID]);
                }
        }

        // first the local labels 
        for( ULONG i = 0; i < m_messages[rank].size(); i++) {
                if( i + 8 < m_messages[rank].size() ) {
                        label_mapping.prefetch( m_messages[rank][i + 8] );
                        label_mapping_to_cnode.prefetch( m_messages[rank][i + 8] );
                }
                label_mapping.ref( m_messages[rank][i] ) = label_mapping_to_cnode.ref( m_messages[rank][i] );
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
                        if( i + 8 < message_length ) label_mapping.prefetch( m_messages[peID][i + 8] );
                        label_mapping.ref( m_messages[peID][i] ) = incmessage[i];
                }
        }
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}


void parallel_contraction::get_nodes_to_cnodes_ghost_nodes( MPI_Comm communicator, parallel_graph_access & G ) {
        PEID rank, size;
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);
        
        // PE_packed_by[peID] == node + 1 iff node was already sent to peID (replaces resetting flags in a second edge scan)
        std::vector< NodeID > PE_packed_by( size, 0 );
        m_send_buffers.resize( size );

        forall_local_nodes(G, node) {
                if(G.is_interface_node(node)) {
                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                if( !G.is_local_node(target)  ) {
                                        PEID peID = G.getTargetPE(target);
                                        if( PE_packed_by[peID] != node + 1 ) { // make sure a node is sent at most once
                                                m_send_buffers[peID].push_back(G.getGlobalID(node));
                                                m_send_buffers[peID].push_back(G.getCNode(node));
                                                PE_packed_by[peID] = node + 1;
                                        }
                                }
                        } endfor
                }
        } endfor

        //send all neighbors their packages using Isends
        //a neighbor that does not receive something gets a specific token
        std::vector<MPI_Request> requests(m_send_buffers.size(), MPI_REQUEST_NULL);
        for( PEID peID = 0; peID < (PEID)m_send_buffers.size(); peID++) {
                if( G.is_adjacent_PE(peID) ) {
                        //now we have to send a message
                        if( m_send_buffers[peID].size() == 0 ){
                                // length 1 encode no message
                                m_send_buffers[peID].push_back(0);
                        }

                        MPI_Isend( &m_send_buffers[peID][0], 
                                    m_send_buffers[peID].size(), MPI_UNSIGNED_LONG_LONG, peID, peID+6*size, communicator, &requests[peID]);
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
                        if( i + 16 < message_length-1 ) G.prefetch_local_id( message[i + 16] );
                        NodeID global_id = message[i];
                        NodeID cnode     = message[i+1];

                        G.setCNode( G.getLocalID(global_id), cnode);
                }
        }
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}


void parallel_contraction::build_quotient_graph_locally( parallel_graph_access & G, 
                                                         NodeID number_of_distinct_labels, 
                                                         ordered_edge_map & hG, 
                                                         open_hash_map< NodeWeight > & node_weights) {
        forall_local_nodes(G, node) {
                NodeID cur_cnode = G.getCNode( node );
                node_weights.ref(cur_cnode) += G.getNodeWeight( node );

                forall_out_edges(G, e, node) {
                        NodeID target       = G.getEdgeTarget(e);
                        NodeID target_cnode = G.getCNode(target);
                        if( cur_cnode != target_cnode ) {
                                // update the edge
                                hG.weight(cur_cnode, target_cnode) += G.getEdgeWeight(e);
                        }
                } endfor
        } endfor
}



void parallel_contraction::redistribute_hased_graph_and_build_graph_locally( MPI_Comm communicator, ordered_edge_map & hG, 
                                                                             open_hash_map< NodeWeight > & node_weights,
                                                                             NodeID number_of_cnodes, 
                                                                             parallel_graph_access & Q  ) {
        PEID rank, size;
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);
        
        NodeID divisor          = ceil( number_of_cnodes/(double)size);

        //std::vector< std::vector< NodeID > >  messages;
        m_messages.resize(size);

        //build messages
        hG.for_each([&](NodeID source, NodeID target, NodeWeight weight) {
                PEID peID = source / divisor;
                m_messages[ peID ].push_back( source );
                m_messages[ peID ].push_back( target );
                m_messages[ peID ].push_back( weight );

                peID = target / divisor;
                m_messages[ peID ].push_back( target );
                m_messages[ peID ].push_back( source );
                m_messages[ peID ].push_back( weight );
        });

        // now flood the network
        std::vector<MPI_Request> requests(size, MPI_REQUEST_NULL);
        for( PEID peID = 0; peID < size; peID++) {
                if( peID != rank ) {
                        if( m_messages[peID].size() == 0 ){
                                m_messages[peID].push_back(std::numeric_limits<NodeID>::max());
                        }

                        MPI_Isend( &m_messages[peID][0], 
                                    m_messages[peID].size(), 
                                    MPI_UNSIGNED_LONG_LONG, 
                                    peID, peID+7*size, communicator, &requests[peID]);
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
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

        ordered_edge_map local_graph(number_of_cnodes);
        for( PEID peID = 0; peID < size; peID++) {
                if(local_msg_byPE[peID].size() > 0) {
                for( ULONG i = 0; i < local_msg_byPE[peID].size()-2; i+=3) {
                        local_graph.weight(local_msg_byPE[peID][i], local_msg_byPE[peID][i+1]) += local_msg_byPE[peID][i+2];
                }}
        }


        ULONG from = rank     * ceil(number_of_cnodes / (double)size);
        ULONG to   = (rank+1) * ceil(number_of_cnodes / (double)size) - 1;
        // handle the case where we dont have local edges
        from = std::min(from, number_of_cnodes);
        to   = std::min(to, number_of_cnodes - 1);
        ULONG local_num_cnodes = to - from + 1;

        // the adjacency of the local coarse nodes as one array: the edges are collected in iteration
        // order together with their local row, then a stable counting sort by row gives every row
        // its edges in the same order as appending them to one vector per row
        std::vector< NodeID > edge_row;
        std::vector< std::pair< NodeID, NodeWeight > > edge_list;
        edge_row.reserve(2 * local_graph.size());
        edge_list.reserve(2 * local_graph.size());

        local_graph.for_each([&](NodeID source, NodeID target, NodeWeight weight) {
                if( from <= target && target <= to) {
                        edge_row.push_back( target - from );
                        edge_list.push_back( std::make_pair( source, weight/4 ) );
                        edge_row.push_back( source - from );
                        edge_list.push_back( std::make_pair( target, weight/4 ) );
                } else {
                        edge_row.push_back( source - from );
                        edge_list.push_back( std::make_pair( target, weight/2 ) );
                }
        });
        EdgeID edge_counter = edge_list.size();

        std::vector< EdgeID > row_start( local_num_cnodes + 1, 0 );
        for( EdgeID i = 0; i < edge_counter; i++) row_start[ edge_row[i] + 1 ]++;
        for( NodeID r = 0; r < local_num_cnodes; r++) row_start[ r + 1 ] += row_start[ r ];
        std::vector< std::pair< NodeID, NodeWeight > > sorted_edges( edge_counter );
        {
                std::vector< EdgeID > fill( row_start.begin(), row_start.end() - 1 );
                for( EdgeID i = 0; i < edge_counter; i++) sorted_edges[ fill[ edge_row[i] ]++ ] = edge_list[i];
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
               
                for( EdgeID e = row_start[node]; e < row_start[node + 1]; e++) {
                        NodeID target = sorted_edges[e].first;
                        EdgeID e_bar = Q.new_edge(node, target);
                        Q.setEdgeWeight(e_bar, sorted_edges[e].second);
                }
        }

        Q.finish_construction();

        for( PEID peID = 0; peID < size; peID++) {
                m_messages[peID].clear();
        }
        //now distribute the node weights
        //pack messages
        // the order of the (node, weight) pairs is irrelevant: the receiver adds them up
        node_weights.for_each([&](NodeID node, NodeWeight weight) {
                PEID peID         = node / divisor;

                m_messages[ peID ].push_back( node );
                m_messages[ peID ].push_back( weight );
        });

        for( PEID peID = 0; peID < size; peID++) {
                if( peID != rank ) {
                        if( m_messages[peID].size() == 0 ){
                                m_messages[peID].push_back(std::numeric_limits<NodeID>::max());
                        }

                        MPI_Isend( &m_messages[peID][0], 
                                    m_messages[peID].size(), 
                                    MPI_UNSIGNED_LONG_LONG, 
                                    peID, peID+8*size, communicator, &requests[peID]);
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
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
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
        std::vector<MPI_Request> requests(m_send_buffers.size(), MPI_REQUEST_NULL);
        for( PEID peID = 0; peID < (PEID)m_send_buffers.size(); peID++) {
                if( G.is_adjacent_PE(peID) ) {
                        //now we have to send a message
                        if( m_send_buffers[peID].size() == 0 ){
                                // length 1 encode no message
                                m_send_buffers[peID].push_back(0);
                        }

                        MPI_Isend( &m_send_buffers[peID][0], 
                                    m_send_buffers[peID].size(), MPI_UNSIGNED_LONG_LONG, peID, peID+9*size, communicator, &requests[peID]);
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
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

}
