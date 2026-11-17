/******************************************************************************
 * parallel_projection.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <algorithm>

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

        // nodes whose coarse node lives on another PE, as (cnode, node) in node order
        std::vector< std::pair< NodeID, NodeID > > remote;
        forall_local_nodes(finer, node) {
                NodeID cnode = finer.getCNode(node);
                if( coarser.is_local_node_from_global_id(cnode) ) {
                        NodeID new_label = coarser.getNodeLabel(coarser.getLocalID(cnode));
                        finer.setNodeLabel(node, new_label);
                } else {
                        remote.push_back(std::make_pair(cnode, node));
                }
        } endfor

        // group by cnode (nodes stay in increasing order within a group); request every cnode
        // once, in the order of its first node, as the map based version did
        std::stable_sort(remote.begin(), remote.end(),
                         [](const std::pair< NodeID, NodeID > & a, const std::pair< NodeID, NodeID > & b) {
                                return a.first < b.first; });
        std::vector< size_t > group_start;
        for( size_t i = 0; i < remote.size(); i++) {
                if( i == 0 || remote[i].first != remote[i-1].first ) group_start.push_back(i);
        }
        group_start.push_back(remote.size());
        std::vector< size_t > group_order(group_start.size() - 1);
        for( size_t g = 0; g < group_order.size(); g++) group_order[g] = g;
        std::sort(group_order.begin(), group_order.end(), [&](size_t a, size_t b) {
                return remote[group_start[a]].second < remote[group_start[b]].second; });

        std::vector< std::vector< size_t > > groups_of_pe(size);
        for( size_t gi = 0; gi < group_order.size(); gi++) {
                size_t g     = group_order[gi];
                NodeID cnode = remote[group_start[g]].first;
                PEID peID    = cnode / divisor;
                m_messages[peID].push_back(cnode); // we are requesting the label of this node 
                groups_of_pe[peID].push_back(g);
        }

        std::vector<MPI_Request> requests(size, MPI_REQUEST_NULL);
        for( PEID peID = 0; peID < size; peID++) {
                if( peID != rank ) {
                        if( m_messages[peID].size() == 0 ){
                                m_messages[peID].push_back(std::numeric_limits<NodeID>::max());
                        }

                        MPI_Isend( &m_messages[peID][0], 
                                                m_messages[peID].size(), 
                                                MPI_UNSIGNED_LONG_LONG, 
                                                peID, peID+size, communicator, &requests[peID]);
                }
        }

        std::vector< std::vector< NodeID > > out_messages;
        out_messages.resize(size);
        std::vector<MPI_Request> replies(size, MPI_REQUEST_NULL);

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
                        MPI_Isend( &out_messages[peID][0], 
                                        out_messages[peID].size(), 
                                        MPI_UNSIGNED_LONG_LONG, 
                                        peID, peID+2*size, communicator, &replies[peID]);

                        continue; // nothing to do
                }


                for( int i = 0; i < message_length; i++) {
                        NodeID cnode = coarser.getLocalID(incmessage[i]);
                        out_messages[peID].push_back(coarser.getNodeLabel(cnode));
                }

                MPI_Isend( &out_messages[peID][0], 
                                out_messages[peID].size(), 
                                MPI_UNSIGNED_LONG_LONG, 
                                peID, peID+2*size, communicator, &replies[peID]);

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
                        size_t g     = groups_of_pe[peID][i];
                        NodeID label = incmessage[i];

                        for( size_t j = group_start[g]; j < group_start[g+1]; j++) {
                                finer.setNodeLabel(remote[j].second, label);
                        }
                }
        }
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
        MPI_Waitall(replies.size(), replies.data(), MPI_STATUSES_IGNORE);

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
