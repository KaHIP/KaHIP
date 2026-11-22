/******************************************************************************
 * parallel_contraction_send_lists.cpp
 * Reuse of the ghost cnode send lists across contractions of the same graph. Kept in its own translation unit,
 * so that the inlining in parallel_contraction.cpp does not change.
 *****************************************************************************/

#include "parallel_contraction.h"

// the adjacency of G is fixed after its construction, so the nodes sent to each PE and their order are the same
// on every contraction of G; only their cnodes change
bool parallel_contraction::refill_ghost_send_buffers( parallel_graph_access & G ) {
        std::vector< std::vector< NodeID > > & send_lists = G.ghost_send_lists();
        if( send_lists.empty() ) return false;

        const NodeID from = G.get_from_range();
        for( PEID peID = 0; peID < (PEID)send_lists.size(); peID++) {
                const std::vector< NodeID > & nodes = send_lists[peID];
                std::vector< NodeID > & buffer = m_send_buffers[peID];
                buffer.resize( 2 * nodes.size() );
                for( NodeID i = 0; i < nodes.size(); i++) {
                        buffer[2*i]   = from + nodes[i]; // global id of a local node
                        buffer[2*i+1] = G.getCNode(nodes[i]);
                }
        }
        return true;
}

void parallel_contraction::record_ghost_send_lists( parallel_graph_access & G ) {
        std::vector< std::vector< NodeID > > & send_lists = G.ghost_send_lists();
        const NodeID from = G.get_from_range();
        send_lists.resize( m_send_buffers.size() );
        for( PEID peID = 0; peID < (PEID)m_send_buffers.size(); peID++) {
                const std::vector< NodeID > & buffer = m_send_buffers[peID];
                send_lists[peID].resize( buffer.size() / 2 );
                for( NodeID i = 0; i < send_lists[peID].size(); i++) send_lists[peID][i] = buffer[2*i] - from;
        }
}
