/******************************************************************************
 * push_relabel.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef MAX_FLOW_MIN_CUT_Q5EJKHNS
#define MAX_FLOW_MIN_CUT_Q5EJKHNS

#include <iostream>
#include <limits>
#include "definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/flow_graph.h"
#include "tools/timer.h"

const int    WORK_OP_RELABEL    = 9;
const double GLOBAL_UPDATE_FRQ  = 0.51;
const int    WORK_NODE_TO_EDGES = 4;

class push_relabel {
public:
        push_relabel( );
        virtual ~push_relabel();

        void init( flow_graph & G, NodeID source, NodeID sink ) {
                m_st.resize(G.number_of_nodes(), node_state());
                m_count.resize(2*G.number_of_nodes()+2,0); // labels go up to 2n
                // FIFO of active nodes as a ring (a node is queued at most once at a time)
                m_ring.assign(std::max<NodeID>(G.number_of_nodes(), 1), 0);
                m_ring_head  = 0;
                m_ring_count = 0;

                m_count[0] = G.number_of_nodes()-1;
                m_count[G.number_of_nodes()] = 1;
                
                m_st[source].distance = G.number_of_nodes();
                m_st[source].active   = true;
                m_st[sink].active     = true;

                forall_out_edges(G, e, source) {
                        m_st[source].excess += G.getEdgeCapacity(source, e);
                        push(source, e);
                } endfor
        
        }

        // perform a backward bfs in the residual starting at the sink
        // to update distance labels
        void global_relabeling( NodeID source, NodeID sink ) {
                global_relabeling_from(sink, source, 0);
        }

        // labels: base + residual distance to root for the nodes that reach root (other excluded),
        // at least base + n for all others
        void global_relabeling_from( NodeID root, NodeID other, NodeID base ) {
                // FIFO in an array (each node is queued at most once), same order as std::queue
                const NodeID n = m_G->number_of_nodes();
                m_bfsqueue.resize(n);
                NodeID * Q = m_bfsqueue.data();
                size_t q_head = 0, q_tail = 0;
                forall_nodes((*m_G), node) {
                        m_st[node].distance   = std::max(m_st[node].distance, base + n);
                        m_st[node].touched = false;
                } endfor
                
                Q[q_tail++] = root;
                m_st[root].touched   = true;
                m_st[other].touched  = true;
                m_st[root].distance  = base;

                NodeID node = 0;
                while( q_head < q_tail ) {
                       node = Q[q_head++];

                       const rEdge * node_edges = m_G->edges_of(node);
                       const EdgeID  end        = m_G->get_first_invalid_edge(node);
                       for(EdgeID e = 0; e < end; ++e) {
                               NodeID target = node_edges[e].target;
                               if(m_st[target].touched) continue;

                               const rEdge & rev = m_G->edges_of(target)[node_edges[e].reverse_edge_index];
                               // capacities are read as NodeID, as getEdgeCapacity returns them
                               if( (NodeID) rev.capacity - rev.flow > 0 ) {
                                        m_count[ m_st[target].distance ] --;
                                        m_st[target].distance = m_st[node].distance+1;
                                        m_count[ m_st[target].distance ] ++;
                                        Q[q_tail++] = target;
                                        m_st[target].touched = true;
                               }
                       }
                }
        }

        // push flow from source to target if possible
        void push( NodeID source, EdgeID e) {
                m_pushes++;
                rEdge & edge      = m_G->edges_of(source)[e];
                FlowType capacity = (NodeID) edge.capacity; // as returned by getEdgeCapacity
                FlowType flow     = edge.flow;
                FlowType amount   = std::min((long long)(capacity-flow), m_st[source].excess);
                NodeID target     = edge.target;

                if( m_st[source].distance <= m_st[target].distance || amount == 0) return;

                edge.flow = flow+amount;
 
                rEdge & rev       = m_G->edges_of(target)[edge.reverse_edge_index];
                FlowType rev_flow = rev.flow;
                rev.flow          = rev_flow-amount;

                m_st[source].excess -= amount;
                m_st[target].excess += amount;

                enqueue(target);
        }

        // put a vertex in the FIFO queue of the algorithm 
        void enqueue( NodeID target ) {
                if( m_st[target].active ) return;
                if( m_st[target].excess > 0) {
                        m_st[target].active = true;
                        size_t pos = m_ring_head + m_ring_count;
                        if(pos >= m_ring.size()) pos -= m_ring.size();
                        m_ring[pos] = target;
                        m_ring_count++;
                }
        }

        // try to push as much excess as possible out of the node node
        void discharge( NodeID node ) {
                // the pushes of push(node, e) with the node's distance (unchanged by pushes) and
                // excess held in locals; a push never goes to node itself (equal distance)
                EdgeID end = m_G->get_first_invalid_edge(node);
                rEdge * node_edges     = m_G->edges_of(node);
                const NodeID dist_node = m_st[node].distance;
                long long excess       = m_st[node].excess;
                for(EdgeID e = m_G->get_first_edge(node); e < end && excess > 0; ++e) {
                        m_pushes++;
                        rEdge & edge      = node_edges[e];
                        FlowType capacity = (NodeID) edge.capacity; // as returned by getEdgeCapacity
                        FlowType flow     = edge.flow;
                        FlowType amount   = std::min((long long)(capacity-flow), excess);
                        NodeID target     = edge.target;

                        if( dist_node <= m_st[target].distance || amount == 0) continue;

                        edge.flow = flow+amount;
                        rEdge & rev = m_G->edges_of(target)[edge.reverse_edge_index];
                        rev.flow    = rev.flow-amount;

                        excess            -= amount;
                        m_st[target].excess  += amount;

                        enqueue(target);
                }
                m_st[node].excess = excess;

                if( m_st[node].excess > 0 ) {
                        if( m_count[ m_st[node].distance ] == 1 && m_st[node].distance < m_G->number_of_nodes()) {
                                // hence this layer will be empty after the relabel step
                                gap_heuristic(m_st[node].distance);
                        } else {
                                relabel(node);
                        }
                }
        }


        // gap heuristic
        void gap_heuristic( NodeID level ) {
                m_gaps++;
                forall_nodes((*m_G), node) {
                        if(m_st[node].distance < level ) continue;
                        m_count[m_st[node].distance]--;
                        m_st[node].distance = std::max(m_st[node].distance, m_G->number_of_nodes());
                        m_count[m_st[node].distance]++;
                        enqueue(node);
                } endfor
        }

        // relabel a node with respect to its 
        // neighboring nodes
        void relabel( NodeID node ) {
                m_work += WORK_OP_RELABEL;
                m_num_relabels++;

                m_count[m_st[node].distance]--;
                m_st[node].distance = 2*m_G->number_of_nodes();

                const rEdge * node_edges = m_G->edges_of(node);
                const EdgeID  end        = m_G->get_first_invalid_edge(node);
                for(EdgeID e = 0; e < end; ++e) {
                        if( (NodeID) node_edges[e].capacity - node_edges[e].flow > 0) {
                                NodeID target    = node_edges[e].target;
                                m_st[node].distance = std::min(m_st[node].distance, m_st[target].distance+1);
                        }
                        m_work++;
                }

                m_count[m_st[node].distance]++;
                enqueue(node);
        }

        FlowType solve_max_flow_min_cut( flow_graph & G, 
                                         NodeID source, 
                                         NodeID sink, 
                                         bool compute_source_set, 
                                         std::vector< NodeID > & source_set,
                                         FlowType stop_flow = std::numeric_limits<FlowType>::max()) {
                // stop_flow: an upper bound on the maximum flow known to the caller; once the sink
                // has received that much the flow value is final and the solver returns at once
                // (source_set is then left empty)
                m_G                  = & G;
                m_work               = 0;
                m_num_relabels       = 0;
                m_gaps               = 0;
                m_pushes             = 0;
                m_global_updates     = 1;

                init(G, source, sink);
                global_relabeling( source, sink );
         
                int work_todo = WORK_NODE_TO_EDGES*G.number_of_nodes() + G.number_of_edges();
                const NodeID n = G.number_of_nodes();
                // phase 1: maximum preflow. Only nodes with label < n are discharged; a label >= n
                // means that the node cannot reach the sink in the residual graph (the labels are
                // valid and the gap heuristic fires only on true gaps), so its excess cannot reach it
                // either and stays where it is until phase 2.
                while(m_ring_count > 0) {
                        NodeID v    = m_ring[m_ring_head];
                        if(++m_ring_head == m_ring.size()) m_ring_head = 0;
                        m_ring_count--;
                        m_st[v].active = false;
                        if( m_st[v].distance >= n ) continue;
                        discharge(v);
                        if( m_st[sink].excess >= stop_flow ) return m_st[sink].excess;

                        if( m_work > GLOBAL_UPDATE_FRQ*work_todo) {
                                global_relabeling( source,  sink );
                                m_work = 0;
                                m_global_updates++;
                        }
                }

                // phase 2: return the remaining excess to the source, with labels n + residual
                // distance to the source (the sink gets 2n and is never pushed to). The flow value
                // is final after phase 1; the result is a maximum flow.
                global_relabeling_from(source, sink, n);
                m_work = 0;
                forall_nodes(G, node) {
                        if(node != source && node != sink) enqueue(node);
                } endfor
                while(m_ring_count > 0) {
                        NodeID v    = m_ring[m_ring_head];
                        if(++m_ring_head == m_ring.size()) m_ring_head = 0;
                        m_ring_count--;
                        m_st[v].active = false;
                        discharge(v);

                        if( m_work > GLOBAL_UPDATE_FRQ*work_todo) {
                                global_relabeling_from(source, sink, n);
                                m_work = 0;
                        }
                }

                if(compute_source_set) {
                        // perform bfs starting from source set 
                        source_set.clear();

                        forall_nodes(G, node) {
                                m_st[node].touched = false;
                        } endfor

                        std::queue< NodeID > Q;
                        Q.push(source);
                        m_st[source].touched = true;

                        while( !Q.empty() ) {
                                NodeID node = Q.front();
                                Q.pop();
                                source_set.push_back(node);

                                forall_out_edges(G, e, node) {
                                        NodeID target = G.getEdgeTarget(node, e);
                                        FlowType resCap = G.getEdgeCapacity(node, e) - G.getEdgeFlow(node, e);
                                        if(resCap > 0 && !m_st[target].touched) {
                                                Q.push(target);
                                                m_st[target].touched = true;
                                        }
                                } endfor
                        }
                }

                //return value of flow
                return m_st[sink].excess;
        }
private:
        // per-node state in one record, so that a push or BFS step reads one cache line per target
        struct node_state {
                long long excess   = 0;
                NodeID    distance = 0;
                char      active   = 0; // node is in the queue already
                char      touched  = 0; // BFS mark
        };
        std::vector<node_state> m_st;
	std::vector<int>       m_count;
        std::vector<NodeID>    m_ring;
        size_t                 m_ring_head  = 0;
        size_t                 m_ring_count = 0;
        std::vector<NodeID>    m_bfsqueue;
        //highest_label_queue    m_Q;
        int m_num_relabels;
        int m_gaps;
        int m_global_updates;
        int m_pushes;
        int m_work;
        flow_graph * m_G;
};


#endif /* end of include guard: MAX_FLOW_MIN_CUT_Q5EJKHNS */
