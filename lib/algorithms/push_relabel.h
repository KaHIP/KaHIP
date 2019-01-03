/******************************************************************************
 * push_relabel.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef MAX_FLOW_MIN_CUT_Q5EJKHNS
#define MAX_FLOW_MIN_CUT_Q5EJKHNS

#include <iostream>
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
                m_excess.resize(G.number_of_nodes(),0);
                m_distance.resize(G.number_of_nodes(),0);
                m_active.resize(G.number_of_nodes(), false);
                m_count.resize(2*G.number_of_nodes(),0);
                m_bfstouched.resize(G.number_of_nodes());

                m_count[0] = G.number_of_nodes()-1;
                m_count[G.number_of_nodes()] = 1;
                
                m_distance[source] = G.number_of_nodes();
                m_active[source]   = true;
                m_active[sink]     = true;

                forall_out_edges(G, e, source) {
                        m_excess[source] += G.getEdgeCapacity(source, e);
                        push(source, e);
                } endfor
        
        }

        // perform a backward bfs in the residual starting at the sink
        // to update distance labels
        void global_relabeling( NodeID source, NodeID sink ) {
                std::queue< NodeID > Q; 
                forall_nodes((*m_G), node) {
                        m_distance[node]   = std::max(m_distance[node], m_G->number_of_nodes());
                        m_bfstouched[node] = false;
                } endfor
                
                Q.push(sink);
                m_bfstouched[sink]   = true;
                m_bfstouched[source] = true;
                m_distance[sink]     = 0;

                NodeID node = 0;
                while( !Q.empty() ) {
                       node = Q.front();
                       Q.pop();  

                       forall_out_edges((*m_G), e, node) {
                               NodeID target = m_G->getEdgeTarget(node, e);
                               if(m_bfstouched[target]) continue;

                               EdgeID rev_e = m_G->getReverseEdge(node, e);
                               if( m_G->getEdgeCapacity( target, rev_e) - m_G->getEdgeFlow( target, rev_e) > 0 ) {
                                        m_count[ m_distance[target] ] --;
                                        m_distance[target] = m_distance[node]+1;
                                        m_count[ m_distance[target] ] ++;
                                        Q.push(target);
                                        m_bfstouched[target] = true;
                               }
                       } endfor
                }
        }

        // push flow from source to target if possible
        void push( NodeID source, EdgeID e) {
                m_pushes++;
                FlowType capacity = m_G->getEdgeCapacity(source, e);
                FlowType flow     = m_G->getEdgeFlow(source, e);
                FlowType amount   = std::min((long long)(capacity-flow), m_excess[source]);
                NodeID target     = m_G->getEdgeTarget(source, e);

                if( m_distance[source] <= m_distance[target] || amount == 0) return;

                m_G->setEdgeFlow(source, e, flow+amount);
 
                EdgeID   rev_e    = m_G->getReverseEdge(source, e);
                FlowType rev_flow = m_G->getEdgeFlow(target, rev_e);
                m_G->setEdgeFlow(target, rev_e, rev_flow-amount);

                m_excess[source] -= amount;
                m_excess[target] += amount;

                enqueue(target);
        }

        // put a vertex in the FIFO queue of the algorithm 
        void enqueue( NodeID target ) {
                if( m_active[target] ) return;
                if( m_excess[target] > 0) {
                        m_active[target] = true;
                        //m_Q.push(target, m_distance[target]);
                        m_Q.push(target);
                }
        }

        // try to push as much excess as possible out of the node node
        void discharge( NodeID node ) {
                EdgeID end = m_G->get_first_invalid_edge(node);
                for(EdgeID e = m_G->get_first_edge(node); e < end && m_excess[node] > 0; ++e) {
                        push( node, e );
                }

                if( m_excess[node] > 0 ) {
                        if( m_count[ m_distance[node] ] == 1 && m_distance[node] < m_G->number_of_nodes()) {
                                // hence this layer will be empty after the relabel step
                                gap_heuristic(m_distance[node]);
                        } else {
                                relabel(node);
                        }
                }
        }


        // gap heuristic
        void gap_heuristic( NodeID level ) {
                m_gaps++;
                forall_nodes((*m_G), node) {
                        if(m_distance[node] < level ) continue;
                        m_count[m_distance[node]]--;
                        m_distance[node] = std::max(m_distance[node], m_G->number_of_nodes());
                        m_count[m_distance[node]]++;
                        enqueue(node);
                } endfor
        }

        // relabel a node with respect to its 
        // neighboring nodes
        void relabel( NodeID node ) {
                m_work += WORK_OP_RELABEL;
                m_num_relabels++;

                m_count[m_distance[node]]--;
                m_distance[node] = 2*m_G->number_of_nodes();

                forall_out_edges((*m_G), e, node) {
                        if( m_G->getEdgeCapacity( node, e) - m_G->getEdgeFlow( node, e) > 0) {
                                NodeID target    = m_G->getEdgeTarget( node, e);
                                m_distance[node] = std::min(m_distance[node], m_distance[target]+1);
                        }
                        m_work++;
                } endfor

                m_count[m_distance[node]]++;
                enqueue(node);
        }

        FlowType solve_max_flow_min_cut( flow_graph & G, 
                                         NodeID source, 
                                         NodeID sink, 
                                         bool compute_source_set, 
                                         std::vector< NodeID > & source_set) {
                m_G                  = & G;
                m_work               = 0;
                m_num_relabels       = 0;
                m_gaps               = 0;
                m_pushes             = 0;
                m_global_updates     = 1;

                init(G, source, sink);
                global_relabeling( source, sink );
         
                int work_todo = WORK_NODE_TO_EDGES*G.number_of_nodes() + G.number_of_edges();
                // main loop
                while(!m_Q.empty()) {
                        NodeID v    = m_Q.front(); m_Q.pop();
                        m_active[v] = false;
                        discharge(v);

                        if( m_work > GLOBAL_UPDATE_FRQ*work_todo) {
                                global_relabeling( source,  sink );
                                m_work = 0;
                                m_global_updates++;
                        }
                }

                if(compute_source_set) {
                        // perform bfs starting from source set 
                        source_set.clear();

                        forall_nodes(G, node) {
                                m_bfstouched[node] = false;
                        } endfor

                        std::queue< NodeID > Q;
                        Q.push(source);
                        m_bfstouched[source] = true;

                        while( !Q.empty() ) {
                                NodeID node = Q.front();
                                Q.pop();
                                source_set.push_back(node);

                                forall_out_edges(G, e, node) {
                                        NodeID target = G.getEdgeTarget(node, e);
                                        FlowType resCap = G.getEdgeCapacity(node, e) - G.getEdgeFlow(node, e);
                                        if(resCap > 0 && !m_bfstouched[target]) {
                                                Q.push(target);
                                                m_bfstouched[target] = true;
                                        }
                                } endfor
                        }
                }

                //return value of flow
                return m_excess[sink];
        }
private:
        std::vector<long long> m_excess;
        std::vector<NodeID>    m_distance;
	std::vector<bool>      m_active; // store which nodes are in the queue already
	std::vector<int>       m_count;
        std::queue<NodeID>     m_Q;
	std::vector<bool>      m_bfstouched; 
        //highest_label_queue    m_Q;
        int m_num_relabels;
        int m_gaps;
        int m_global_updates;
        int m_pushes;
        int m_work;
        flow_graph * m_G;
};


#endif /* end of include guard: MAX_FLOW_MIN_CUT_Q5EJKHNS */
