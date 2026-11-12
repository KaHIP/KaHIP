/******************************************************************************
 * bk_maxflow.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 * Maximum flow by the algorithm of Boykov and Kolmogorov: a search tree grows
 * from the source and one from the sink, a path is augmented where they meet,
 * and the nodes cut off by saturated tree edges look for a new parent in their
 * tree (adoption, with timestamps and distances to the terminal). The result is
 * a flow (no excesses) with the usual partner convention of flow_graph: the flow
 * of an edge is the negated flow of its reverse edge.
 *
 * When solve() returns without stopping early, the source tree is exactly the
 * set of nodes reachable from the source in the residual graph: every tree node
 * has no residual edge to a free node (it would have been grown) and none to the
 * sink tree (that would be a path).
 *****************************************************************************/

#ifndef BK_MAXFLOW_H
#define BK_MAXFLOW_H

#include <algorithm>
#include <limits>
#include <vector>
#include "definitions.h"
#include "data_structure/flow_graph.h"

class bk_maxflow {
public:
        // Returns the value of a maximum flow, or the first flow value >= stop_flow.
        FlowType solve( flow_graph & G, NodeID s, NodeID t,
                        FlowType stop_flow = std::numeric_limits<FlowType>::max()) {
                m_G = &G;
                const NodeID n = G.number_of_nodes();
                TERMINAL = n; ORPHAN = n + 1; NONE = n + 2;
                m_st.assign(n, node_state{NONE, 0, 0, FREE, 0});
                m_pedge.resize(n);
                m_active.clear(); m_active_head = 0;
                m_orphans.clear();
                m_time = 0;

                m_st[s].tree = S; m_st[s].parent = TERMINAL;
                m_st[t].tree = T; m_st[t].parent = TERMINAL;
                set_active(s); set_active(t);

                FlowType flow = 0;
                while(true) {
                        rEdge * mid = nullptr; NodeID ps = 0, pt = 0;
                        grow(mid, ps, pt);
                        if(!mid) break;

                        m_time++;
                        long long b = res(*mid);
                        for(NodeID v = ps; m_st[v].parent != TERMINAL; v = m_st[v].parent) b = std::min(b, res(*m_pedge[v]));
                        for(NodeID v = pt; m_st[v].parent != TERMINAL; v = m_st[v].parent) b = std::min(b, res(*m_pedge[v]));

                        push(*mid, b);
                        for(NodeID v = ps; m_st[v].parent != TERMINAL; ) {
                                NodeID u = m_st[v].parent;
                                push(*m_pedge[v], b);
                                if(res(*m_pedge[v]) == 0) { m_st[v].parent = ORPHAN; m_orphans.push_back(v); }
                                v = u;
                        }
                        for(NodeID v = pt; m_st[v].parent != TERMINAL; ) {
                                NodeID u = m_st[v].parent;
                                push(*m_pedge[v], b);
                                if(res(*m_pedge[v]) == 0) { m_st[v].parent = ORPHAN; m_orphans.push_back(v); }
                                v = u;
                        }
                        flow += b;
                        if(flow >= stop_flow) return flow;

                        while(!m_orphans.empty()) {
                                NodeID p = m_orphans.back(); m_orphans.pop_back();
                                adopt(p);
                        }
                }
                return flow;
        }

        // after solve(): is v in the source tree?
        bool in_source_tree(NodeID v) const { return m_st[v].tree == S; }

private:
        enum : char { FREE = 0, S = 1, T = 2 };
        struct node_state {
                NodeID parent;
                int    ts;
                int    dist;
                char   tree;
                char   in_active;
        };

        static long long res(const rEdge & e) { return (long long)(NodeID)e.capacity - e.flow; }

        // add b to the flow of e and subtract it from the flow of its reverse edge
        void push(rEdge & e, long long b) {
                e.flow += b;
                m_G->edges_of(e.target)[e.reverse_edge_index].flow -= b;
        }

        void set_active(NodeID v) {
                if(!m_st[v].in_active) { m_st[v].in_active = 1; m_active.push_back(v); }
        }

        // Grows the trees from the active nodes in FIFO order until an edge from the source tree
        // to the sink tree is found (mid, from ps to pt) or no active node is left. The node that
        // found the edge stays active.
        void grow(rEdge * & mid, NodeID & ps, NodeID & pt) {
                flow_graph & G = *m_G;
                while(m_active_head < m_active.size()) {
                        const NodeID p = m_active[m_active_head];
                        node_state & sp = m_st[p];
                        if(sp.tree != FREE) {
                                rEdge * pe = G.edges_of(p);
                                const EdgeID end = G.get_first_invalid_edge(p);
                                if(sp.tree == S) {
                                        for(EdgeID e = 0; e < end; e++) {
                                                if(res(pe[e]) <= 0) continue;
                                                const NodeID q = pe[e].target;
                                                node_state & sq = m_st[q];
                                                if(sq.tree == FREE) {
                                                        sq.tree = S; sq.parent = p; m_pedge[q] = &pe[e];
                                                        sq.ts = sp.ts; sq.dist = sp.dist + 1;
                                                        set_active(q);
                                                } else if(sq.tree == T) {
                                                        mid = &pe[e]; ps = p; pt = q; return;
                                                }
                                        }
                                } else {
                                        for(EdgeID e = 0; e < end; e++) {
                                                const NodeID q = pe[e].target;
                                                node_state & sq = m_st[q];
                                                if(sq.tree == T) continue;
                                                rEdge & qp = G.edges_of(q)[pe[e].reverse_edge_index];
                                                if(res(qp) <= 0) continue;
                                                if(sq.tree == FREE) {
                                                        sq.tree = T; sq.parent = p; m_pedge[q] = &qp;
                                                        sq.ts = sp.ts; sq.dist = sp.dist + 1;
                                                        set_active(q);
                                                } else {
                                                        mid = &qp; ps = q; pt = p; return;
                                                }
                                        }
                                }
                        }
                        sp.in_active = 0; m_active_head++;
                }
        }

        // the edge that would make q the parent of p in tree tr (q -> p in S, p -> q in T)
        rEdge * tree_edge(rEdge & pq, char tr) {
                return tr == S ? &m_G->edges_of(pq.target)[pq.reverse_edge_index] : &pq;
        }

        void adopt(NodeID p) {
                const int INF = std::numeric_limits<int>::max();
                rEdge * pe = m_G->edges_of(p);
                const EdgeID end = m_G->get_first_invalid_edge(p);
                const char tr = m_st[p].tree;
                NodeID best = NONE; rEdge * best_edge = nullptr; int dmin = INF;
                for(EdgeID e = 0; e < end; e++) {
                        const NodeID q = pe[e].target;
                        if(m_st[q].tree != tr) continue;
                        rEdge * cand = tree_edge(pe[e], tr);
                        if(res(*cand) <= 0) continue;
                        // distance of q to the terminal, if q still has an origin
                        int d = 0; NodeID j = q;
                        while(true) {
                                if(m_st[j].ts == m_time) { d += m_st[j].dist; break; }
                                const NodeID pj = m_st[j].parent;
                                if(pj == TERMINAL) { m_st[j].ts = m_time; m_st[j].dist = 0; break; }
                                if(pj == ORPHAN) { d = INF; break; }
                                d++; j = pj;
                        }
                        if(d == INF) continue;
                        if(d < dmin) { dmin = d; best = q; best_edge = cand; }
                        int dd = d;
                        for(NodeID k = q; m_st[k].ts != m_time; k = m_st[k].parent) {
                                m_st[k].ts = m_time; m_st[k].dist = dd; dd--;
                        }
                }
                if(best != NONE) {
                        m_st[p].parent = best; m_pedge[p] = best_edge;
                        m_st[p].ts = m_time; m_st[p].dist = dmin + 1;
                        return;
                }
                for(EdgeID e = 0; e < end; e++) {
                        const NodeID q = pe[e].target;
                        if(m_st[q].tree != tr) continue;
                        if(res(*tree_edge(pe[e], tr)) > 0) set_active(q);
                        if(m_st[q].parent == p) { m_st[q].parent = ORPHAN; m_orphans.push_back(q); }
                }
                m_st[p].tree = FREE; m_st[p].parent = NONE;
        }

        flow_graph *            m_G;
        NodeID                  TERMINAL, ORPHAN, NONE;
        std::vector<node_state> m_st;
        std::vector<rEdge*>     m_pedge;
        std::vector<NodeID>     m_active;
        size_t                  m_active_head;
        std::vector<NodeID>     m_orphans;
        int                     m_time;
};

#endif
