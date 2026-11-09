/******************************************************************************
 * flow_graph.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef FLOW_GRAPH_636S5L2S
#define FLOW_GRAPH_636S5L2S

#include <vector>
#include "definitions.h"

struct rEdge {
    NodeID     target;
    EdgeID     reverse_edge_index;
    FlowType   capacity;
    FlowType   flow;
};

// Residual graph in CSR form. new_edge() records the edge (and its reverse edge with capacity 0);
// finish_construction() lays the edges out per node in insertion order, so the local edge indices and
// reverse indices are exactly those of the former adjacency-list implementation (one std::vector per
// node, appended in insertion order). finish_construction() must be called before the graph is used.
class flow_graph {
public:
        flow_graph() {
                m_num_edges = 0;
                m_num_nodes = 0;
        };
        virtual ~flow_graph() {};

        void start_construction(NodeID nodes, EdgeID edges = 0) {
                m_num_nodes = nodes;
                m_num_edges = edges;
                m_offset.assign(nodes + 1, 0);
                m_pending.clear();
        }

        void finish_construction() {
                // m_offset[v+1] holds the degree of v; turn it into offsets
                for(NodeID v = 0; v < m_num_nodes; v++) m_offset[v+1] += m_offset[v];
                m_edges.resize(m_offset[m_num_nodes]);
                std::vector<EdgeID> fill(m_num_nodes, 0);
                for(size_t i = 0; i < m_pending.size(); i++) {
                        const pending_edge & pe = m_pending[i];
                        EdgeID ls = fill[pe.source]++;
                        EdgeID lt = fill[pe.target]++;
                        rEdge & f = m_edges[m_offset[pe.source] + ls];
                        f.target = pe.target; f.reverse_edge_index = lt; f.capacity = pe.capacity; f.flow = 0;
                        rEdge & r = m_edges[m_offset[pe.target] + lt];
                        r.target = pe.source; r.reverse_edge_index = ls; r.capacity = 0; r.flow = 0;
                }
                std::vector<pending_edge>().swap(m_pending);
        };

        NodeID number_of_nodes() {return m_num_nodes;};
        EdgeID number_of_edges() {return m_num_edges;};

        NodeID getEdgeTarget(NodeID source, EdgeID e);
        NodeID getEdgeCapacity(NodeID source, EdgeID e);
        FlowType getEdgeFlow(NodeID source, EdgeID e);
        void setEdgeFlow(NodeID source, EdgeID e, FlowType flow);
        EdgeID getReverseEdge(NodeID source, EdgeID e);

        void new_edge(NodeID source, NodeID target, FlowType capacity) {
               pending_edge pe; pe.source = source; pe.target = target; pe.capacity = capacity;
               m_pending.push_back(pe);
               m_offset[source+1]++;
               m_offset[target+1]++;
               m_num_edges += 2;
        };

        EdgeID get_first_edge(NodeID node) {return 0;};
        EdgeID get_first_invalid_edge(NodeID node) {return m_offset[node+1] - m_offset[node];};

private:
        struct pending_edge {
                NodeID source;
                NodeID target;
                FlowType capacity;
        };
        std::vector<EdgeID> m_offset;
        std::vector<rEdge> m_edges;
        std::vector<pending_edge> m_pending;
        NodeID m_num_nodes;
        NodeID m_num_edges;
};

inline
NodeID flow_graph::getEdgeCapacity(NodeID source, EdgeID e) {
        return m_edges[m_offset[source] + e].capacity;
};

inline
void flow_graph::setEdgeFlow(NodeID source, EdgeID e, FlowType flow) {
        m_edges[m_offset[source] + e].flow = flow;
};

inline
FlowType flow_graph::getEdgeFlow(NodeID source, EdgeID e) {
        return m_edges[m_offset[source] + e].flow;
};

inline
NodeID flow_graph::getEdgeTarget(NodeID source, EdgeID e) {
        return m_edges[m_offset[source] + e].target;
};

inline
EdgeID flow_graph::getReverseEdge(NodeID source, EdgeID e) {
        return m_edges[m_offset[source] + e].reverse_edge_index;
}

#endif /* end of include guard: FLOW_GRAPH_636S5L2S */
