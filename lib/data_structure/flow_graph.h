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

// Residual graph in CSR form. new_edge() records the edge and its reverse edge (capacity rev_capacity, 0
// by default);
// finish_construction() lays the edges out per node in insertion order, so the local edge indices and
// reverse indices are exactly those of the former adjacency-list implementation (one std::vector per
// node, appended in insertion order). finish_construction() must be called before the graph is used.
class flow_graph {
public:
        flow_graph() {
                m_num_edges = 0;
                m_num_nodes = 0;
        };
        // the edge arrays go back to the spare buffers, so the next flow problem reuses their capacity
        virtual ~flow_graph() {
                if(spare_edges().capacity() < m_edges.capacity()) spare_edges().swap(m_edges);
                if(spare_pending().capacity() < m_pending.capacity()) spare_pending().swap(m_pending);
        };

        void start_construction(NodeID nodes, EdgeID edges = 0) {
                m_num_nodes = nodes;
                m_num_edges = edges;
                m_offset.assign(nodes + 1, 0);
                if(m_pending.capacity() == 0) m_pending.swap(spare_pending());
                m_pending.clear();
                if(m_edges.capacity() == 0) m_edges.swap(spare_edges());
                m_edges.clear();
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
                        r.target = pe.source; r.reverse_edge_index = ls; r.capacity = pe.rev_capacity; r.flow = 0;
                }
                m_pending.clear();
                if(spare_pending().capacity() < m_pending.capacity()) spare_pending().swap(m_pending);
        };

        NodeID number_of_nodes() {return m_num_nodes;};
        EdgeID number_of_edges() {return m_num_edges;};

        NodeID getEdgeTarget(NodeID source, EdgeID e);
        NodeID getEdgeCapacity(NodeID source, EdgeID e);
        FlowType getEdgeFlow(NodeID source, EdgeID e);
        void setEdgeFlow(NodeID source, EdgeID e, FlowType flow);
        EdgeID getReverseEdge(NodeID source, EdgeID e);

        void new_edge(NodeID source, NodeID target, FlowType capacity, FlowType rev_capacity = 0) {
               pending_edge pe; pe.source = source; pe.target = target; pe.capacity = capacity;
               pe.rev_capacity = rev_capacity;
               m_pending.push_back(pe);
               m_offset[source+1]++;
               m_offset[target+1]++;
               m_num_edges += 2;
        };

        // the edges of node (local index e is edges_of(node)[e])
        rEdge * edges_of(NodeID node) {return m_edges.data() + m_offset[node];};

        EdgeID get_first_edge(NodeID node) {return 0;};
        EdgeID get_first_invalid_edge(NodeID node) {return m_offset[node+1] - m_offset[node];};

private:
        struct pending_edge {
                NodeID source;
                NodeID target;
                FlowType capacity;
                FlowType rev_capacity;
        };
        // buffers kept between flow problems (taken by start_construction, returned when done)
        static std::vector<pending_edge> & spare_pending() { static std::vector<pending_edge> b; return b; }
        static std::vector<rEdge> & spare_edges() { static std::vector<rEdge> b; return b; }

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
