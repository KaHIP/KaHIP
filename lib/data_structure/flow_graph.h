/******************************************************************************
 * flow_graph.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning 
 * 
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef FLOW_GRAPH_636S5L2S
#define FLOW_GRAPH_636S5L2S

#include "definitions.h"

struct rEdge {
    NodeID     source;
    NodeID     target;
    FlowType   capacity;
    FlowType   flow;
    EdgeID     reverse_edge_index;

    rEdge( NodeID source, NodeID target, FlowType capacity, FlowType flow, EdgeID reverse_edge_index) {
        this->source             = source;
        this->target             = target;
        this->capacity           = capacity;
        this->flow               = flow;
        this->reverse_edge_index = reverse_edge_index;
    }

};

// this is a adjacency list implementation of the residual graph
// for each edge we create, we create a rev edge with cap 0
// zero capacity edges are residual edges
class flow_graph {
public:
        flow_graph() {
                m_num_edges = 0;
                m_num_nodes = 0;
        };

        virtual ~flow_graph() {};

        void start_construction(NodeID nodes, EdgeID edges = 0) {
                m_adjacency_lists.resize(nodes);
                m_num_nodes = nodes;
                m_num_edges = edges;
        }
 
        void finish_construction() {};

        NodeID number_of_nodes() {return m_num_nodes;};
        EdgeID number_of_edges() {return m_num_edges;};

        NodeID getEdgeTarget(NodeID source, EdgeID e);
        NodeID getEdgeCapacity(NodeID source, EdgeID e);

        FlowType getEdgeFlow(NodeID source, EdgeID e);
        void setEdgeFlow(NodeID source, EdgeID e, FlowType flow);

        EdgeID getReverseEdge(NodeID source, EdgeID e);
        
        void new_edge(NodeID source, NodeID target, FlowType capacity) {
               m_adjacency_lists[source].push_back(rEdge(source, target, capacity, 0, m_adjacency_lists[target].size()));
               // for each edge we add a reverse edge
               m_adjacency_lists[target].push_back(rEdge(target, source, 0, 0, m_adjacency_lists[source].size() - 1));
               m_num_edges += 2;
        };

        EdgeID get_first_edge(NodeID node) {return 0;};
        EdgeID get_first_invalid_edge(NodeID node) {return m_adjacency_lists[node].size();};


private:
        std::vector< std::vector<rEdge> > m_adjacency_lists;
        NodeID m_num_nodes;
        NodeID m_num_edges;
};

inline
NodeID flow_graph::getEdgeCapacity(NodeID source, EdgeID e) {
#ifdef NDEBUG
        return m_adjacency_lists[source][e].capacity;        
#else
        return m_adjacency_lists.at(source).at(e).capacity;        
#endif
};

inline
void flow_graph::setEdgeFlow(NodeID source, EdgeID e, FlowType flow) {
#ifdef NDEBUG
        m_adjacency_lists[source][e].flow = flow;        
#else
        m_adjacency_lists.at(source).at(e).flow = flow;        
#endif
};

inline
FlowType flow_graph::getEdgeFlow(NodeID source, EdgeID e) {
#ifdef NDEBUG
        return m_adjacency_lists[source][e].flow;        
#else
        return m_adjacency_lists.at(source).at(e).flow;        
#endif
};

inline
NodeID flow_graph::getEdgeTarget(NodeID source, EdgeID e) {
#ifdef NDEBUG
        return m_adjacency_lists[source][e].target;        
#else
        return m_adjacency_lists.at(source).at(e).target;        
#endif
};

inline
EdgeID flow_graph::getReverseEdge(NodeID source, EdgeID e) {
#ifdef NDEBUG
        return m_adjacency_lists[source][e].reverse_edge_index;
#else
        return m_adjacency_lists.at(source).at(e).reverse_edge_index;        
#endif

}

#endif /* end of include guard: FLOW_GRAPH_636S5L2S */
