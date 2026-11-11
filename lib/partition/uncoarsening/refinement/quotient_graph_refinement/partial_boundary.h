/******************************************************************************
 * partial_boundary.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PARTIAL_BOUNDARY_963CRO9F_
#define PARTIAL_BOUNDARY_963CRO9F_

#include <unordered_map>
#include <vector>
#include "definitions.h"
#include "data_structure/node_hash_map.h"

struct compare_nodes_contains {
        bool operator()(const NodeID lhs, const NodeID rhs) const {
                return (lhs == rhs);
        }
};


struct is_boundary {
       bool contains;
       is_boundary() {
                contains = false;
       }
};


struct hash_boundary_nodes {
       size_t operator()(const NodeID idx) const {
                return idx;
       }
};

// Boundary nodes in a dense array with an open-addressing position index: insert, delete and
// contains are O(1) without node allocation; deletion swaps the last node into the freed slot.
// Iteration runs over the array (rule (b): the order differs from the former std::unordered_map).
class PartialBoundary {
        public:
                PartialBoundary( );
                virtual ~PartialBoundary();

                bool contains(NodeID node);
                void insert(NodeID node);
                void clear();
                void deleteNode(NodeID node);
                NodeID size();

                std::vector<NodeID>   m_nodes;
                node_hash_map<NodeID> m_pos;
};

inline bool PartialBoundary::contains(NodeID node) {
        return m_pos.find(node) != m_pos.end();
}

inline void PartialBoundary::insert(NodeID node) {
        if(m_pos.find(node) != m_pos.end()) return;
        m_pos[node] = m_nodes.size();
        m_nodes.push_back(node);
}

inline void PartialBoundary::deleteNode(NodeID node) {
        NodeID * p = m_pos.find(node);
        if(p == m_pos.end()) return;
        NodeID pos  = *p;
        NodeID last = m_nodes.back();
        m_nodes[pos] = last;
        m_nodes.pop_back();
        if(last != node) *m_pos.find(last) = pos;
        m_pos.erase(node);
}

inline NodeID PartialBoundary::size() {
        return m_nodes.size();
}

inline void PartialBoundary::clear() {
        m_nodes.clear();
        m_pos = node_hash_map<NodeID>();
}



//iterator for
#define forall_boundary_nodes(boundary, n) { NodeID n; for(size_t bnd_idx_ = 0; bnd_idx_ < boundary.m_nodes.size(); bnd_idx_++ ) { n = boundary.m_nodes[bnd_idx_];

#endif /* end of include guard: PARTIAL_BOUNDARY_963CRO9F */
