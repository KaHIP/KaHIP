/******************************************************************************
 * partial_boundary.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PARTIAL_BOUNDARY_963CRO9F_
#define PARTIAL_BOUNDARY_963CRO9F_

#include "definitions.h"

struct compare_nodes_contains {
        bool operator()(const NodeID& lhs, const NodeID& rhs) const {
                return lhs == rhs;
        }
};

struct is_boundary {
        bool contains;
        is_boundary()
        : contains(false) {
        }
};

struct hash_boundary_nodes {
       size_t operator()(const NodeID idx) const {
                return idx;
       }
};

typedef extlib::unordered_map_with_custom_hash_and_comparator<NodeID, is_boundary, hash_boundary_nodes, compare_nodes_contains> is_boundary_node_hashtable;

class PartialBoundary {
        public:
                PartialBoundary( );
                virtual ~PartialBoundary();

                bool contains(NodeID node);
                void insert(NodeID node);
                void deleteNode(NodeID node);
                NodeID size();

                is_boundary_node_hashtable internal_boundary;
};

inline bool PartialBoundary::contains(NodeID node) {
        return internal_boundary.find(node) != internal_boundary.end(); 
}

inline void PartialBoundary::insert(NodeID node) {
        internal_boundary[node].contains = true; 
}

inline void PartialBoundary::deleteNode(NodeID node) {
        internal_boundary.erase(node);
}

inline NodeID PartialBoundary::size() {
        return internal_boundary.size();
}


//iterator for
#define forall_boundary_nodes(boundary, n) { is_boundary_node_hashtable::iterator iter; NodeID n; for(iter = boundary.internal_boundary.begin(); iter != boundary.internal_boundary.end(); iter++ ) { n = iter->first;

#endif /* end of include guard: PARTIAL_BOUNDARY_963CRO9F */
