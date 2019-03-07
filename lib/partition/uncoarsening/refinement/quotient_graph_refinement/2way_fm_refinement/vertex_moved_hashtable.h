/******************************************************************************
 * vertex_moved_hashtable.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef VMOVEDHT_4563r97820954
#define VMOVEDHT_4563r97820954

#include <unordered_map>

#include "definitions.h"
#include "limits.h"

struct compare_nodes {
        bool operator()(const NodeID lhs, const NodeID rhs) const {
                return (lhs == rhs);
        }
};


const NodeID NOT_MOVED = std::numeric_limits<NodeID>::max();
const NodeID MOVED = 0;

struct moved_index {
       NodeID index;
       moved_index() {
                index = NOT_MOVED;
       }
};

struct hash_nodes {
       size_t operator()(const NodeID idx) const {
                return idx;
       }
};

typedef std::unordered_map<const NodeID, moved_index, hash_nodes, compare_nodes> vertex_moved_hashtable;

#endif
