/******************************************************************************
 * vertex_moved_hashtable.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Graph Partitioning 
 ******************************************************************************
 * Copyright (C) 2017 Christian Schulz 
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

#ifndef VMOVEDHT_4563r97820954
#define VMOVEDHT_4563r97820954

#include <unordered_map>

#include "definitions.h"
#include "limits.h"

using namespace __gnu_cxx;

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
