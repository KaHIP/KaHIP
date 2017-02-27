/******************************************************************************
 * hashed_graph.h
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

#ifndef HASHED_GRAPH_DG1JG7O0
#define HASHED_GRAPH_DG1JG7O0

#include <unordered_map>

#include "definitions.h"
#include "limits.h"

struct hashed_edge {
        NodeID k;
        NodeID source;
        NodeID target;
};


struct compare_hashed_edge {
        bool operator()(const hashed_edge e_1, const hashed_edge e_2) const {
                bool eq = (e_1.source == e_2.source && e_1.target == e_2.target);
                     eq = eq || (e_1.source == e_2.target && e_1.target == e_2.source); 
                return eq;
        }
};

struct data_hashed_edge{
        NodeWeight weight;

        data_hashed_edge() {
                weight = 0;
        }
};

struct hash_hashed_edge {
       ULONG operator()(const hashed_edge e) const {
                if(e.source < e.target) 
                        return e.source*e.k + e.target;
                else 
                        return e.target*e.k + e.source;
       }
};

typedef std::unordered_map<const hashed_edge, data_hashed_edge, hash_hashed_edge, compare_hashed_edge> hashed_graph;


#endif /* end of include guard: HASHED_GRAPH_DG1JG7O0 */
