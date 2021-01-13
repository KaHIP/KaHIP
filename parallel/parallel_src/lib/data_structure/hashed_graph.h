/******************************************************************************
 * hashed_graph.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef HASHED_GRAPH_DG1JG7O0
#define HASHED_GRAPH_DG1JG7O0

#include "definitions.h"
#include "limits.h"

struct hashed_edge {
        NodeID k;
        NodeID source;
        NodeID target;
};

struct compare_hashed_edge {
        bool operator()(const hashed_edge& e_1, const hashed_edge& e_2) const {
                return (e_1.source == e_2.source && e_1.target == e_2.target) ||
                       (e_1.source == e_2.target && e_1.target == e_2.source);
        }
};

struct data_hashed_edge {
        NodeWeight weight;

        data_hashed_edge()
        : weight(0) {
        }
};

struct hash_hashed_edge {
       ULONG operator()(const hashed_edge& e) const {
                if(e.source < e.target) 
                        return e.source*e.k + e.target;
                else 
                        return e.target*e.k + e.source;
       }
};

typedef extlib::unordered_map_with_custom_hash_and_comparator<hashed_edge, data_hashed_edge, hash_hashed_edge, compare_hashed_edge> hashed_graph;


#endif /* end of include guard: HASHED_GRAPH_DG1JG7O0 */
