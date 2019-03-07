/******************************************************************************
 * boundary_lookup.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef BOUNDARY_LOOKUP_2JMSKBSI
#define BOUNDARY_LOOKUP_2JMSKBSI

#include <unordered_map>

#include "definitions.h"
#include "limits.h"
#include "partial_boundary.h"

struct boundary_pair {
        PartitionID k;
        PartitionID lhs;
        PartitionID rhs;
};


struct compare_boundary_pair {
        bool operator()(const boundary_pair pair_a, const boundary_pair pair_b) const {
                bool eq = (pair_a.lhs == pair_b.lhs && pair_a.rhs == pair_b.rhs);
                     eq = eq || (pair_a.lhs == pair_b.rhs && pair_a.rhs == pair_b.lhs); 
                return eq;
        }
};

struct compare_boundary_pair_directed {
        bool operator()(const boundary_pair pair_a, const boundary_pair pair_b) const {
                bool eq = (pair_a.lhs == pair_b.lhs && pair_a.rhs == pair_b.rhs);
                return eq;
        }
};

struct data_boundary_pair {
        PartialBoundary pb_lhs;
        PartialBoundary pb_rhs;
        PartitionID lhs;
        PartitionID rhs;
        EdgeWeight edge_cut;

        bool initialized;

        data_boundary_pair() {
                edge_cut = 0;
                lhs = std::numeric_limits<PartitionID>::max();
                rhs = std::numeric_limits<PartitionID>::max();
                initialized = false;
        }
};

struct hash_boundary_pair_directed{
       size_t operator()(const boundary_pair pair) const {
                return pair.lhs*pair.k + pair.rhs;
       }
};

struct hash_boundary_pair{
       size_t operator()(const boundary_pair pair) const {
                if(pair.lhs < pair.rhs) 
                        return pair.lhs*pair.k + pair.rhs;
                else 
                        return pair.rhs*pair.k + pair.lhs;
       }
};

typedef std::unordered_map<const boundary_pair, data_boundary_pair, hash_boundary_pair, compare_boundary_pair> block_pairs;



#endif /* end of include guard: BOUNDARY_LOOKUP_2JMSKBSI */

