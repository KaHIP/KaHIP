/******************************************************************************
 * cycle_definitions.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef CYCLE_DEFINITIONS_4GQMW8PZ
#define CYCLE_DEFINITIONS_4GQMW8PZ

#include <unordered_map>

#include "uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"

struct undo_struct {
        NodeID node;
        PartitionID to;
};

struct data_qgraph_edge {
        NodeID to_move;
        Gain gain;

        data_qgraph_edge() {
                to_move = std::numeric_limits<PartitionID>::max();
                gain    = std::numeric_limits<PartitionID>::min();
        }
};

typedef std::unordered_map<const boundary_pair, data_qgraph_edge, hash_boundary_pair_directed, compare_boundary_pair_directed> edge_movements;


#endif /* end of include guard: DEFINITIONS_4GQMW8PZ */
