/******************************************************************************
 * cycle_definitions.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
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
