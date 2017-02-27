/******************************************************************************
 * multitry_kway_fm.h
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

#ifndef MULTITRY_KWAYFM_PVGY97EW
#define MULTITRY_KWAYFM_PVGY97EW

#include <vector>

#include "definitions.h"
#include "kway_graph_refinement_commons.h"
#include "uncoarsening/refinement/refinement.h"

class multitry_kway_fm {
        public:
                multitry_kway_fm( );
                virtual ~multitry_kway_fm();

                int perform_refinement(PartitionConfig & config, graph_access & G, 
                                       complete_boundary & boundary, unsigned rounds, 
                                       bool init_neighbors, unsigned alpha);

                int perform_refinement_around_parts(PartitionConfig & config, graph_access & G, 
                                                    complete_boundary & boundary, bool init_neighbors, 
                                                    unsigned alpha, 
                                                    PartitionID & lhs, PartitionID & rhs,
                                                    std::unordered_map<PartitionID, PartitionID> & touched_blocks);


        private:
                int start_more_locallized_search(PartitionConfig & config, graph_access & G, 
                                                 complete_boundary & boundary, 
                                                 bool init_neighbors, 
                                                 bool compute_touched_blocks, 
                                                 std::unordered_map<PartitionID, PartitionID> & touched_blocks, 
                                                 std::vector<NodeID> & todolist);

                kway_graph_refinement_commons* commons;
};

#endif /* end of include guard: MULTITRY_KWAYFM_PVGY97EW  */


