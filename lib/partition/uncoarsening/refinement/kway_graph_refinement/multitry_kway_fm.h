/******************************************************************************
 * multitry_kway_fm.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
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


