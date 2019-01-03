/******************************************************************************
 * kway_graph_refinement.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef KWAY_GRAPH_REFINEMENT_PVGY97EW
#define KWAY_GRAPH_REFINEMENT_PVGY97EW

#include <vector>

#include "data_structure/priority_queues/priority_queue_interface.h"
#include "definitions.h"
#include "kway_graph_refinement_commons.h"
#include "random_functions.h"
#include "uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/vertex_moved_hashtable.h"
#include "uncoarsening/refinement/refinement.h"

class kway_graph_refinement : public refinement {
        public:
                kway_graph_refinement( );
                virtual ~kway_graph_refinement();

                EdgeWeight perform_refinement(PartitionConfig & config, 
                                              graph_access & G, 
                                              complete_boundary & boundary);

                void setup_start_nodes(PartitionConfig & config, 
                                       graph_access & G, 
                                       complete_boundary & boundary,  
                                       boundary_starting_nodes & start_nodes);
};

#endif /* end of include guard: KWAY_GRAPH_REFINEMENT_PVGY97EW */

