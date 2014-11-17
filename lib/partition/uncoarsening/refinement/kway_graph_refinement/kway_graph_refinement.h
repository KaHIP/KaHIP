/******************************************************************************
 * kway_graph_refinement.h 
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
                
        private:
                
                kway_graph_refinement_commons* commons;
};

#endif /* end of include guard: KWAY_GRAPH_REFINEMENT_PVGY97EW */

