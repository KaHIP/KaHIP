/******************************************************************************
 * two_way_refinement.h 
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

#ifndef TWO_WAY_REFINEMENT_INTERFACE_1ZWCSI0J
#define TWO_WAY_REFINEMENT_INTERFACE_1ZWCSI0J

#include "definitions.h"

class two_way_refinement{
        public:
                two_way_refinement( ) {};
                virtual ~two_way_refinement() {};

                virtual EdgeWeight perform_refinement(PartitionConfig & config, 
                                                      graph_access & G,
                                                      complete_boundary & boundary, 
                                                      std::vector<NodeID> & lhs_pq_start_nodes, 
                                                      std::vector<NodeID> & rhs_pq_start_nodes,
                                                      boundary_pair * refinement_pair,        
                                                      NodeWeight & lhs_part_weight,
                                                      NodeWeight & rhs_part_weight,
                                                      EdgeWeight & cut,
                                                      bool & something_changed) = 0;

};


#endif /* end of include guard: TWO_WAY_REFINEMENT_INTERFACE_1ZWCSI0J */
