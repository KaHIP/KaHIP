/******************************************************************************
 * quotient_graph_refinement.h
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

#ifndef QUOTIENT_GRAPH_REFINEMENT_A0Y1Y6LL
#define QUOTIENT_GRAPH_REFINEMENT_A0Y1Y6LL

#include "definitions.h"
#include "uncoarsening/refinement/refinement.h"

class quotient_graph_refinement : public refinement {
        public:
                quotient_graph_refinement( );
                virtual ~quotient_graph_refinement();

                EdgeWeight perform_refinement(PartitionConfig & config, graph_access & G, complete_boundary & boundary);

                void setup_start_nodes(graph_access & G, 
                                       PartitionID partition, 
                                       boundary_pair & bp, 
                                       complete_boundary & boundary,  
                                       boundary_starting_nodes & start_nodes);

        private:
                EdgeWeight perform_a_two_way_refinement(PartitionConfig & config, 
                                                        graph_access & G,
                                                        complete_boundary & boundary, 
                                                        boundary_pair & bp,
                                                        PartitionID & lhs, 
                                                        PartitionID & rhs,
                                                        NodeWeight & lhs_part_weight,
                                                        NodeWeight & rhs_part_weight,
                                                        EdgeWeight & cut,
                                                        bool & something_changed); 

};


#endif /* end of include guard: QUOTIENT_GRAPH_REFINEMENT_A0Y1Y6LL */
