/******************************************************************************
 * most_balanced_minimum_cuts.h 
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

#ifndef MOST_BALANCED_MINIMUM_CUTS_SBD5CS
#define MOST_BALANCED_MINIMUM_CUTS_SBD5CS

#include "data_structure/graph_access.h"
#include "partition_config.h"

class most_balanced_minimum_cuts {
        public:
                most_balanced_minimum_cuts();
                virtual ~most_balanced_minimum_cuts();

                void compute_good_balanced_min_cut( graph_access & residualGraph,
                                                    const PartitionConfig & config,
                                                    NodeWeight & perfect_rhs_weight, 
                                                    std::vector< NodeID > & new_rhs_node ); 

        private:
                void build_internal_scc_graph( graph_access & residualGraph,  
                                               std::vector<int> & components, 
                                               int comp_count, 
                                               graph_access & scc_graph);

                void compute_new_rhs( graph_access & scc_graph, 
                                      const PartitionConfig & config,
                                      std::vector< NodeWeight > & comp_weights,
                                      int comp_of_s,
                                      int comp_of_t,
                                      NodeWeight optimal_rhs_weight,
                                      std::vector<int> & comp_for_rhs); 
};


#endif /* end of include guard: MOST_BALANCED_MINIMUM_CUTS_SBD5CS */
