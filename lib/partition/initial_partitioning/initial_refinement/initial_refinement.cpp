/******************************************************************************
 * initial_refinement.cpp 
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

#include "initial_refinement.h"
#include "coarsening/coarsening.h"
#include "uncoarsening/uncoarsening.h"

initial_refinement::initial_refinement() {
                
}

initial_refinement::~initial_refinement() {
                
}

int initial_refinement::optimize( const PartitionConfig & config, graph_access & G, EdgeWeight & initial_cut) {

        PartitionConfig partition_config                      = config;
        partition_config.graph_allready_partitioned           = true;
        partition_config.stop_rule                            = STOP_RULE_STRONG;
        partition_config.fm_search_limit                      = partition_config.initial_partition_optimize_fm_limits;
        partition_config.kway_fm_search_limit                 = partition_config.initial_partition_optimize_fm_limits;
        partition_config.local_multitry_fm_alpha              = partition_config.initial_partition_optimize_multitry_fm_alpha;
        partition_config.local_multitry_rounds                = partition_config.initial_partition_optimize_multitry_rounds;
        partition_config.matching_type                        = MATCHING_GPA;
        partition_config.gpa_grow_paths_between_blocks        = false;
        partition_config.kaffpa_perfectly_balanced_refinement = false; // for runtime reasons
        
        graph_hierarchy hierarchy;

        coarsening coarsen;
        coarsen.perform_coarsening(partition_config, G, hierarchy);

        //ommit initial partitioning since we have the partition allread given
        uncoarsening uncoarsen;
        int improvement = 0;
        improvement     = uncoarsen.perform_uncoarsening(partition_config, hierarchy);
        initial_cut    -= improvement;

        return improvement;
} 
