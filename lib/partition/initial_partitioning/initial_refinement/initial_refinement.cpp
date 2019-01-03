/******************************************************************************
 * initial_refinement.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
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
