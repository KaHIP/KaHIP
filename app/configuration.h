/******************************************************************************
 * configuration.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/


#ifndef CONFIGURATION_3APG5V7Z
#define CONFIGURATION_3APG5V7Z

#include "partition/partition_config.h"

class configuration {
        public:
                configuration() {} ;
                virtual ~configuration() {};

                void strong( PartitionConfig & config );
                void eco( PartitionConfig & config );
                void fast( PartitionConfig & config );

                void strong_separator( PartitionConfig & config );
                void eco_separator( PartitionConfig & config );
                void fast_separator( PartitionConfig & config );

                void standard( PartitionConfig & config );
                void standardsnw( PartitionConfig & config );

                void fastsocial( PartitionConfig & config );
                void ecosocial( PartitionConfig & config );
                void strongsocial( PartitionConfig & config ); 

                void fastsocial_separator( PartitionConfig & config );
                void ecosocial_separator( PartitionConfig & config );
                void strongsocial_separator( PartitionConfig & config ); 

};

inline void configuration::strong( PartitionConfig & partition_config ) {
        standard(partition_config);
        partition_config.matching_type                          = MATCHING_GPA;
        partition_config.permutation_quality                    = PERMUTATION_QUALITY_GOOD;
        partition_config.permutation_during_refinement          = PERMUTATION_QUALITY_GOOD;
        partition_config.edge_rating_tiebreaking                = true;
        partition_config.fm_search_limit                        = 5;
        partition_config.bank_account_factor                    = 3;
        partition_config.edge_rating                            = EXPANSIONSTAR2;
        partition_config.refinement_scheduling_algorithm        = REFINEMENT_SCHEDULING_ACTIVE_BLOCKS_REF_KWAY;
        partition_config.refinement_type                        = REFINEMENT_TYPE_FM_FLOW;
        partition_config.global_cycle_iterations                = 2;
        partition_config.flow_region_factor                     = 8;
        partition_config.corner_refinement_enabled              = true;
        partition_config.kway_stop_rule                         = KWAY_ADAPTIVE_STOP_RULE;
        partition_config.kway_adaptive_limits_alpha             = 10;
        partition_config.kway_rounds                            = 10;
        partition_config.rate_first_level_inner_outer           = true;
        partition_config.use_wcycles                            = false; 
        partition_config.no_new_initial_partitioning            = true;
        partition_config.use_fullmultigrid                      = true;
        partition_config.most_balanced_minimum_cuts             = true;
        partition_config.local_multitry_fm_alpha                = 10;
        partition_config.local_multitry_rounds                  = 10;

        partition_config.mh_initial_population_fraction         = 10;
        partition_config.mh_flip_coin                           = 1; 
#ifndef MODE_NODESEP
        partition_config.epsilon                                = 3; 
        partition_config.imbalance                              = 3;
#else
        partition_config.epsilon                                = 20; 
        partition_config.imbalance                              = 20;
#endif

        partition_config.initial_partitioning_type              = INITIAL_PARTITIONING_RECPARTITION;
        partition_config.bipartition_tries                      = 4;
        partition_config.minipreps                              = 4;
        partition_config.initial_partitioning_repetitions       = 64;
        partition_config.strong = true;

}

inline void configuration::eco( PartitionConfig & partition_config ) {
        standard(partition_config);
        partition_config.eco                      = true;
        partition_config.aggressive_random_levels = std::max(2, (int)(7 - log2(partition_config.k)));
        
        partition_config.kway_rounds                            = std::min(5, (int)log2(partition_config.k));
        partition_config.matching_type                          = MATCHING_RANDOM_GPA;
        partition_config.permutation_quality                    = PERMUTATION_QUALITY_NONE;
        partition_config.permutation_during_refinement          = PERMUTATION_QUALITY_GOOD;
        partition_config.edge_rating                            = EXPANSIONSTAR2;
        partition_config.fm_search_limit                        = 1;
        partition_config.refinement_type                        = REFINEMENT_TYPE_FM_FLOW;
        partition_config.flow_region_factor                     = 2;
        partition_config.corner_refinement_enabled              = true;
        partition_config.kway_stop_rule                         = KWAY_SIMPLE_STOP_RULE;
        partition_config.kway_fm_search_limit                   = 1;
        partition_config.mh_initial_population_fraction         = 50;
        partition_config.mh_flip_coin                           = 1; 

        partition_config.initial_partitioning_type              = INITIAL_PARTITIONING_RECPARTITION;
        partition_config.bipartition_tries                      = 4;
        partition_config.minipreps                              = 4;
        partition_config.initial_partitioning_repetitions       = 16;
}

inline void configuration::fast( PartitionConfig & partition_config ) {
        standard(partition_config);

        partition_config.fast = true;
        if(partition_config.k > 8) {
                partition_config.quotient_graph_refinement_disabled     = true;
                partition_config.kway_fm_search_limit                   = 0; 
                partition_config.kway_stop_rule                         = KWAY_SIMPLE_STOP_RULE; 
                partition_config.corner_refinement_enabled              = true; 
        } else {
                partition_config.corner_refinement_enabled              = false; 
        }
        partition_config.permutation_quality                    = PERMUTATION_QUALITY_FAST;
        partition_config.permutation_during_refinement          = PERMUTATION_QUALITY_NONE;
        partition_config.matching_type                          = MATCHING_RANDOM_GPA;
        partition_config.aggressive_random_levels               = 4;
        partition_config.refinement_scheduling_algorithm        = REFINEMENT_SCHEDULING_FAST;
        partition_config.edge_rating                            = EXPANSIONSTAR2;
        partition_config.fm_search_limit                        = 0;
        partition_config.bank_account_factor                    = 1;

        partition_config.initial_partitioning_type              = INITIAL_PARTITIONING_RECPARTITION;
        partition_config.bipartition_tries                      = 4;
        partition_config.minipreps                              = 1;
        partition_config.initial_partitioning_repetitions       = 0;

}

inline void configuration::strong_separator( PartitionConfig & partition_config ) {
        standard(partition_config);
        partition_config.matching_type                          = MATCHING_GPA;
        partition_config.permutation_quality                    = PERMUTATION_QUALITY_GOOD;
        partition_config.permutation_during_refinement          = PERMUTATION_QUALITY_GOOD;
        partition_config.edge_rating_tiebreaking                = true;
        partition_config.fm_search_limit                        = 5;
        partition_config.bank_account_factor                    = 3;
        partition_config.edge_rating                            = EXPANSIONSTAR2;
        partition_config.refinement_scheduling_algorithm        = REFINEMENT_SCHEDULING_ACTIVE_BLOCKS_REF_KWAY;
        partition_config.refinement_type                        = REFINEMENT_TYPE_FM_FLOW;
        partition_config.global_cycle_iterations                = 3;
        partition_config.flow_region_factor                     = 8;
        partition_config.corner_refinement_enabled              = true;
        partition_config.kway_stop_rule                         = KWAY_ADAPTIVE_STOP_RULE;
        partition_config.kway_adaptive_limits_alpha             = 10;
        partition_config.kway_rounds                            = 10;
        partition_config.rate_first_level_inner_outer           = true;
        partition_config.use_wcycles                            = false; 
        partition_config.no_new_initial_partitioning            = true;
        partition_config.use_fullmultigrid                      = true;
        partition_config.most_balanced_minimum_cuts             = true;
        partition_config.most_balanced_minimum_cuts_node_sep    = false;
        partition_config.local_multitry_fm_alpha                = 10;
        partition_config.local_multitry_rounds                  = 10;

        partition_config.mh_initial_population_fraction         = 10;
        partition_config.mh_flip_coin                           = 1; 

#ifndef MODE_NODESEP
        partition_config.epsilon                                = 3; 
        partition_config.imbalance                              = 3;
#else
        partition_config.epsilon                                = 20; 
        partition_config.imbalance                              = 20;
#endif

        partition_config.initial_partitioning_type              = INITIAL_PARTITIONING_RECPARTITION;
        partition_config.bipartition_tries                      = 4;
        partition_config.minipreps                              = 4;
        partition_config.initial_partitioning_repetitions       = 64;
        partition_config.strong = true;

        partition_config.mode_node_separators = true;
        partition_config.use_fullmultigrid    = false;
        partition_config.use_wcycles          = false;
        partition_config.matching_type        = MATCHING_GPA;
        partition_config.global_cycle_iterations = 3;
        partition_config.region_factor_node_separators = 1;

}

inline void configuration::eco_separator( PartitionConfig & partition_config ) {
        standard(partition_config);
        partition_config.eco                      = true;
        partition_config.aggressive_random_levels = std::max(2, (int)(7 - log2(partition_config.k)));
        
        partition_config.kway_rounds                            = std::min(5, (int)log2(partition_config.k));
        partition_config.matching_type                          = MATCHING_RANDOM_GPA;
        partition_config.permutation_quality                    = PERMUTATION_QUALITY_NONE;
        partition_config.permutation_during_refinement          = PERMUTATION_QUALITY_GOOD;
        partition_config.edge_rating                            = EXPANSIONSTAR2;
        partition_config.fm_search_limit                        = 1;
        partition_config.refinement_type                        = REFINEMENT_TYPE_FM_FLOW;
        partition_config.flow_region_factor                     = 2;
        partition_config.corner_refinement_enabled              = true;
        partition_config.kway_stop_rule                         = KWAY_SIMPLE_STOP_RULE;
        partition_config.kway_fm_search_limit                   = 1;
        partition_config.mh_initial_population_fraction         = 50;
        partition_config.mh_flip_coin                           = 1; 

        partition_config.initial_partitioning_type              = INITIAL_PARTITIONING_RECPARTITION;
        partition_config.bipartition_tries                      = 4;
        partition_config.minipreps                              = 4;
        partition_config.initial_partitioning_repetitions       = 16;

        partition_config.mode_node_separators = true;
        partition_config.use_fullmultigrid    = false;
        partition_config.use_wcycles          = false;
        partition_config.matching_type        = MATCHING_GPA;
	partition_config.sep_loc_fm_disabled  = true;
	partition_config.sep_flows_disabled   = false;
	partition_config.sep_fm_disabled      = false;
	partition_config.sep_greedy_disabled  = true;
	partition_config.max_simplicial_degree = 12;
        partition_config.region_factor_node_separators = 0.5;
        partition_config.global_cycle_iterations = 2;

}

inline void configuration::fast_separator( PartitionConfig & partition_config ) {
        standard(partition_config);

        partition_config.fast = true;
        if(partition_config.k > 8) {
                partition_config.quotient_graph_refinement_disabled     = true;
                partition_config.kway_fm_search_limit                   = 0; 
                partition_config.kway_stop_rule                         = KWAY_SIMPLE_STOP_RULE; 
                partition_config.corner_refinement_enabled              = true; 
        } else {
                partition_config.corner_refinement_enabled              = false; 
        }
        partition_config.permutation_quality                    = PERMUTATION_QUALITY_FAST;
        partition_config.permutation_during_refinement          = PERMUTATION_QUALITY_NONE;
        partition_config.matching_type                          = MATCHING_RANDOM_GPA;
        partition_config.aggressive_random_levels               = 4;
        partition_config.refinement_scheduling_algorithm        = REFINEMENT_SCHEDULING_FAST;
        partition_config.edge_rating                            = EXPANSIONSTAR2;
        partition_config.fm_search_limit                        = 0;
        partition_config.bank_account_factor                    = 1;

        partition_config.initial_partitioning_type              = INITIAL_PARTITIONING_RECPARTITION;
        partition_config.bipartition_tries                      = 4;
        partition_config.minipreps                              = 1;
        partition_config.initial_partitioning_repetitions       = 0;

        partition_config.mode_node_separators    = true;
        partition_config.use_fullmultigrid       = false;
        partition_config.use_wcycles             = false;
        partition_config.matching_type           = MATCHING_GPA;
        partition_config.sep_loc_fm_disabled     = true;
        partition_config.sep_flows_disabled      = true;
        partition_config.sep_fm_disabled         = false;
        partition_config.sep_greedy_disabled     = true;
        partition_config.global_cycle_iterations = 1;
}



inline void configuration::standard( PartitionConfig & partition_config ) {
        partition_config.filename_output                        = "";
        partition_config.use_mmap_io = false;
        partition_config.seed                                   = 0;
        partition_config.fast                                   = false;
        partition_config.mode_node_separators                   = false;
        partition_config.eco                                    = false;
        partition_config.strong                                 = false;
        partition_config.first_level_random_matching            = false;
        partition_config.initial_partitioning_repetitions       = 5;
        partition_config.edge_rating_tiebreaking                = false;
        partition_config.edge_rating                            = WEIGHT;
        partition_config.matching_type                          = MATCHING_RANDOM;
        partition_config.permutation_quality                    = PERMUTATION_QUALITY_FAST;
        partition_config.initial_partitioning                   = false;
        partition_config.initial_partitioning_type              = INITIAL_PARTITIONING_RECPARTITION;
        partition_config.bipartition_tries                      = 9;
        partition_config.minipreps                              = 10;
        partition_config.enable_omp                             = false;
        partition_config.combine                                = false;
#ifndef MODE_NODESEP
        partition_config.epsilon                                = 3; 
        partition_config.imbalance                              = 3;
#else
        partition_config.epsilon                                = 20; 
        partition_config.imbalance                              = 20;
#endif

        partition_config.buffoon                                = false;
        partition_config.balance_edges                          = false;

        partition_config.time_limit 				= 0; 
        partition_config.mh_pool_size                           = 5;
        partition_config.mh_plain_repetitions                   = false;
        partition_config.no_unsuc_reps				= 10;
        partition_config.local_partitioning_repetitions 	= 1;

        partition_config.mh_disable_nc_combine                  = false;
        partition_config.mh_disable_cross_combine               = false;
        partition_config.mh_disable_combine                     = false;
        partition_config.mh_enable_quickstart                   = false; 
        partition_config.mh_disable_diversify_islands           = false;
        partition_config.mh_diversify                           = true;
        partition_config.mh_diversify_best                      = false;
        partition_config.mh_cross_combine_original_k            = false;
        partition_config.mh_enable_tournament_selection         = true;
        partition_config.mh_initial_population_fraction         = 10;
        partition_config.mh_flip_coin                           = 1;
        partition_config.mh_print_log                           = false;
        partition_config.mh_penalty_for_unconnected             = false;
        partition_config.mh_no_mh                               = false;
        partition_config.mh_optimize_communication_volume       = false; 
        partition_config.use_bucket_queues                      = true; 
        partition_config.walshaw_mh_repetitions                 = 50;
        partition_config.scaleing_factor                        = 1;
        partition_config.scale_back                             = false;
        partition_config.initial_partition_optimize_fm_limits   = 20;
        partition_config.initial_partition_optimize_multitry_fm_alpha = 20;
        partition_config.initial_partition_optimize_multitry_rounds   = 100;

        partition_config.suppress_partitioner_output		= false;

        if( partition_config.k <= 4 ) { 
                partition_config.bipartition_post_fm_limits             = 30;
                partition_config.bipartition_post_ml_limits             = 6;
        } else {
                partition_config.bipartition_post_fm_limits             = 25;
                partition_config.bipartition_post_ml_limits             = 5;
        }

        partition_config.disable_max_vertex_weight_constraint   = false;
        partition_config.permutation_during_refinement          = PERMUTATION_QUALITY_GOOD;
        partition_config.fm_search_limit                        = 5;
        partition_config.use_bucket_queues                      = false;
        partition_config.bank_account_factor                    = 1.5;
        partition_config.refinement_scheduling_algorithm        = REFINEMENT_SCHEDULING_ACTIVE_BLOCKS;
        partition_config.rate_first_level_inner_outer           = false;
        partition_config.match_islands                          = false;
        partition_config.refinement_type                        = REFINEMENT_TYPE_FM;
        partition_config.flow_region_factor                     = 4.0;
        partition_config.aggressive_random_levels               = 3;
        partition_config.refined_bubbling                       = true;
        partition_config.corner_refinement_enabled              = false;
        partition_config.bubbling_iterations                    = 1;
        partition_config.kway_rounds                            = 1;
        partition_config.quotient_graph_refinement_disabled     = false;
        partition_config.kway_fm_search_limit                   = 3;
        partition_config.global_cycle_iterations                = 1;
        partition_config.softrebalance                          = false;
        partition_config.rebalance                              = false;
        partition_config.use_wcycles                            = false;
        partition_config.stop_rule                              = STOP_RULE_SIMPLE;
        partition_config.num_vert_stop_factor                   = 20;
        partition_config.level_split                            = 2;
        partition_config.no_new_initial_partitioning            = true;
        partition_config.omit_given_partitioning                = false;
        partition_config.use_fullmultigrid                      = false;
        partition_config.kway_stop_rule                         = KWAY_SIMPLE_STOP_RULE;
        partition_config.kway_adaptive_limits_alpha             = 1.0;
        partition_config.max_flow_iterations                    = 10;
        partition_config.no_change_convergence                  = false;
        partition_config.compute_vertex_separator               = false;
        partition_config.toposort_iterations                    = 4;
        partition_config.initial_partition_optimize             = false;
        partition_config.most_balanced_minimum_cuts             = false;
        partition_config.most_balanced_minimum_cuts_node_sep    = false;
        partition_config.gpa_grow_paths_between_blocks          = true; //used for the kaffpa paper

        partition_config.bipartition_algorithm                  = BIPARTITION_BFS;
        partition_config.local_multitry_rounds                  = 1;
        partition_config.local_multitry_fm_alpha                = 10;

        partition_config.only_first_level                       = false;
        partition_config.use_balance_singletons                 = true;

        partition_config.disable_hard_rebalance                 = false;
        partition_config.amg_iterations                         = 5;
        partition_config.kaffpa_perfectly_balance               = false;
        partition_config.kaffpaE                                = false;

        partition_config.remove_negative_cycles			= false;
        partition_config.cycle_refinement_algorithm             = CYCLE_REFINEMENT_ALGORITHM_ULTRA_MODEL;
        partition_config.kabaE_internal_bal                     = 0.01;
        partition_config.kaba_packing_iterations                = 20;
        partition_config.kaba_flip_packings                     = false;
        partition_config.kaba_lsearch_p                         = NOCOIN_RNDTIE;
        partition_config.kaffpa_perfectly_balanced_refinement   = false;
        partition_config.kaba_enable_zero_weight_cycles         = true;
        partition_config.mh_enable_gal_combine			= false;
        partition_config.mh_easy_construction                   = false;
        partition_config.faster_ns                              = false;

        partition_config.maxT                                   = 100;
        partition_config.maxIter                                = 500000;

        if( partition_config.k <= 8 ) {
               partition_config.kaba_internal_no_aug_steps_aug = 15; 
        } else {
               partition_config.kaba_internal_no_aug_steps_aug = 7; 
        }

        partition_config.kaba_unsucc_iterations = 6;
	partition_config.initial_bipartitioning = false;
        partition_config.kabapE = false;


#ifdef MODE_BUFFOON

        partition_config.mh_diversify         = true;
        partition_config.kway_fm_search_limit = 10;
        partition_config.fm_search_limit      = 10;
        partition_config.initial_partition_optimize = true;
        partition_config.nc_div               = 2;
        partition_config.nc_b                 = 1;
#endif

        // social networking parameters
        partition_config.cluster_coarsening_factor    = 18;
        partition_config.ensemble_clusterings         = false;
        partition_config.label_iterations             = 10;
        partition_config.label_iterations_refinement  = 25;
        partition_config.number_of_clusterings        = 1;
        partition_config.label_propagation_refinement = false;
        partition_config.balance_factor               = 0;
        partition_config.cluster_coarsening_during_ip = false;
        partition_config.set_upperbound               = true;
        partition_config.repetitions                  = 1;
        partition_config.node_ordering                = DEGREE_NODEORDERING;

        //node separator parameters
        partition_config.max_flow_improv_steps         = 5;
        partition_config.max_initial_ns_tries          = 25;
        partition_config.region_factor_node_separators = 0.5;
	partition_config.sep_flows_disabled  	       = false;
	partition_config.sep_fm_disabled  	       = false;
	partition_config.sep_greedy_disabled  	       = true;

	partition_config.sep_fm_unsucc_steps           = 2000;
        partition_config.sep_num_fm_reps               = 200;

	partition_config.sep_loc_fm_disabled  	       = false;
	partition_config.sep_loc_fm_no_snodes          = 20;
	partition_config.sep_loc_fm_unsucc_steps       = 50;
        partition_config.sep_num_loc_fm_reps           = 25;
        partition_config.sep_num_vert_stop             = 8000;
        partition_config.sep_full_boundary_ip          = false;
        partition_config.sep_edge_rating_during_ip     = SEPARATOR_MULTX;

        partition_config.enable_mapping                    = false;
        partition_config.ls_neighborhood                   = COMMUNICATIONGRAPH;
        partition_config.communication_neighborhood_dist   = 10;
        partition_config.construction_algorithm            = MAP_CONST_FASTHIERARCHY_TOPDOWN;
        partition_config.distance_construction_algorithm   = DIST_CONST_HIERARCHY;
        partition_config.search_space_s                    = 64;
        partition_config.preconfiguration_mapping          = PRE_CONFIG_MAPPING_ECO;
        partition_config.max_recursion_levels_construction = std::numeric_limits< int >::max();


        partition_config.group_sizes.push_back(4);
        partition_config.group_sizes.push_back(8);
        partition_config.group_sizes.push_back(8);

        partition_config.distances.push_back(1);
        partition_config.distances.push_back(10);
        partition_config.distances.push_back(100);

        partition_config.convergence_factor = 1;
        partition_config.reduction_order = {simplicial_nodes,
                                            degree_2_nodes};

        //partition_config.reduction_order = {simplicial_nodes,
                //indistinguishable_nodes,
                ////twins,
                ////path_compression,
                //degree_2_nodes,
                //triangle_contraction};

        partition_config.dissection_rec_limit = 120;
        partition_config.disable_reductions = false;

        partition_config.ilp_min_gain = -1;
        partition_config.ilp_bfs_depth = 2;
        partition_config.ilp_overlap_presets = OverlapPresets::NOEQUAL;

        partition_config.ilp_limit_nonzeroes = 5000000;
        partition_config.ilp_overlap_runs = 3;
        partition_config.ilp_timeout = 7200;
}

inline void configuration::standardsnw( PartitionConfig & partition_config ) {
        partition_config.matching_type        = CLUSTER_COARSENING;
        partition_config.stop_rule            = STOP_RULE_MULTIPLE_K;
        partition_config.num_vert_stop_factor = 5000;

        if(2 <= partition_config.k && partition_config.k <= 3) {
                partition_config.number_of_clusterings = 18;
        } else if(4 <= partition_config.k && partition_config.k <= 7) {
                partition_config.number_of_clusterings = 17;
        } else if(8 <= partition_config.k && partition_config.k <= 15) {
                partition_config.number_of_clusterings = 15;
        } else if(16 <= partition_config.k && partition_config.k <= 31) {
                partition_config.number_of_clusterings = 7;
        } else { 
                partition_config.number_of_clusterings = 3;
        }

        partition_config.balance_factor               = 0.016;
        if( partition_config.k <= 8 ) {
                partition_config.balance_factor       = 0.00;
        }


}

inline void configuration::fastsocial_separator( PartitionConfig & partition_config ) {
        eco(partition_config);
        standardsnw(partition_config);
        
        partition_config.label_propagation_refinement = true;
        partition_config.cluster_coarsening_during_ip = true;
        partition_config.balance_factor               = 0;

        partition_config.mode_node_separators    = true;
        partition_config.use_fullmultigrid       = false;
        partition_config.use_wcycles             = false;
        partition_config.matching_type           = MATCHING_GPA;
        partition_config.sep_loc_fm_disabled     = true;
        partition_config.sep_flows_disabled      = true;
        partition_config.sep_fm_disabled         = false;
        partition_config.sep_greedy_disabled     = true;
        partition_config.global_cycle_iterations = 1;
}

inline void configuration::ecosocial_separator( PartitionConfig & partition_config ) {
        eco(partition_config);
        standardsnw(partition_config);
        partition_config.label_propagation_refinement = false;
        partition_config.global_cycle_iterations      = 3;
        partition_config.use_wcycles                  = false; 
        partition_config.no_new_initial_partitioning  = true;
        partition_config.balance_factor               = 0.016;
        partition_config.cluster_coarsening_during_ip = true;
}

inline void configuration::strongsocial_separator( PartitionConfig & partition_config ) {
        strong(partition_config);
        standardsnw(partition_config);

        partition_config.label_propagation_refinement = false;
        partition_config.cluster_coarsening_during_ip = true;
        partition_config.ensemble_clusterings         = true;

        partition_config.mode_node_separators = true;
        partition_config.use_fullmultigrid    = false;
        partition_config.use_wcycles          = false;
        partition_config.matching_type        = MATCHING_GPA;
        partition_config.global_cycle_iterations = 3;
}


inline void configuration::fastsocial( PartitionConfig & partition_config ) {
        eco(partition_config);
        standardsnw(partition_config);
        partition_config.label_propagation_refinement = true;
        partition_config.cluster_coarsening_during_ip = true;
        partition_config.balance_factor               = 0;
}

inline void configuration::ecosocial( PartitionConfig & partition_config ) {
        eco(partition_config);
        standardsnw(partition_config);
        partition_config.label_propagation_refinement = false;
        partition_config.global_cycle_iterations      = 3;
        partition_config.use_wcycles                  = false; 
        partition_config.no_new_initial_partitioning  = true;
        partition_config.balance_factor               = 0.016;
        partition_config.cluster_coarsening_during_ip = true;
}

inline void configuration::strongsocial( PartitionConfig & partition_config ) {
        strong(partition_config);
        standardsnw(partition_config);

        partition_config.label_propagation_refinement = false;
        partition_config.cluster_coarsening_during_ip = true;
        partition_config.ensemble_clusterings         = true;

}

#endif /* end of include guard: CONFIGURATION_3APG5V7Z */
