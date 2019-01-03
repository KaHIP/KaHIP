/******************************************************************************
 * partition_config.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PARTITION_CONFIG_DI1ES4T0
#define PARTITION_CONFIG_DI1ES4T0

#include "definitions.h"

// Configuration for the partitioning.
struct PartitionConfig
{
        PartitionConfig() {}


        //============================================================
        //=======================MATCHING=============================
        //============================================================
        bool edge_rating_tiebreaking;

        EdgeRating edge_rating;
        
        PermutationQuality permutation_quality;

        MatchingType matching_type;
        
        bool match_islands;

        bool first_level_random_matching;
        
        bool rate_first_level_inner_outer;

        NodeWeight max_vertex_weight; 
        
        NodeWeight largest_graph_weight; 

	NodeWeight work_load;

        unsigned aggressive_random_levels;
        
        bool disable_max_vertex_weight_constraint;

        //============================================================
        //===================INITIAL PARTITIONING=====================
        //============================================================
        unsigned int initial_partitioning_repetitions;

        unsigned int minipreps;

        bool refined_bubbling; 

        InitialPartitioningType initial_partitioning_type;

        bool initial_partition_optimize;

        BipartitionAlgorithm bipartition_algorithm;

        bool initial_partitioning;

        int bipartition_tries;

        int bipartition_post_fm_limits;

        int bipartition_post_ml_limits;

        //============================================================
        //====================REFINEMENT PARAMETERS===================
        //============================================================
        bool corner_refinement_enabled;

        bool use_bucket_queues;

        RefinementType refinement_type;

        PermutationQuality permutation_during_refinement;

        ImbalanceType imbalance;

        unsigned bubbling_iterations;
       
        unsigned kway_rounds; 
       
        bool quotient_graph_refinement_disabled; 

        KWayStopRule kway_stop_rule;

        double kway_adaptive_limits_alpha;

        double kway_adaptive_limits_beta;

        unsigned max_flow_iterations;

        unsigned local_multitry_rounds;
        
        unsigned local_multitry_fm_alpha;

        bool graph_allready_partitioned;

        unsigned int fm_search_limit;
        
        unsigned int kway_fm_search_limit;

        NodeWeight upper_bound_partition;

        double bank_account_factor;

        RefinementSchedulingAlgorithm refinement_scheduling_algorithm; 

        bool most_balanced_minimum_cuts;

        bool most_balanced_minimum_cuts_node_sep;
        
        unsigned toposort_iterations;

        bool softrebalance;

        bool rebalance;

        double flow_region_factor;

        bool gpa_grow_paths_between_blocks;

        //=======================================
        //==========GLOBAL SEARCH PARAMETERS=====
        //=======================================
        unsigned global_cycle_iterations;

        bool use_wcycles;

        bool use_fullmultigrid;

        unsigned level_split;

        bool no_new_initial_partitioning; 

        bool omit_given_partitioning; 

        StopRule stop_rule;

        int num_vert_stop_factor;
        
        bool no_change_convergence;

        //=======================================
        //===PERFECTLY BALANCED PARTITIONING ====
        //=======================================
	bool remove_negative_cycles;

        bool kaba_include_removal_of_paths;

        bool kaba_enable_zero_weight_cycles;

        double kabaE_internal_bal;

        CycleRefinementAlgorithm cycle_refinement_algorithm;

        int kaba_internal_no_aug_steps_aug;

        unsigned kaba_packing_iterations;

        bool kaba_flip_packings;

        MLSRule kaba_lsearch_p; // more localized search pseudo directed

        bool kaffpa_perfectly_balanced_refinement;

        unsigned kaba_unsucc_iterations;

        
        //=======================================
        //============PAR_PSEUDOMH / MH =========
        //=======================================
	double time_limit;

        double epsilon;

	unsigned no_unsuc_reps;

	unsigned local_partitioning_repetitions;

        bool mh_plain_repetitions;
        
        bool mh_easy_construction;

        bool mh_enable_gal_combine;

        bool mh_no_mh;

        bool mh_print_log;

        int  mh_flip_coin;

        int  mh_initial_population_fraction;

        bool mh_disable_cross_combine;

        bool mh_cross_combine_original_k;

        bool mh_disable_nc_combine;

        bool mh_disable_combine;

        bool mh_enable_quickstart;

        bool mh_disable_diversify_islands;

        bool mh_diversify;

        bool mh_diversify_best;

        bool mh_enable_tournament_selection;

        bool mh_optimize_communication_volume;

        unsigned mh_num_ncs_to_compute;

        unsigned mh_pool_size;

        bool combine; // in this case the second index is filled and edges between both partitions are not contracted

        unsigned initial_partition_optimize_fm_limits;

        unsigned initial_partition_optimize_multitry_fm_alpha;

        unsigned initial_partition_optimize_multitry_rounds;

        unsigned walshaw_mh_repetitions;

        unsigned scaleing_factor;

        bool scale_back;

	bool suppress_partitioner_output;

        unsigned maxT; 
        
        unsigned maxIter;
        //=======================================
        //===============BUFFOON=================
        //=======================================
        bool disable_hard_rebalance;

        bool buffoon;

        bool kabapE;
        
        bool mh_penalty_for_unconnected;
        //=======================================
        //===============MISC====================
        //=======================================
        std::string input_partition;

        int seed;

        bool fast;

        bool eco;

        bool strong;

        bool kaffpaE;

	bool balance_edges;

        // number of blocks the graph should be partitioned in
        PartitionID k;

        bool compute_vertex_separator;

        bool only_first_level;

        bool use_balance_singletons;

        int amg_iterations;

        std::string graph_filename;

        std::string filename_output;

        bool kaffpa_perfectly_balance;

        bool mode_node_separators;

        //=======================================
        //===========SNW PARTITIONING============
        //=======================================
        NodeOrderingType node_ordering;

        int cluster_coarsening_factor; 

        bool ensemble_clusterings; 

        int label_iterations;

        int label_iterations_refinement;

        int number_of_clusterings;

        bool label_propagation_refinement;

        double balance_factor;

        bool cluster_coarsening_during_ip;

        bool set_upperbound;

        int repetitions;

        //=======================================
        //===========NODE SEPARATOR==============
        //=======================================
        int max_flow_improv_steps;

        int max_initial_ns_tries;

        double region_factor_node_separators;

	bool sep_flows_disabled;

	bool sep_fm_disabled;

	bool sep_loc_fm_disabled;

        int sep_loc_fm_no_snodes;

	bool sep_greedy_disabled;

	int sep_fm_unsucc_steps;

	int sep_loc_fm_unsucc_steps;

	int sep_num_fm_reps;

	int sep_num_loc_fm_reps;

        int sep_num_vert_stop;

        bool sep_full_boundary_ip;

        bool faster_ns;

        EdgeRating sep_edge_rating_during_ip;

        //=======================================
        //=========LABEL PROPAGATION=============
        //=======================================
        NodeWeight cluster_upperbound;

        //=======================================
        //=========INITIAL PARTITIONING==========
        //=======================================

        // variables controling the size of the blocks during 
        // multilevel recursive bisection
        // (for the case where k is not a power of 2)
        std::vector<int> target_weights;

        bool initial_bipartitioning;

        int grow_target;

        //=======================================
        //===============QAP=====================
        //=======================================

        int communication_neighborhood_dist;

        LsNeighborhoodType ls_neighborhood;

        ConstructionAlgorithm construction_algorithm;

        DistanceConstructionAlgorithm distance_construction_algorithm;

        std::vector< int > group_sizes;

        std::vector< int > distances;

	int search_space_s;

        PreConfigMapping preconfiguration_mapping;

        int max_recursion_levels_construction; 

        bool enable_mapping;

        //=======================================
        //===============Shared Mem OMP==========
        //=======================================
        bool enable_omp;

        void LogDump(FILE *out) const {
        }
};


#endif /* end of include guard: PARTITION_CONFIG_DI1ES4T0 */
