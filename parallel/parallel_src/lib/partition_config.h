/******************************************************************************
 * partition_config.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PARTITION_CONFIG_DI1ES4T0A
#define PARTITION_CONFIG_DI1ES4T0A

#include "definitions.h"

// Configuration for the partitioning.
struct PPartitionConfig
{
        PPartitionConfig() {}

        //=======================================
        //============ Graph Gen=================
        //=======================================
        int log_num_verts;

        long edge_factor;

        bool generate_rgg;

        bool generate_ba;

        //=======================================
        //============ Communication ============
        //=======================================

        ULONG comm_rounds;

        //=======================================
        //============ Global Data===============
        //=======================================

        NodeID number_of_overall_nodes;

        //=======================================
        //===============MISC====================
        //=======================================

        PermutationQuality permutation_quality;

        unsigned int label_iterations;
        
        unsigned int label_iterations_coarsening;

        unsigned int label_iterations_refinement;

        double cluster_coarsening_factor;

        double time_limit;

        unsigned epsilon;

        unsigned inbalance;

        std::string input_partition;

        int seed;

        PartitionID k;

        std::string graph_filename;

        std::string input_partition_filename;

        int evolutionary_time_limit;

	NodeWeight upper_bound_partition;

	NodeWeight upper_bound_cluster;

        NodeID total_num_labels;

        InitialPartitioningAlgorithm initial_partitioning_algorithm;

        int stop_factor;

        bool vcycle;

        int num_vcycles;

        int num_tries; // number of repetitions to perform

        NodeOrderingType node_ordering;

        bool no_refinement_in_last_iteration;

        double ht_fill_factor;

        bool eco;

	int binary_io_window_size;

        ULONG barabasi_albert_mindegree;

        bool compute_degree_sequence_ba;

        bool compute_degree_sequence_k_first;

        bool kronecker_internal_only;

        ULONG k_deg;

        bool generate_ba_32bit;

        ULONG n;

	bool save_partition;

	bool save_partition_binary;

        bool vertex_degree_weights;

        bool converter_evaluate;

        //=======================================
        //===============Shared Mem OMP==========
        //=======================================
        void LogDump(FILE *out) const {
        }
};


#endif /* end of include guard: PARTITION_CONFIG_DI1ES4T0 */
