/******************************************************************************
 * parse_parameters.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/


#ifndef PARSE_PARAMETERS_GPJMGSM8A
#define PARSE_PARAMETERS_GPJMGSM8A

#include <regex.h>
#include <string.h>

#include "configuration.h"

// #include "data_structure/matrix/online_distance_matrix.h"
// #include "data_structure/matrix/online_precalc_matrix.h"
// #include "data_structure/matrix/online_binary_matrix.h"
// #include "data_structure/matrix/full_matrix.h"

int parse_parameters(int argn, char **argv, 
                     PPartitionConfig & partition_config, 
                     std::string & graph_filename) {

        const char *progname = argv[0];

        // Setup argtable parameters.
        struct arg_lit *help                           = arg_lit0(NULL, "help","Print help.");
        struct arg_str *filename                       = arg_str1(NULL, NULL, "FILE", "Path to graph file to partition.");
        struct arg_str *filename_output                = arg_str0(NULL, "filename_output", NULL, "Specify the name of the output file (that contains the partition).");
        struct arg_str *input_partition_filename       = arg_str1(NULL, "input_partition", "FILE", "Path to partition file to convert.");
        struct arg_int *user_seed                      = arg_int0(NULL, "seed", NULL, "Seed to use for the PRNG.");
        struct arg_int *k                              = arg_int1(NULL, "k", NULL, "Number of blocks to partition the graph.");
        struct arg_int *k_opt                          = arg_int0(NULL, "k", NULL, "Number of blocks to partition the graph.");
        struct arg_int *inbalance                      = arg_int0(NULL, "imbalance", NULL, "Desired balance. Default: 3 (%).");
        struct arg_int *comm_rounds                    = arg_int0(NULL, "comm_rounds", NULL, "Number of communication rounds per complete graph iteration.");
        struct arg_dbl *cluster_coarsening_factor      = arg_dbl0(NULL, "cluster_coarsening_factor", NULL, "The coarsening factor basically involes a bound on the block weights.");
        struct arg_int *stop_factor                    = arg_int0(NULL, "stop_factor", NULL, "Stop factor l to stop coarsening if total num vert <= lk.");
        struct arg_int *evolutionary_time_limit        = arg_int0(NULL, "evolutionary_time_limit", NULL, "Time limit for the evolutionary algorithm.");
#ifndef TOOLBOX
#endif
        struct arg_int *label_iterations_coarsening    = arg_int0(NULL, "label_iterations_coarsening", NULL, "Number of label propagation iterations during coarsening.");
        struct arg_int *label_iterations_refinement    = arg_int0(NULL, "label_iterations_refinement", NULL, "Number of label propagation iterations during refinement.");
        struct arg_int *num_tries                      = arg_int0(NULL, "num_tries", NULL, "Number of repetitions to perform.");
        struct arg_int *binary_io_window_size 	       = arg_int0(NULL, "binary_io_window_size", NULL, "Binary IO window size.");
        struct arg_rex *initial_partitioning_algorithm = arg_rex0(NULL, "initial_partitioning_algorithm", "^(kaffpaEstrong|kaffpaEeco|kaffpaEfast|fastsocial|ecosocial|strongsocial|random)$", "PARTITIONER", REG_EXTENDED, "Initial partitioning algorithm to use. One of {kaffpaEstrong, kaffpaEeco, kaffpaEfast, fastsocial, ecosocial, strongsocial, random}." );
        struct arg_int *num_vcycles                    = arg_int0(NULL, "num_vcycles", NULL, "Number of vcycles to perform.");
        struct arg_lit *no_refinement_in_last_iteration= arg_lit0(NULL, "no_refinement_in_last_iteration","No local search during last v-cycle.");
        struct arg_lit *converter_evaluate             = arg_lit0(NULL, "evaluate","Enable this tag the partition to be evaluated.");
        struct arg_lit *save_partition		       = arg_lit0(NULL, "save_partition","Enable this tag if you want to store the partition to disk.");
        struct arg_lit *save_partition_binary	       = arg_lit0(NULL, "save_partition_binary","Enable this tag if you want to store the partition to disk in a binary format.");
        struct arg_lit *vertex_degree_weights          = arg_lit0(NULL, "vertex_degree_weights","Use 1+deg(v) as vertex weights.");
        struct arg_rex *node_ordering                  = arg_rex0(NULL, "node_ordering", "^(random|degree|leastghostnodesfirst_degree|degree_leastghostnodesfirst)$", "VARIANT", REG_EXTENDED, "Type of node ordering to use for the clustering algorithm. (Default: degree) [random|degree|leastghostnodesfirst_degree|degree_leastghostnodesfirst]." );
        struct arg_rex *preconfiguration               = arg_rex1(NULL, "preconfiguration", "^(ecosocial|fastsocial|ultrafastsocial|ecomesh|fastmesh|ultrafastmesh)$", "VARIANT", REG_EXTENDED, "Use a preconfiguration. (Default: fast) [ecosocial|fastsocial|ultrafastsocial|ecomesh|fastmesh|ultrafastmesh]." );
        struct arg_dbl *ht_fill_factor                 = arg_dbl0(NULL, "ht_fill_factor", NULL, "");
        struct arg_int *n                              = arg_int0(NULL, "n", NULL, "");
        struct arg_end *end                            = arg_end(100);

        //integrated mapping
        //

        //struct arg_lit *integrated_mapping                   = arg_lit0(NULL, "integrated_mapping", "Enable integrated mapping algorithms to map quotient graph onto processor graph defined by hierarchy and distance options. (Default: disabled)");
        struct arg_str *hierarchy_parameter_string           = arg_str0(NULL, "hierarchy_parameter_string", NULL, "Specify as 4:8:8 for 4 cores per PE, 8 PEs per rack, ... and so forth; in total 4x8x8=256 PEs.");
        struct arg_str *distance_parameter_string            = arg_str0(NULL, "distance_parameter_string", NULL, "Specify as 1:10:100 if cores on the same chip have distance 1, PEs in the same rack have distance 10, ... and so forth.");
        struct arg_lit *only_boundary                        = arg_lit0(NULL, "only_boundary", "when refinement, move boundary vertices only" );

        // Define argtable.
        void* argtable[] = {
#ifdef PARALLEL_LABEL_COMPRESSION
	  help, filename, filename_output, user_seed, k, inbalance, preconfiguration, vertex_degree_weights,
                save_partition, save_partition_binary, hierarchy_parameter_string, distance_parameter_string,
                only_boundary, num_vcycles, label_iterations_refinement, label_iterations_coarsening, stop_factor,
                no_refinement_in_last_iteration,
#elif defined TOOLBOX 
                help, filename, k_opt, input_partition_filename, save_partition, save_partition_binary, converter_evaluate,
#endif 
                 end
        };

        // Parse arguments.
        int nerrors = arg_parse(argn, argv, argtable);

        // Catch case that help was requested.
        if(help->count > 0) {
                int rank;
                MPI_Comm_rank( MPI_COMM_WORLD, &rank);

                if( rank == ROOT ) {
                        printf("Usage: %s", progname);
                        arg_print_syntax(stdout, argtable, "\n");
                        arg_print_glossary(stdout, argtable,"  %-40s %s\n");
                        printf("This is the experimental parallel partitioner program.\n");
                        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
                }
                return 1;
        }


        if(nerrors > 0) {
                int rank;
                MPI_Comm_rank( MPI_COMM_WORLD, &rank);
                if( rank == ROOT ) {
                        arg_print_errors(stderr, end, progname);
                        printf("Try '%s --help' for more information.\n",progname);
                        arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
                }
                return 1; 
        }

        configuration cfg;
        cfg.standard(partition_config);

        if(k->count > 0) {
                partition_config.k = k->ival[0];
        }

        if(k_opt->count > 0) {
                partition_config.k = k_opt->ival[0];
        }


#ifdef PARALLEL_LABEL_COMPRESSION
        if(filename->count > 0) {
                graph_filename = filename->sval[0];
        } else {
                if(partition_config.generate_rgg == false && partition_config.generate_ba == false) {
                        printf("You must specify a filename or enable the graph generator tag.\n");
                        return 1;
                }
        }

#else 
        if(filename->count > 0) {
                graph_filename = filename->sval[0];
        }
#endif
        if(input_partition_filename->count > 0) {
                partition_config.input_partition_filename = input_partition_filename->sval[0];
        }

        if(preconfiguration->count > 0) {
                if (strcmp("ecosocial", preconfiguration->sval[0]) == 0) {
                        cfg.eco(partition_config);
                } else if (strcmp("fastsocial", preconfiguration->sval[0]) == 0) {
                        cfg.fast(partition_config);
                } else if (strcmp("ultrafastsocial", preconfiguration->sval[0]) == 0) {
                        cfg.ultrafast(partition_config);
                } else if (strcmp("ecomesh", preconfiguration->sval[0]) == 0) {
                        cfg.eco(partition_config);
                        partition_config.cluster_coarsening_factor = 20000;
                } else if (strcmp("fastmesh", preconfiguration->sval[0]) == 0) {
                        cfg.fast(partition_config);
                        partition_config.cluster_coarsening_factor = 20000;
                } else if (strcmp("ultrafastmesh", preconfiguration->sval[0]) == 0) {
                        cfg.ultrafast(partition_config);
                        partition_config.cluster_coarsening_factor = 20000;
                } else {
                        fprintf(stderr, "Invalid preconfconfiguration variant: \"%s\"\n", preconfiguration->sval[0]);
                        exit(0);
                }
        }

        if(vertex_degree_weights->count > 0) {
                partition_config.vertex_degree_weights = true;
        }

	if(converter_evaluate->count > 0) {
		partition_config.converter_evaluate = true;
	}

	if(save_partition->count > 0) {
		partition_config.save_partition = true;
	}

	if(filename_output->count > 0) {
	        partition_config.filename_output = filename_output->sval[0];
        }

	if(save_partition_binary->count > 0) {
		partition_config.save_partition_binary = true;
	}

        if(n->count > 0) {
                partition_config.n = pow(10,n->ival[0]);
        }

        if(inbalance->count > 0) {
                partition_config.epsilon = inbalance->ival[0];
        }

        if(num_tries->count > 0) {
                partition_config.num_tries = num_tries->ival[0];
        }

        if (user_seed->count > 0) {
                partition_config.seed = user_seed->ival[0];
        }

        if (binary_io_window_size->count > 0) {
                partition_config.binary_io_window_size = binary_io_window_size->ival[0];
        }

        if (num_vcycles->count > 0) {
                partition_config.num_vcycles = num_vcycles->ival[0];
        }

        if (comm_rounds->count > 0) {
                partition_config.comm_rounds = comm_rounds->ival[0];
        }

        if (inbalance->count > 0) {
                partition_config.inbalance = inbalance->ival[0];
        }

        if (ht_fill_factor->count > 0) {
                partition_config.ht_fill_factor = ht_fill_factor->dval[0];
        }


        if (evolutionary_time_limit->count > 0) {
                int size;
                MPI_Comm_size( MPI_COMM_WORLD, &size);
                partition_config.evolutionary_time_limit = evolutionary_time_limit->ival[0]/size;
        }

        if( label_iterations_coarsening->count > 0)  {
                partition_config.label_iterations_coarsening = label_iterations_coarsening->ival[0];
        }

        if( label_iterations_refinement->count > 0)  {
                partition_config.label_iterations_refinement = label_iterations_refinement->ival[0];
        }

        if( cluster_coarsening_factor->count > 0)  {
                partition_config.cluster_coarsening_factor = cluster_coarsening_factor->dval[0];
        }


        if(stop_factor->count > 0) {
                partition_config.stop_factor = stop_factor->ival[0];
        }

        if(no_refinement_in_last_iteration->count > 0) {
                partition_config.no_refinement_in_last_iteration = true;
        }else{
                partition_config.no_refinement_in_last_iteration = false;
        }

        if(initial_partitioning_algorithm->count > 0) {
                if(strcmp("kaffpaEstrong", initial_partitioning_algorithm->sval[0]) == 0) {
                        partition_config.initial_partitioning_algorithm = KAFFPAESTRONG;
                } else if (strcmp("kaffpaEeco",initial_partitioning_algorithm->sval[0]) == 0) {
                        partition_config.initial_partitioning_algorithm = KAFFPAEECO;
                } else if (strcmp("kaffpaEfast", initial_partitioning_algorithm->sval[0]) == 0) {
                        partition_config.initial_partitioning_algorithm = KAFFPAEFAST;
                } else if (strcmp("fastsocial", initial_partitioning_algorithm->sval[0]) == 0) {
                        partition_config.initial_partitioning_algorithm = KAFFPAEFASTSNW;
                } else if (strcmp("ecosocial", initial_partitioning_algorithm->sval[0]) == 0) {
                        partition_config.initial_partitioning_algorithm = KAFFPAEECOSNW;
                } else if (strcmp("strongsocial", initial_partitioning_algorithm->sval[0]) == 0) {
                        partition_config.initial_partitioning_algorithm = KAFFPAESTRONGSNW;
                } else if (strcmp("random", initial_partitioning_algorithm->sval[0]) == 0) {
                        partition_config.initial_partitioning_algorithm = RANDOMIP;
                } else {
                        fprintf(stderr, "Invalid initial partitioning algorithm: \"%s\"\n", initial_partitioning_algorithm->sval[0]);
                        exit(0);
                }
        }

        if(node_ordering->count > 0) {
                if(strcmp("random", node_ordering->sval[0]) == 0) {
                        partition_config.node_ordering = RANDOM_NODEORDERING;
                } else if (strcmp("degree", node_ordering->sval[0]) == 0) {
                        partition_config.node_ordering = DEGREE_NODEORDERING;
                } else if (strcmp("leastghostnodesfirst_degree", node_ordering->sval[0]) == 0) {
                        partition_config.node_ordering = LEASTGHOSTNODESFIRST_DEGREE_NODEODERING;
                } else if (strcmp("degree_leastghostnodesfirst", node_ordering->sval[0]) == 0) {
                        partition_config.node_ordering = DEGREE_LEASTGHOSTNODESFIRST_NODEODERING;
                } else {
                        fprintf(stderr, "Invalid node ordering variant: \"%s\"\n", node_ordering->sval[0]);
                        exit(0);
                }
        }

        //
        // integrated mapping
        //
        partition_config.integrated_mapping = false;

        //check and store the hierarchy levels in partition_config.group_sizes
        if(hierarchy_parameter_string->count) {
                std::istringstream f(hierarchy_parameter_string->sval[0]);
                std::string s;    
                partition_config.group_sizes.clear();
                while (getline(f, s, ':')) {
                        partition_config.group_sizes.push_back(stoi(s));
                }       

                PartitionID old_k = partition_config.k;
                partition_config.k = 1; // recompute k 
                for( unsigned int i = 0; i < partition_config.group_sizes.size(); i++) {
                        partition_config.k *= partition_config.group_sizes[i];
                }
                if( old_k != partition_config.k ) {
                        std::cout <<  "number of processors defined through specified hierarchy does not match k!"  << std::endl;
                        std::cout <<  "please specify k as " << partition_config.k  << std::endl;
                        exit(0);
                }
        }

        //store the PE tree distances in partition_config.distances
        if(distance_parameter_string->count) {
                std::istringstream f(distance_parameter_string->sval[0]);
                std::string s;    
                partition_config.distances.clear();
                while (getline(f, s, ':')) {
                        partition_config.distances.push_back(stoi(s));
                }       
        }

        if( partition_config.distances.size()!=partition_config.group_sizes.size() ){
            std::cout << "ERROR: distances and hierarchy for the processor tree should have the same size " << std::endl;
            std::cout<< partition_config.distances.size() << " vs " << partition_config.group_sizes.size() <<std::endl;
            exit(0);
        }else{
            if( partition_config.distances.size()>0 ){
                partition_config.integrated_mapping = true;
            }
        }

        if( only_boundary->count ){
            partition_config.only_boundary = true;
        }

//next lines appear in main() at the SEA_mapping code; not sure if (and why) we need them
//see for example SEA_mapping/app/fastmesh.cpp
/*
        
        std::vector< NodeID > *perm_rank = NULL;
        if (partition_config.enable_mapping || partition_config.integrated_mapping) {
                perm_rank = new std::vector< NodeID >(partition_config.k);
                for( unsigned i = 0; i < perm_rank->size(); i++) {
                        (*perm_rank)[i] = i;
                }
                partition_config.perm_rank = perm_rank;
        }
*/

//next lines are also taken from SEA_mapping/app/fastmesh.cpp

/*
        matrix* D=NULL;
        
        if ( partition_config.integrated_mapping ){
                //commenting this out so I do not need to add construction_algorithm to partition config
                //bool power_of_two = (partition_config.k & (partition_config.k-1)) == 0;
                // if (power_of_two  // && !partition_config.enable_mapping && !partition_config.multisection 
                // ) {
                //         partition_config.construction_algorithm = MAP_CONST_IDENTITY;
                // }
                

                if (partition_config.use_bin_id) {
                        D = new online_precalc_matrix(partition_config.k, partition_config.k);
                        D->setPartitionConfig(partition_config);
                } else if (partition_config.use_compact_bin_id) {
                        D = new online_binary_matrix(partition_config.k, partition_config.k);
                        D->setPartitionConfig(partition_config);
                } else if (partition_config.full_matrix) {
                        D = new full_matrix(partition_config.k, partition_config.k);
                        D->setPartitionConfig(partition_config);
                } else  
                        if( partition_config.distance_construction_algorithm != DIST_CONST_HIERARCHY_ONLINE) {
                                D = new normal_matrix(partition_config.k, partition_config.k);
                        } else {
                                D = new online_distance_matrix(partition_config.k, partition_config.k);
                                D->setPartitionConfig(partition_config);
                        }
                partition_config.D = D;
        }
*/

        return 0;
}

#endif /* end of include guard: PARSE_PARAMETERS_GPJMGSM8 */
