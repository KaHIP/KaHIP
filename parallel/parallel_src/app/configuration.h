/******************************************************************************
 * configuration.h
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


#ifndef CONFIGURATION_3APG5V7ZA
#define CONFIGURATION_3APG5V7ZA

#include "partition_config.h"

class configuration {
        public:
                configuration() {} ;
                virtual ~configuration() {};

                void standard( PPartitionConfig & config );
                void ultrafast( PPartitionConfig & config );
                void fast( PPartitionConfig & config );
                void eco( PPartitionConfig & config );
                void strong( PPartitionConfig & config );
};

inline void configuration::ultrafast( PPartitionConfig & partition_config ) {
        partition_config.initial_partitioning_algorithm  = KAFFPAEULTRAFASTSNW;
        partition_config.no_refinement_in_last_iteration = true;
        partition_config.stop_factor                     = 18000;
        partition_config.num_vcycles                     = 1;
}


inline void configuration::fast( PPartitionConfig & partition_config ) {
        partition_config.initial_partitioning_algorithm  = KAFFPAEULTRAFASTSNW;
        partition_config.no_refinement_in_last_iteration = true;
        partition_config.stop_factor                     = 18000;
}

inline void configuration::eco( PPartitionConfig & partition_config ) {
        partition_config.initial_partitioning_algorithm  = KAFFPAEFASTSNW;
        partition_config.no_refinement_in_last_iteration = true;
        partition_config.stop_factor                     = 18000;
        int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
        partition_config.evolutionary_time_limit         = 2048/size;
        partition_config.eco                             = true;
        partition_config.num_vcycles                     = 6;
}

inline void configuration::strong( PPartitionConfig & partition_config ) {
        partition_config.initial_partitioning_algorithm = KAFFPAESTRONGSNW;

}
inline void configuration::standard( PPartitionConfig & partition_config ) {
        partition_config.seed                                   = 0;
        partition_config.k                                      = 2;
        partition_config.inbalance                              = 3;
        partition_config.epsilon                                = 3;
        partition_config.time_limit 				= 0; 
        partition_config.evolutionary_time_limit 	        = 0; 
        partition_config.log_num_verts                          = 16;
        partition_config.edge_factor                            = 16;
        partition_config.generate_rgg                           = false; 
        partition_config.generate_ba                            = false; 
        partition_config.comm_rounds                            = 128; 
        partition_config.label_iterations                       = 4;
        partition_config.label_iterations_coarsening            = 3;
        partition_config.label_iterations_refinement            = 6;
        partition_config.cluster_coarsening_factor              = 14;
        partition_config.initial_partitioning_algorithm         = KAFFPAEFASTSNW;
        partition_config.stop_factor                            = 14000;
        partition_config.vcycle                                 = false;
        partition_config.num_vcycles                            = 2;
        partition_config.num_tries                              = 10;
        partition_config.node_ordering                          = DEGREE_NODEORDERING;
        partition_config.no_refinement_in_last_iteration        = false;
        partition_config.ht_fill_factor                         = 1.6;
        partition_config.eco                                    = false;
	partition_config.binary_io_window_size                  = 64;
        partition_config.barabasi_albert_mindegree              = 5;
        partition_config.compute_degree_sequence_ba             = true;
        partition_config.compute_degree_sequence_k_first        = false;
        partition_config.k_deg                                  = 1000000;
        partition_config.kronecker_internal_only                = false;
        partition_config.generate_ba_32bit                      = false;
        partition_config.n                                      = 0;
	partition_config.save_partition 			= false;
	partition_config.save_partition_binary 			= false;
        partition_config.vertex_degree_weights                  = false;
        partition_config.converter_evaluate                     = false;
}

#endif /* end of include guard: CONFIGURATION_3APG5V7Z */
