/******************************************************************************
 * initial_partitioning.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "distributed_evolutionary_partitioning.h"
#include "initial_partitioning.h"
#include "random_initial_partitioning.h"

initial_partitioning_algorithm::initial_partitioning_algorithm() {
                
}

initial_partitioning_algorithm::~initial_partitioning_algorithm() {
                
}


void initial_partitioning_algorithm::perform_partitioning( MPI_Comm communicator, PPartitionConfig & config, parallel_graph_access & Q) {
        if( config.initial_partitioning_algorithm == RANDOMIP) {
                random_initial_partitioning dist_rpart;
                dist_rpart.perform_partitioning( communicator, config, Q );
        } else {
                distributed_evolutionary_partitioning dist_epart;
                dist_epart.perform_partitioning( communicator, config, Q);
        }
}

