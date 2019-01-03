/******************************************************************************
 * random_initial_partitioning.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef RANDOM_INITIAL_PARTITIONING_FM8LJSI0
#define RANDOM_INITIAL_PARTITIONING_FM8LJSI0

#include <mpi.h>
#include "partition_config.h"

class parallel_graph_access;

class random_initial_partitioning {
public:
        random_initial_partitioning();
        virtual ~random_initial_partitioning();
        
        void perform_partitioning( MPI_Comm communicator, PPartitionConfig & config, parallel_graph_access & G);
};


#endif /* end of include guard: RANDOM_INITIAL_PARTITIONING_FM8LJSI0 */
