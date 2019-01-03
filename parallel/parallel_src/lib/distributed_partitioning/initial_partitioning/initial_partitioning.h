/******************************************************************************
 * initial_partitioning.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef INITIAL_PARTITIONING_SFMCJN2U
#define INITIAL_PARTITIONING_SFMCJN2U

#include "partition_config.h"
#include "data_structure/parallel_graph_access.h"

class initial_partitioning_algorithm {
public:
        initial_partitioning_algorithm();
        virtual ~initial_partitioning_algorithm();

        void perform_partitioning( MPI_Comm communicator, PPartitionConfig & config, parallel_graph_access & G);
};


#endif /* end of include guard: INITIAL_PARTITIONING_SFMCJN2U */
