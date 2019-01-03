/******************************************************************************
 * distributed_evolutionary_partitioning.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef DISTRIBUTED_EVOLUTIONARY_PARTITIONING_OJ2RIKR7
#define DISTRIBUTED_EVOLUTIONARY_PARTITIONING_OJ2RIKR7

#include "partition_config.h"
#include "data_structure/parallel_graph_access.h"

class distributed_evolutionary_partitioning {
public:
        distributed_evolutionary_partitioning();
        virtual ~distributed_evolutionary_partitioning();

        void perform_partitioning( MPI_Comm communicator, PPartitionConfig & config, parallel_graph_access & G);
};

#endif /* end of include guard: DISTRIBUTED_EVOLUTIONARY_PARTITIONING_OJ2RIKR7 */
