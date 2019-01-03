/******************************************************************************
 * initial_partitioner.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef INITIAL_PARTITIONER_TJKC6RWY
#define INITIAL_PARTITIONER_TJKC6RWY

#include "partition_config.h"
#include "data_structure/graph_access.h"

class initial_partitioner {
        public:
                initial_partitioner( );
                virtual ~initial_partitioner();

                virtual void initial_partition( const PartitionConfig & config, const unsigned int seed,  
                                graph_access & G, 
                                int* xadj,
                                int* adjncy, 
                                int* vwgt, 
                                int* adjwgt,
                                int* partition_map) = 0; 

                virtual void initial_partition(const PartitionConfig & config, 
                                               const unsigned int seed,  
                                               graph_access & G, 
                                               int* partition_map) = 0;
};


#endif /* end of include guard: INITIAL_PARTITIONER_TJKC6RWY */
