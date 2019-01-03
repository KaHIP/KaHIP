/******************************************************************************
 * initial_partition_bipartition.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef INITIAL_PARTITION_BIPARTITION_HMA7329W
#define INITIAL_PARTITION_BIPARTITION_HMA7329W

#include "initial_partitioner.h"

class initial_partition_bipartition : public initial_partitioner {
public:
        initial_partition_bipartition();
        virtual ~initial_partition_bipartition();

        void initial_partition( const PartitionConfig & config, const unsigned int seed,  graph_access & G, int* partition_map); 

        void initial_partition( const PartitionConfig & config, const unsigned int seed,  
                                graph_access & G, 
                                int* xadj,
                                int* adjncy, 
                                int* vwgt, 
                                int* adjwgt,
                                int* partition_map); 

};


#endif /* end of include guard: INITIAL_PARTITION_BIPARTITION_HMA7329W */
