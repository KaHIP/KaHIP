/******************************************************************************
 * misc.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef MISC_C6QUUWLI
#define MISC_C6QUUWLI

#include "data_structure/graph_access.h"
#include "partition_config.h"

class misc {
public:
        misc();
        virtual ~misc();

        void balance_singletons(const PartitionConfig & config, graph_access & G);
};


#endif /* end of include guard: MISC_C6QUUWLI */
