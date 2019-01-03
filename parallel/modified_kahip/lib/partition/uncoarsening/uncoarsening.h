/******************************************************************************
 * uncoarsening.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef UNCOARSENING_XSN847F2
#define UNCOARSENING_XSN847F2

#include "data_structure/graph_hierarchy.h"
#include "partition_config.h"

class uncoarsening {
public:
        uncoarsening( );
        virtual ~uncoarsening();
        
        int perform_uncoarsening(const PartitionConfig & config, graph_hierarchy & hierarchy);
};


#endif /* end of include guard: UNCOARSENING_XSN847F2 */
