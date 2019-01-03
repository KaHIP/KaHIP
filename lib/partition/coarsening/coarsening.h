/******************************************************************************
 * coarsening.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef COARSENING_UU97ZBTR
#define COARSENING_UU97ZBTR

#include "data_structure/graph_access.h"
#include "data_structure/graph_hierarchy.h"
#include "partition_config.h"

class coarsening {
public:
        coarsening ();
        virtual ~coarsening ();

        void perform_coarsening(const PartitionConfig & config, graph_access & G, graph_hierarchy & hierarchy);
};

#endif /* end of include guard: COARSENING_UU97ZBTR */
