/******************************************************************************
 * graph_partitioner.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PARTITION_OL9XTLU4
#define PARTITION_OL9XTLU4

#include "coarsening/coarsening.h"
#include "coarsening/stop_rules/stop_rules.h"
#include "data_structure/graph_access.h"
#include "partition_config.h"
#include "uncoarsening/refinement/refinement.h"

class graph_partitioner {
public:
        graph_partitioner();
        virtual ~graph_partitioner();

        void perform_partitioning(PartitionConfig & graph_partitioner_config, graph_access & G);
        void perform_recursive_partitioning(PartitionConfig & graph_partitioner_config, graph_access & G);

private:
        void perform_recursive_partitioning_internal(PartitionConfig & graph_partitioner_config, 
                                                     graph_access & G, 
                                                     PartitionID lb, PartitionID ub);
        void single_run( PartitionConfig & config, graph_access & G);

        unsigned m_global_k;
	int m_global_upper_bound;
        int m_rnd_bal;
};

#endif /* end of include guard: PARTITION_OL9XTLU4 */
