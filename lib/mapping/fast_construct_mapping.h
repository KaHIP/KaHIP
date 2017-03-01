//
// Author: Christian Schulz <christian.schulz@kit.edu>
// 

#ifndef FAST_CONSTRUCT_MAPPING_1MEOBVNJ
#define FAST_CONSTRUCT_MAPPING_1MEOBVNJ

#include "data_structure/graph_access.h"
#include "data_structure/matrix/matrix.h"
#include "partition_config.h"
#include "tools/quality_metrics.h"

class fast_construct_mapping {
public:
        fast_construct_mapping();
        virtual ~fast_construct_mapping();

        void construct_initial_mapping_bottomup( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank);
        void construct_initial_mapping_topdown( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank);

private:

        void construct_initial_mapping_topdown_internal( PartitionConfig & config, graph_access & C,  
                                                         std::vector< int >  group_sizes, int start_id, 
                                                         std::vector< NodeID > & map_to_original, 
                                                         std::vector< NodeID > & perm_rank);

        void construct_initial_mapping_bottomup_internal( PartitionConfig & config, graph_access & C, matrix & D, int idx,  std::vector< NodeID > & perm_rank);


        void partition_C_perfectly_balanced( PartitionConfig & config, graph_access & C, PartitionID blocks);

        int m_tmp_num_nodes;

};

#endif /* end of include guard: FAST_CONSTRUCT_MAPPING_1MEOBVNJ */
