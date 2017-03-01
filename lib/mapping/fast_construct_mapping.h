/******************************************************************************
 * fast_construct_mapping.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Graph Partitioning 
 ******************************************************************************
 * Copyright (C) 2017 Christian Schulz 
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

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
