/******************************************************************************
 * graph_partitioner.h
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
