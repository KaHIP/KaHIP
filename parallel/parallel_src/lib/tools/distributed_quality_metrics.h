/******************************************************************************
 * distributed_quality_metrics.h
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

#ifndef DISTRIBUTED_QUALITY_METRICS_UAVSEXBT
#define DISTRIBUTED_QUALITY_METRICS_UAVSEXBT

#include "definitions.h"
#include "data_structure/parallel_graph_access.h"
#include "partition_config.h"

class distributed_quality_metrics {
public:
        distributed_quality_metrics();
        virtual ~distributed_quality_metrics();

        EdgeWeight local_edge_cut( parallel_graph_access & G, int * partition_map, MPI_Comm communicator );
        EdgeWeight edge_cut( parallel_graph_access & G, MPI_Comm communicator );
        EdgeWeight edge_cut_second( parallel_graph_access & G, MPI_Comm communicator  );
        NodeWeight local_max_block_weight( PPartitionConfig & config, parallel_graph_access & G, int * partition_map, MPI_Comm communicator  );
        double balance( PPartitionConfig & config, parallel_graph_access & G, MPI_Comm communicator  );
        double balance_load( PPartitionConfig & config, parallel_graph_access & G, MPI_Comm communicator  );
        double balance_load_dist( PPartitionConfig & config, parallel_graph_access & G, MPI_Comm communicator  );
        double balance_second( PPartitionConfig & config, parallel_graph_access & G, MPI_Comm communicator  );
	EdgeWeight comm_vol( PPartitionConfig & config, parallel_graph_access & G, MPI_Comm communicator  );
	EdgeWeight comm_vol_dist( parallel_graph_access & G, MPI_Comm communicator );
};


#endif /* end of include guard: DISTRIBUTED_QUALITY_METRICS_UAVSEXBT */
