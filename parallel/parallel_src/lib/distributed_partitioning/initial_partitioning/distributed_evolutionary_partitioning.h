/******************************************************************************
 * distributed_evolutionary_partitioning.h
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

#ifndef DISTRIBUTED_EVOLUTIONARY_PARTITIONING_OJ2RIKR7
#define DISTRIBUTED_EVOLUTIONARY_PARTITIONING_OJ2RIKR7

#include "partition_config.h"
#include "data_structure/parallel_graph_access.h"

class distributed_evolutionary_partitioning {
public:
        distributed_evolutionary_partitioning();
        virtual ~distributed_evolutionary_partitioning();

        void perform_partitioning( MPI_Comm communicator, PPartitionConfig & config, parallel_graph_access & G);
};

#endif /* end of include guard: DISTRIBUTED_EVOLUTIONARY_PARTITIONING_OJ2RIKR7 */
