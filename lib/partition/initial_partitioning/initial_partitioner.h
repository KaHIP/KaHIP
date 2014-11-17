/******************************************************************************
 * initial_partitioner.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
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

#ifndef INITIAL_PARTITIONER_TJKC6RWY
#define INITIAL_PARTITIONER_TJKC6RWY

#include "partition_config.h"
#include "data_structure/graph_access.h"

class initial_partitioner {
        public:
                initial_partitioner( );
                virtual ~initial_partitioner();

                virtual void initial_partition( const PartitionConfig & config, const unsigned int seed,  
                                graph_access & G, 
                                int* xadj,
                                int* adjncy, 
                                int* vwgt, 
                                int* adjwgt,
                                int* partition_map) = 0; 

                virtual void initial_partition(const PartitionConfig & config, 
                                               const unsigned int seed,  
                                               graph_access & G, 
                                               int* partition_map) = 0;
};


#endif /* end of include guard: INITIAL_PARTITIONER_TJKC6RWY */
