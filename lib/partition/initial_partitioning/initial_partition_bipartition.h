/******************************************************************************
 * initial_partition_bipartition.h 
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

#ifndef INITIAL_PARTITION_BIPARTITION_HMA7329W
#define INITIAL_PARTITION_BIPARTITION_HMA7329W

#include "initial_partitioner.h"

class initial_partition_bipartition : public initial_partitioner {
public:
        initial_partition_bipartition();
        virtual ~initial_partition_bipartition();

        void initial_partition( const PartitionConfig & config, const unsigned int seed,  graph_access & G, int* partition_map); 

        void initial_partition( const PartitionConfig & config, const unsigned int seed,  
                                graph_access & G, 
                                int* xadj,
                                int* adjncy, 
                                int* vwgt, 
                                int* adjwgt,
                                int* partition_map); 

};


#endif /* end of include guard: INITIAL_PARTITION_BIPARTITION_HMA7329W */
