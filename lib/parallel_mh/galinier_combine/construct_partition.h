/******************************************************************************
 * construct_partition.h 
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

#ifndef CONSTRUCT_PARTITION_E86DQF5S
#define CONSTRUCT_PARTITION_E86DQF5S

#include "data_structure/graph_access.h"
#include "parallel_mh/population.h"
#include "partition_config.h"

class construct_partition {
public:
        construct_partition();
        virtual ~construct_partition();

        void construct_starting_from_partition( PartitionConfig & config, graph_access & G);
        void createIndividuum( PartitionConfig & config, graph_access & G, 
                               Individuum & ind, bool output); 
};


#endif /* end of include guard: CONSTRUCT_PARTITION_E86DQF5S */
