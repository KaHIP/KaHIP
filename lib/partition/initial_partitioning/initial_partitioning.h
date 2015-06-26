/******************************************************************************
 * initial_partitioning.h 
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

#ifndef INITIAL_PARTITIONING_D7VA0XO9
#define INITIAL_PARTITIONING_D7VA0XO9

#include "data_structure/graph_hierarchy.h"
#include "partition_config.h"

class initial_partitioning {
public:
        initial_partitioning( );
        virtual ~initial_partitioning();
        void perform_initial_partitioning(const PartitionConfig & config, graph_hierarchy & hierarchy);
        void perform_initial_partitioning(const PartitionConfig & config, graph_access &  G);
        void perform_initial_partitioning_separator(const PartitionConfig & config, graph_access &  G);
};


#endif /* end of include guard: INITIAL_PARTITIONING_D7VA0XO9 */
