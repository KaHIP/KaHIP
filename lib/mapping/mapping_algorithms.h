/******************************************************************************
 * mapping_algorithms.h
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

#ifndef MAPPING_ALGORITHMS_W4I4JZHS
#define MAPPING_ALGORITHMS_W4I4JZHS

#include "data_structure/graph_access.h"
#include "data_structure/matrix/normal_matrix.h"
#include "partition_config.h"
#include "tools/quality_metrics.h"
#include "tools/timer.h"

class mapping_algorithms {
public:
        mapping_algorithms();
        virtual ~mapping_algorithms();

        void construct_a_mapping( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank);
        void graph_to_matrix( graph_access & C, matrix & C_bar);

private:
        quality_metrics qm; 
        timer t;
};


#endif /* end of include guard: MAPPING_ALGORITHMS_W4I4JZHS */
