/******************************************************************************
 * uncoarsening.h 
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

#ifndef UNCOARSENING_XSN847F2
#define UNCOARSENING_XSN847F2

#include "data_structure/graph_hierarchy.h"
#include "partition_config.h"

class uncoarsening {
public:
        uncoarsening( );
        virtual ~uncoarsening();
        
        int perform_uncoarsening(const PartitionConfig & config, graph_hierarchy & hierarchy);
        int perform_uncoarsening_cut(const PartitionConfig & config, graph_hierarchy & hierarchy);
        int perform_uncoarsening_nodeseparator(const PartitionConfig & config, graph_hierarchy & hierarchy);
        int perform_uncoarsening_nodeseparator_fast(const PartitionConfig & config, graph_hierarchy & hierarchy);
};


#endif /* end of include guard: UNCOARSENING_XSN847F2 */
