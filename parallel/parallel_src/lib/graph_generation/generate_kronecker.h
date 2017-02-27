/******************************************************************************
 * generate_kronecker.h 
 *
 * Source of the Parallel Partitioning Program
 ******************************************************************************
 * Copyright (C) 2014 Christian Schulz <christian.schulz@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
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

#ifndef GENERATE_KRONECKER_UM0MTHVS
#define GENERATE_KRONECKER_UM0MTHVS

#include "partition_config.h"
#include "data_structure/parallel_graph_access.h"

class generate_kronecker {
public:
        generate_kronecker();
        virtual ~generate_kronecker();

        // G is the output parameter
        void generate_kronecker_graph( PPartitionConfig & config, parallel_graph_access & G);
};




#endif /* end of include guard: GENERATE_KRONECKER_UM0MTHVS */
