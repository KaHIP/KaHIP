/******************************************************************************
 * graph_partition_assertions.h 
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

#ifndef GRAPH_PARTITION_ASSERTIONS_609QZZDM
#define GRAPH_PARTITION_ASSERTIONS_609QZZDM

#include "data_structure/graph_access.h"
#include "partition_config.h"

class graph_partition_assertions {
        public:
                graph_partition_assertions( ) {};
                virtual ~graph_partition_assertions() {};

                static bool assert_graph_has_kway_partition(const PartitionConfig & config, graph_access & G) {
                        bool* allpartsthere = new bool[config.k];
                        for(unsigned int i = 0; i < config.k; i++) {
                                allpartsthere[i] = false;
                        }

                        forall_nodes(G, n) {
                                allpartsthere[G.getPartitionIndex(n)] = true; 
                        } endfor

                        for(unsigned int i = 0; i < config.k; i++) {
                                ASSERT_TRUE(allpartsthere[i]);
                        }

                        delete[] allpartsthere;
                        return true;
                };

};


#endif /* end of include guard: GRAPH_PARTITION_ASSERTIONS_609QZZDM */
