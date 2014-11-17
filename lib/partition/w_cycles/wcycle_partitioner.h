/******************************************************************************
 * wcycle_partitioner.h 
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

#ifndef WCYCLE_PARTITIONER_EPNDQMK
#define WCYCLE_PARTITIONER_EPNDQMK

#include "coarsening/coarsening.h"
#include "coarsening/stop_rules/stop_rules.h"
#include "data_structure/graph_access.h"
#include "partition_config.h"
#include "uncoarsening/refinement/refinement.h"

class wcycle_partitioner {
        public:
                wcycle_partitioner( ) : m_level(0) {};
                virtual ~wcycle_partitioner() {};
                int perform_partitioning( const PartitionConfig & config, 
                                          graph_access & G); 

        private:
                int perform_partitioning_recursive( PartitionConfig & partition_config, 
                                                    graph_access & G, 
                                                    complete_boundary ** c_boundary); 

                unsigned   m_level;
                unsigned   m_deepest_level;
                stop_rule* m_coarsening_stop_rule;

                std::unordered_map<unsigned, bool> m_have_been_level_down;
};

#endif /* end of include guard: WCYCLE_PARTITIONER_EPNDQMK */
