/******************************************************************************
 * parallel_mh_async.h 
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

#ifndef PARALLEL_MH_ASYNC_HF106Y0G
#define PARALLEL_MH_ASYNC_HF106Y0G

#include "data_structure/graph_access.h"
#include "partition_config.h"
#include "population.h"
#include "timer.h"

class parallel_mh_async {
public:
        parallel_mh_async();
        virtual ~parallel_mh_async();

        void perform_partitioning(const PartitionConfig & graph_partitioner_config, graph_access & G);
        void initialize(PartitionConfig & graph_partitioner_config, graph_access & G);
        EdgeWeight perform_local_partitioning(PartitionConfig & graph_partitioner_config, graph_access & G);
        EdgeWeight collect_best_partitioning(graph_access & G);
        void perform_cycle_partitioning(PartitionConfig & graph_partitioner_config, graph_access & G);

private:
        //misc
        const unsigned MASTER;
        timer    m_t;
        int      m_rank;
        int      m_size;
        double   m_time_limit;
        bool     m_termination;
        unsigned m_rounds;

        //the best cut found so far
        PartitionID* m_best_global_map;
        int          m_best_global_objective;
        int          m_best_cycle_objective;

        //island
        population* m_island;
};


#endif /* end of include guard: PARALLEL_MH_ASYNC_HF106Y0G */
